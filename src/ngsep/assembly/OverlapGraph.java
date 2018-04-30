/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.assembly;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import ngsep.alignments.ReadAlignment;
import ngsep.math.Distribution;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.KmersCounter;
import ngsep.sequences.LimitedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.sequences.io.FastqFileReader;

/**
 * @author Jorge Duitama
 *
 */
public class OverlapGraph {
	private List<LimitedSequence> sequences = new ArrayList<>();
	private FMIndex index = new FMIndex();
	static final int SEARCH_KMER_LENGTH = 15;
	 
	private Map<Integer,List<ReadOverlap>> overlapsForward = new HashMap<>();
	private Map<Integer,List<ReadOverlap>> overlapsBackward = new HashMap<>();
	private Map<Integer,ReadOverlap> embeddedOverlaps = new HashMap<>();
	
	public void loadFasta (String filename) throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList seqsQl = handler.loadSequences(filename);
		for(QualifiedSequence seq:seqsQl) {
			DNAMaskedSequence characters = (DNAMaskedSequence)seq.getCharacters();
			addSequence(characters);
		}
	}
	
	public void loadFastq(String filename) throws IOException {
		try (FastqFileReader reader = new FastqFileReader(filename)) {
			Iterator<RawRead> it = reader.iterator();
			while(it.hasNext()) {
				RawRead read = it.next();
				DNAMaskedSequence characters = (DNAMaskedSequence)read.getCharacters();
				addSequence(characters);
			}
		}
	}
	
	public void addSequence (LimitedSequence sequence) {
		sequences.add(sequence);
	}
	
	public void buildFMIndex() {
		index = new FMIndex();
		index.loadUnnamedSequences("", sequences);
	}
	
	public void findOverlaps () {
		for(int s=0;s<sequences.size() && s<1000;s++) {
			if(embeddedOverlaps.containsKey(s)) continue;
			
			CharSequence sequence = sequences.get(s);
			System.out.println("Processing sequence: "+s+" length: "+sequence.length());
			Map<Integer,List<ReadAlignment>> seqHits = findHitsRead(s, sequence);
			System.out.println("Found "+seqHits.size()+" potentially matching sequences for sequence "+s+" of length: "+sequence.length());
			buildOverlapsFromKmerAlignments (s, seqHits);
			int nF = overlapsForward.get(s)!=null?overlapsForward.get(s).size():0;
			int nB = overlapsBackward.get(s)!=null?overlapsBackward.get(s).size():0;
			System.out.println("Sequence "+s+" overlaps forward: "+nF + " overlaps backward: "+nB);
			
		}
	}

	public Map<Integer,List<ReadAlignment>>  findHitsRead(int idRead, CharSequence read ) {
		Map<Integer,List<ReadAlignment>> seqHits = new HashMap<>();
		CharSequence [] kmers = KmersCounter.extractKmers(read, SEARCH_KMER_LENGTH, true);
		
		for(int i=0;i<kmers.length;i+=SEARCH_KMER_LENGTH) {
			if(kmers[i]== null) continue;
			String kmer = kmers[i].toString();
			List<ReadAlignment> regions = index.search(kmer);
			//System.out.println("Hits kmer "+kmer+": "+regions.size());
			for(ReadAlignment aln:regions) {
				//Use mate start to store the kmer start site producing the hit
				aln.setMateFirst(i);
				int k = Integer.parseInt(aln.getSequenceName().substring(1));
				if(k>idRead) {
					//String seqOut = instance.sequences.get(k).subSequence(region.getFirst(), region.getLast()+1).toString();
					//if(!seqOut.equals(kmer)) throw new RuntimeException ("Hit error in sequence: "+k+" pos: "+region.getFirst()+" search "+kmer+" found: "+seqOut);
					List<ReadAlignment> seqAlns = seqHits.get(k);
					if(seqAlns==null) {
						seqAlns = new ArrayList<>();
						seqHits.put(k, seqAlns);
					}
					seqAlns.add(aln);
					
				}
			}
		}
		return seqHits;
	}
	
	
	
	private void buildOverlapsFromKmerAlignments(int searchId, Map<Integer, List<ReadAlignment>> seqHits) {
		LimitedSequence searchSequence = sequences.get(searchId);
		//GenomicRegionPositionComparator cmp = GenomicRegionPositionComparator.getInstance();
		for(int k:seqHits.keySet()) {
			List<ReadAlignment> alns = seqHits.get(k);
			//Collections.sort(alns, cmp);
			ReadOverlap next = null;
			for(ReadAlignment aln:alns) {
				if(next == null || !next.addKmerAlignment (k,aln) ) {
					if(next!=null) processOverlap(next,searchSequence.length(),sequences.get(k).length());
					next = new ReadOverlap(searchId, aln.getMateFirst(), aln.getMateFirst()+SEARCH_KMER_LENGTH-1, k, aln.getFirst(), aln.getLast(), aln.isNegativeStrand());
					
				}
			}
			if(next!=null) processOverlap(next,searchSequence.length(),sequences.get(k).length());
		}
		
	}

	
	/**
	 * PRE: length1 > length2
	 * @param overlap
	 * @param length1
	 * @param length2
	 */
	private void processOverlap(ReadOverlap overlap, int length1, int length2) {
		int l1 = overlap.getLast1() - overlap.getFirst1();
		int l2 = overlap.getLast2() - overlap.getFirst2();
		if(l1 < 0.2*length1 && l2 < 0.2*length2) return;
		int dl1 = overlap.getFirst1();
		int dl2 = overlap.getFirst2();
		int dr1 = length1 - overlap.getLast1();
		int dr2 = length2 - overlap.getLast2();
		int t = 10;
		int diffBorder = 100;
		
		if(overlap.isNegativeStrand()) {
			if(dl1+t >= dr2 && dr1+t >= dl2) {
				if(dl2 > diffBorder || dr2> diffBorder) {
					//Not a true overlap
					return;
				}
				if(!embeddedOverlaps.containsKey(overlap.getIndexSequence2())) {
					ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(), overlap.getLast2(), overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), true);
					embeddedOverlaps.put(overlap.getIndexSequence2(), ov2);
				}
				return;
			}
			if(dl1 <= dr2+t && dr1 <= dl2+t) {
				if(dl1 > diffBorder || dr1> diffBorder) {
					//Not a true overlap
					return;
				}
				if(!embeddedOverlaps.containsKey(overlap.getIndexSequence1()))embeddedOverlaps.put(overlap.getIndexSequence1(), overlap);
				return;
			}
			if(dl1 < dr2 && dr1 > dl2 ) {
				//2 left of 1 negative strand
				if(dl1>diffBorder || dl2 > diffBorder) {
					//Not a true overlap
					return;
				}
				addOverlap(overlapsBackward, overlap);
				ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(), overlap.getLast2(), overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), true);
				addOverlap(overlapsForward, ov2);
				return;
			}
			if (dl1 > dr2 && dr1 < dl2) {
				//1 left of 2
				if(dr2>diffBorder || dr1 > diffBorder) {
					//Not a true overlap
					return;
				}
				addOverlap(overlapsForward, overlap);
				ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(), overlap.getLast2(), overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), true);
				addOverlap(overlapsBackward, ov2);
				return;
			}
		} else {
			if(dl1+t >= dl2 && dr1+t >= dr2) {
				if(dl2 > diffBorder || dr2> diffBorder) {
					//Not a true overlap
					return;
				}
				if(!embeddedOverlaps.containsKey(overlap.getIndexSequence2())) {
					ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(), overlap.getLast2(), overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), false);
					embeddedOverlaps.put(overlap.getIndexSequence2(), ov2);
				}
				return;
			}
			if(dl1 <= dl2+t && dr1 <= dr2+t) {
				if(dl1 > diffBorder || dr1> diffBorder) {
					//Not a true overlap
					return;
				}
				if(!embeddedOverlaps.containsKey(overlap.getIndexSequence1()))embeddedOverlaps.put(overlap.getIndexSequence1(), overlap);
				return;
			}
			if(dl1 < dl2 && dr1 > dr2 ) {
				//2 left of 1
				if(dl1>diffBorder || dr2 > diffBorder) {
					//Not a true overlap
					return;
				}
				addOverlap(overlapsBackward, overlap);
				ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(), overlap.getLast2(), overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), false);
				addOverlap(overlapsForward, ov2);
				return;
			} 
			if (dl1 > dl2 && dr1 < dr2) {
				//1 left of 2
				if(dl2>diffBorder || dr1 > diffBorder) {
					//Not a true overlap
					return;
				}
				addOverlap(overlapsForward, overlap);
				ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(), overlap.getLast2(), overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), false);
				addOverlap(overlapsBackward, ov2);
				return;
			}
		}
		//TODO: Print weird overlap
	}

	private void addOverlap(Map<Integer,List<ReadOverlap>> overlaps, ReadOverlap overlap) {
		List<ReadOverlap> overlapsSeq1 = overlaps.get(overlap.getIndexSequence1());
		if(overlapsSeq1==null) {
			overlapsSeq1 = new ArrayList<>();
			overlaps.put(overlap.getIndexSequence1(), overlapsSeq1);
		}
		overlapsSeq1.add(overlap);
	}
	
	public void printOverlapsDistribution(PrintStream out) {
		Distribution distForward = new Distribution(0, 20, 1);
		Distribution distBackward = new Distribution(0, 20, 1);
		int numEmbedded = embeddedOverlaps.size();
		for(int i=0;i<sequences.size();i++) {
			List<ReadOverlap> ovs = overlapsForward.get(i);
			if(ovs == null) distForward.processDatapoint(0);
			else distForward.processDatapoint(ovs.size());
			ovs = overlapsBackward.get(i);
			if(ovs == null) distBackward.processDatapoint(0);
			else distBackward.processDatapoint(ovs.size());
		}
		double [] vF = distForward.getDistribution();
		double [] vB = distBackward.getDistribution();
		System.out.println("Overlaps distribution: ");
		for (int i = 0;i<vF.length;i++) {
			out.println(""+i+"\t"+vF[i]+"\t"+vB[i]);
		}
		out.println("More:\t"+distForward.getOutliers().size()+"\t"+distBackward.getOutliers().size());
		out.println();
		out.println("Embedded alignments: "+numEmbedded);
		
		
	}

	public static void main(String[] args) throws Exception {
		OverlapGraph instance = new OverlapGraph();
		long totalLength = 0;
		long time1 = System.currentTimeMillis();
		/*Random r = new Random();
		for(int i=0;i<numSequences;i++) {
			StringBuilder randomSequence = new StringBuilder();
			int length = r.nextInt(40000)+10000;
			totalLength+=length;
			for(int j=0;j<length;j++) {
				int bpI = r.nextInt(4);
				randomSequence.append(DNASequence.BASES_STRING.charAt(bpI));
			}
			DNASequence dna = new DNASequence(randomSequence);
			instance.addSequence(dna);
			if((i+1)%1000 == 0) System.out.println("Created "+instance.sequences.size()+" random DNA sequences");
		}*/
		String filename = args[0];
		if(filename.endsWith(".fastq") || filename.endsWith(".fastq.gz") ) instance.loadFastq(filename);
		else instance.loadFasta(args[0]);
		long time2 = System.currentTimeMillis();
		System.out.println("Loaded "+instance.sequences.size()+" sequences from "+args[0]+". Total length: "+totalLength +" time: "+ (time2 - time1));
		//System.out.println("First sequence: "+instance.sequences.get(0));
		time1 = time2;
		Collections.sort(instance.sequences, new Comparator<LimitedSequence>() {
			@Override
			public int compare(LimitedSequence l1, LimitedSequence l2) {
				return l1.length()-l2.length();
			}
		});
		Collections.reverse(instance.sequences);
		time2 = System.currentTimeMillis();
		System.out.println("Sorting time: "+ (time2 - time1));
		time1=time2;
		instance.buildFMIndex();
		time2 = System.currentTimeMillis();
		System.out.println("Created FMIndex in time: "+ (time2 - time1));
		time1=time2;
		instance.findOverlaps();
		time2 = System.currentTimeMillis();
		System.out.println("Built overlaps in time: "+ (time2 - time1));
		instance.printOverlapsDistribution (System.out);
		
	}

	
}
class ReadOverlap {
	private int indexSequence1;
	private int indexSequence2;
	private int first1;
	private int last1;
	private int first2;
	private int last2;
	private boolean negativeStrand;
	public ReadOverlap(int indexSequence1, int first1, int last1, int indexSequence2, int first2, int last2, boolean negativeStrand) {
		super();
		this.indexSequence1 = indexSequence1;
		this.indexSequence2 = indexSequence2;
		this.first1 = first1;
		this.last1 = last1;
		this.first2 = first2;
		this.last2 = last2;
		this.negativeStrand = negativeStrand;
	}
	public boolean addKmerAlignment(int idxSeq2, ReadAlignment aln) {
		if(aln==null) return false;
		if(idxSeq2!=indexSequence2) return false;
		if(aln.isNegativeStrand()!=negativeStrand) return false;
		int kmerFirst = aln.getMateFirst();
		int kmerLast = kmerFirst + OverlapGraph.SEARCH_KMER_LENGTH;
		int expectedAlnFirst =  last2 + (kmerFirst-last1);
		int expectedAlnLast =  last2 + (kmerLast-last1);
		if (negativeStrand) {
			expectedAlnFirst =  first2 - (kmerLast-last1);
			expectedAlnLast =  first2 - (kmerFirst-last1);
		}
		if(Math.abs(aln.getFirst()-expectedAlnFirst)>5) return false;
		if(Math.abs(aln.getLast()-expectedAlnLast)>5) return false;
		last1 = kmerLast;
		if(negativeStrand) {
			first2 = aln.getFirst();
		} else {
			last2 = aln.getLast();
		}
		return true;
	}
	public int getIndexSequence1() {
		return indexSequence1;
	}
	public int getIndexSequence2() {
		return indexSequence2;
	}
	public int getFirst1() {
		return first1;
	}
	public int getLast1() {
		return last1;
	}
	public int getFirst2() {
		return first2;
	}
	public int getLast2() {
		return last2;
	}
	public boolean isNegativeStrand() {
		return negativeStrand;
	}
	
}
