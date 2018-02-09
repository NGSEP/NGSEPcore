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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.DNASequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.KmersCounter;
import ngsep.sequences.LimitedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;

/**
 * @author Jorge Duitama
 *
 */
public class OverlapGraph {
	private List<LimitedSequence> sequences = new ArrayList<>();
	private List<List<ReadOverlap>> overlaps;
	
	public void addSequence (LimitedSequence sequence) {
		sequences.add(sequence);
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
		FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList seqsQl = handler.loadSequences(args[0]);
		for(QualifiedSequence seq:seqsQl) {
			totalLength+=seq.getLength();
			instance.addSequence((LimitedSequence) seq.getCharacters());
		}
		
		long time2 = System.currentTimeMillis();
		
		//System.out.println("Created "+instance.sequences.size()+" random DNA sequences. Total length: "+totalLength +" time: "+ (time2 - time1));
		System.out.println("Loaded "+instance.sequences.size()+" sequences from "+args[0]+". Total length: "+totalLength +" time: "+ (time2 - time1));
		//System.out.println("First sequence: "+instance.sequences.get(0));
		/*time1 = time2;
		Collections.sort(instance.sequences);
		time2 = System.currentTimeMillis();
		System.out.println("Sorting time: "+ (time2 - time1));*/
		time1=time2;
		FMIndex index = new FMIndex();
		index.loadUnnamedSequences("", instance.sequences);
		time2 = System.currentTimeMillis();
		System.out.println("Created FMIndex in time: "+ (time2 - time1));
		
		for(int s=0;s<100;s++) {
			time1=time2;
			CharSequence sequence = instance.sequences.get(s);
			System.out.println("Processing sequence: "+s+" length: "+sequence.length());
			CharSequence [] kmers = KmersCounter.extractKmers(sequence, 15, true);
			Map<Integer,Integer> seqHits = new HashMap<>();
			for(int i=0;i<kmers.length;i++) {
				String kmer = kmers[i].toString();
				List<ReadAlignment> regions = index.search(kmer);
				//System.out.println("Hits kmer "+kmer+": "+regions.size());
				for(ReadAlignment aln:regions) {
					int k = Integer.parseInt(aln.getSequenceName().substring(1));
					//String seqOut = instance.sequences.get(k).subSequence(region.getFirst(), region.getLast()+1).toString();
					//if(!seqOut.equals(kmer)) throw new RuntimeException ("Hit error in sequence: "+k+" pos: "+region.getFirst()+" search "+kmer+" found: "+seqOut);
					Integer count = seqHits.get(k);
					if(count==null) count=0;
					count++;
					seqHits.put(k, count);
				}
			}
			Set<Integer> selected = new HashSet<>();
			for(int k:seqHits.keySet()) {
				if(seqHits.get(k)>kmers.length/5) selected.add(k);
			}
			time2 = System.currentTimeMillis();
			System.out.println("Found "+selected.size()+" matching sequences for sequence "+s+" of length: "+sequence.length()+". search time: "+ (time2 - time1));
		}
	}
}
class ReadOverlap {
	private int indexSequence2;
	private int first1;
	private int last1;
	private int first2;
	private int last2;
	private boolean negativeStrand;
	public ReadOverlap(int indexSequence2, int first1, int last1, int first2, int last2, boolean negativeStrand) {
		super();
		this.indexSequence2 = indexSequence2;
		this.first1 = first1;
		this.last1 = last1;
		this.first2 = first2;
		this.last2 = last2;
		this.negativeStrand = negativeStrand;
	}
}
