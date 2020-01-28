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
package ngsep.alignments;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.alignments.io.ReadAlignmentFileWriter;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.ReferenceGenomeFMIndex;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndexUngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;

/**
 * Program to align reads to a reference genome
 * @author German Andrade
 * @author Jorge Duitama
 */
public class ReadsAligner {

	// Constants for default values
	
	
	public static final byte INPUT_FORMAT_FASTQ=KmersExtractor.INPUT_FORMAT_FASTQ;
	public static final byte INPUT_FORMAT_FASTA=KmersExtractor.INPUT_FORMAT_FASTA;
	public static final int DEF_KMER_LENGTH = KmersExtractor.DEF_KMER_LENGTH;
	public static final double DEF_MIN_PROPORTION_KMERS = 0.5;
	public static final int DEF_MIN_INSERT_LENGTH=0;
	public static final int DEF_MAX_INSERT_LENGTH=1000;
	
	public static final int DEF_MAX_ALIGNMENTS=3;
	public static final int MAX_SPACE_BETWEEN_KMERS = 50;
	
	// Logging and progress
	private Logger log = Logger.getLogger(ReadsAligner.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	// Parameters
	private String inputFile = null;
	private String inputFile2 = null;
	private String outputFile = null;
	private String fmIndexFile = null;
	private String knownSTRsFile = null;
	private byte inputFormat = INPUT_FORMAT_FASTQ;
	private int kmerLength = DEF_KMER_LENGTH;
	private double minProportionKmers = DEF_MIN_PROPORTION_KMERS;
	private int minInsertLength = DEF_MIN_INSERT_LENGTH;
	private int maxInsertLength = DEF_MAX_INSERT_LENGTH;
	
	// Model attributes
	private Map<String, List<GenomicRegion>> knownSTRs;

	private boolean onlyPositiveStrand = false;
	
	private boolean runFullAlignment = true;

	private ReferenceGenomeFMIndex fMIndex;
	
	private Set<String> repetitiveKmers = new HashSet<String>();
	
	// Get and set methods
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) { 
		this.progressNotifier = progressNotifier;
	}
	
	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}
	
	public String getInputFile2() {
		return inputFile2;
	}
	public void setInputFile2(String inputFile2) {
		this.inputFile2 = inputFile2;
	}
	
	public String getFmIndexFile() {
		return fmIndexFile;
	}
	public void setFmIndexFile(String fmIndexFile) {
		this.fmIndexFile = fmIndexFile;
	}
	
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	
	public String getKnownSTRsFile() {
		return knownSTRsFile;
	}
	public void setKnownSTRsFile(String knownSTRsFile) {
		this.knownSTRsFile = knownSTRsFile;
	}
	public byte getInputFormat() {
		return inputFormat;
	}
	public void setInputFormat(byte inputFormat) {
		this.inputFormat = inputFormat;
	}
	
	public int getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
	}
	public void setKmerLength(String value) {
		setKmerLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getMinInsertLength() {
		return minInsertLength;
	}
	public void setMinInsertLength(int minInsertLength) {
		this.minInsertLength = minInsertLength;
	}
	public void setMinInsertLength(String value) {
		setMinInsertLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getMaxInsertLength() {
		return maxInsertLength;
	}
	public void setMaxInsertLength(int maxInsertLength) {
		this.maxInsertLength = maxInsertLength;
	}
	public void setMaxInsertLength(String value) {
		setMaxInsertLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public double getMinProportionKmers() {
		return minProportionKmers;
	}
	public void setMinProportionKmers(double minProportionKmers) {
		this.minProportionKmers = minProportionKmers;
	}
	public void setMinProportionKmers(String value) {
		this.setMinProportionKmers((double)OptionValuesDecoder.decode(value, Double.class));
	}
	
	public Map<String, List<GenomicRegion>> getKnownSTRs() {
		return knownSTRs;
	}
	public void setKnownSTRs(Map<String, List<GenomicRegion>> knownSTRs) {
		this.knownSTRs = knownSTRs;
	}
	
	public ReadsAligner(String fMIndexFile) throws IOException {
		fMIndex = ReferenceGenomeFMIndex.loadFromBinaries(fMIndexFile);
	}

	public ReadsAligner() {	
		
	}

	public static void main(String[] args) throws Exception 
	{
		ReadsAligner instance = new ReadsAligner();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	
	public void run () throws IOException {
		
		if(fmIndexFile!=null) {
			log.info("Loading reference index from file: "+fmIndexFile);
			fMIndex = ReferenceGenomeFMIndex.loadFromBinaries(fmIndexFile);
		} else if (fMIndex!=null) {
			log.info("Aligning reads using built index with "+fMIndex.getSequencesMetadata().size()+" sequences");
		} else {
			throw new IOException("The genome index file is a required parameter");
		}
		
		QualifiedSequenceList sequences = fMIndex.getSequencesMetadata();
		PrintStream out = System.out; 
		if(outputFile!=null) out = new PrintStream(outputFile); 
		try (ReadAlignmentFileWriter writer = new ReadAlignmentFileWriter(sequences, out)){
			if(inputFile!=null && inputFile2!=null) {
				log.info("Aligning paired end reads from files: "+inputFile + " and "+inputFile2);
				alignReads(inputFile,inputFile2, writer);
			} else if (inputFile!=null) {
				log.info("Aligning single reads from file: "+inputFile);
				alignReads(inputFile, writer);
			} else if (inputFile2!=null ) {
				throw new IOException("The first input file is required for paired end alignment");
			} else {
				log.info("Aligning single reads from standard input");
				try (FastqFileReader reader = new FastqFileReader(System.in)) {
					reader.setSequenceType(DNAMaskedSequence.class);
					alignReads(reader, writer);
				}
				
			}
		}
	}
	
	public void loadSTRsFile(String strsFile) throws IOException {
		SimpleGenomicRegionFileHandler handler = new SimpleGenomicRegionFileHandler();
		knownSTRs=handler.loadRegionsAsMap(strsFile);
		flat();
	}

	private void flat() {
		Set<String>keys=knownSTRs.keySet();
		Iterator<String>it=keys.iterator();
		while(it.hasNext()) {
			String key = it.next();
			List<GenomicRegion> l =knownSTRs.get(key);
			List<GenomicRegion> newList = flat(l);
			while(isOverlappging(newList)) {
				newList=flat(newList);
			}
			knownSTRs.put(key, newList);
		}
	}

	private boolean isOverlappging(List<GenomicRegion>l) {
		for (int i = 0; i < l.size()-1; i++) {
			GenomicRegion current = l.get(i);
			GenomicRegion next = l.get(i+1);
			if(isOverlappgingSorted(current, next)) {
				return true;
			}
		}
		return false;
	}

	private boolean isOverlappgingSorted(GenomicRegion a,GenomicRegion b) {
		return b.getFirst()>=a.getFirst()&&b.getFirst()<=a.getLast();
	}

	private List<GenomicRegion> flat(List<GenomicRegion> l) {
		List<GenomicRegion> newList = new ArrayList<GenomicRegion>();
		for (int i = 0; i < l.size(); i++) {
			if(i+1==l.size())newList.add(l.get(i));
			else {
				GenomicRegion current =l.get(i);
				GenomicRegion next =l.get(i+1);
				if(isOverlappgingSorted(current, next)) {
					newList.add(new GenomicRegionImpl(current.getSequenceName(), current.getFirst(), Math.max(current.getLast(), next.getLast())));
					i++;
				}
				else {
					newList.add(current);
				}
			}	
		}

		return newList;
	}



	/**
	 * Aligns readsFile with the fMIndexFile
	 * @param fMIndexFile Binary file with the serialization of an FMIndex
	 * @param readsFile Fastq file with the reads to align
	 * @param out
	 * @throws IOException
	 */
	public void alignReads( String readsFile, ReadAlignmentFileWriter writer) throws IOException {
		
		if(inputFormat == INPUT_FORMAT_FASTQ) {
			try (FastqFileReader reader = new FastqFileReader(readsFile)) {
				//Load as DNAMaskedSequence to allow reverse complement
				reader.setSequenceType(DNAMaskedSequence.class);
				alignReads(reader,writer);
			}
		} else {
			//TODO: Implement fasta
		}
		
		
	}
	
	private void alignReads(FastqFileReader reader, ReadAlignmentFileWriter writer) throws IOException {
		if(knownSTRsFile!=null && !knownSTRsFile.isEmpty()) loadSTRsFile(knownSTRsFile);
		int totalReads = 0;
		int readsAligned = 0;
		int uniqueAlignments=0;
		long time = System.currentTimeMillis();
		Iterator<RawRead> it = reader.iterator();
		while(it.hasNext()) {
			RawRead read = it.next();
			List<ReadAlignment> alns = alignRead(read, true);
			//System.out.println("Alignments for: "+read.getName()+" "+alns.size());
			for(ReadAlignment aln:alns) writer.write(aln);
			if(alns.size()==0) {
				ReadAlignment alnNoMap = createUnmappedAlignment(read, false, false);
				writer.write(alnNoMap);
			}
			int numAlns = alns.size();
			totalReads++;
			if(numAlns>0) readsAligned++;
			if(numAlns==1) uniqueAlignments++;
			if(totalReads%100000==0) log.info("Processed "+totalReads+" reads. Aligned: "+readsAligned);
		}
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		log.info("Time: "+seconds+" seconds");
		
		log.info("Total reads: "+totalReads);
		log.info("Reads aligned: "+readsAligned);
		log.info("Unique alignments: "+uniqueAlignments);
		log.info("Overall alignment rate: "+(100.0*readsAligned/(double)totalReads)+"%");
	}

	/**
	 * Aligns readsFile with the fMIndexFile for pairend
	 * @param fMIndexFile Binary file with the serialization of an FMIndex
	 * @param readsFile1 Fastq file with the reads to align
	 * @param readsFile2 Fastq file with the reads to align
	 * @param out
	 * @throws IOException
	 */
	public void alignReads( String readsFile1, String readsFile2, ReadAlignmentFileWriter writer) throws IOException {
		if(knownSTRsFile!=null && !knownSTRsFile.isEmpty())loadSTRsFile(knownSTRsFile);
		int totalReads = 0;
		int readsAligned = 0;
		int proper = 0;
		int notProper = 0;
		int single = 0;
		int uniqueAlignments=0;
		long time = System.currentTimeMillis();
		try (FastqFileReader reader1 = new FastqFileReader(readsFile1); FastqFileReader reader2 = new FastqFileReader(readsFile2)) {
			reader1.setSequenceType(DNAMaskedSequence.class);
			reader2.setSequenceType(DNAMaskedSequence.class);
			//Load as DNAMaskedSequence to allow reverse complement
			Iterator<RawRead> it1 = reader1.iterator();
			Iterator<RawRead> it2 = reader2.iterator();
			while(it1.hasNext() && it2.hasNext()) {
				RawRead read1 = it1.next();
				RawRead read2 = it2.next();
				List<ReadAlignment> alns1 = alignRead(read1,false);
				for(ReadAlignment aln:alns1) aln.setFirstOfPair(true);
				List<ReadAlignment> alns2 = alignRead(read2,false);
				for(ReadAlignment aln:alns2) aln.setSecondOfPair(true);
				totalReads++;
				if(alns1.size()==0 || alns2.size()==0) {
					ArrayList<ReadAlignment> unMapped = processUnMapped(read1, alns1,read2,alns2);
					boolean mappedFound = false;
					for (int i = 0; i < Math.min(unMapped.size(),DEF_MAX_ALIGNMENTS+1); i++) {
						ReadAlignment aln = unMapped.get(i);
						if(!aln.isReadUnmapped()) mappedFound=true;
						writer.write(aln);
						
					}
					if(mappedFound) {
						readsAligned++;
						single++;
					}
					
				} else {
					boolean onlyProper=true;
					List<ReadAlignment> alns = new ArrayList<ReadAlignment>();
					List<ReadAlignmentPair> pairAlns = findPairs(alns1, alns2,onlyProper);
					if(pairAlns.isEmpty()) {
						pairAlns = findPairs(alns1, alns2,false);
						if(pairAlns.isEmpty()) {
							single++;
							for(int i=0;i<alns1.size();i++) {
								ReadAlignment aln1 = alns1.get(i);
								aln1.setPaired(true);
								aln1.setMateDifferentSequence(true);
								setMateInfo(aln1, alns2.get(0));
								if(i>0) aln1.setSecondary(true);
								alns.add(aln1);
							}
							for(int i=0;i<alns2.size();i++) {
								ReadAlignment aln2 = alns2.get(i);
								aln2.setPaired(true);
								aln2.setMateDifferentSequence(true);
								setMateInfo(aln2, alns1.get(0));
								if(i>0) aln2.setSecondary(true);
								alns.add(aln2);
							}
							
						}
						else {
							notProper++;
							addPairAlignments(alns, pairAlns);
							if(pairAlns.size()==1) uniqueAlignments++;
						}

					} else {
						proper++;
						addPairAlignments(alns, pairAlns);
						if(pairAlns.size()==1) uniqueAlignments++;
					}
					for(ReadAlignment aln:alns) writer.write(aln);
					readsAligned++;
					if(totalReads%100000==0) {
						log.info("Processed "+totalReads+" reads. Aligned: "+readsAligned);
						log.info("Reads aligned proper: "+proper);
						log.info("Reads aligned notProper: "+notProper);
						log.info("Reads aligned single: "+single);	
					} 
				}
			}

		}
		log.info("Total reads: "+totalReads);
		log.info("Reads aligned proper: "+proper);
		log.info("Reads aligned notProper: "+notProper);
		log.info("Reads aligned single: "+single);
		log.info("Reads aligned: "+readsAligned);
		log.info("Unique alignments: "+uniqueAlignments);
		log.info("Overall pairend alignment rate: "+(100.0*(proper+notProper)/(double)totalReads)+"%");
		log.info("Overall alignment rate: "+(100.0*readsAligned/(double)totalReads)+"%");
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		log.info("Time: "+seconds+" seconds");

	}

	private void addPairAlignments(List<ReadAlignment> alns, List<ReadAlignmentPair> pairAlns) {
		for (int i = 0; i < Math.min(pairAlns.size(),DEF_MAX_ALIGNMENTS); i++) {
			ReadAlignmentPair current = pairAlns.get(i);
			if(i>0) {
				current.getAln1().setSecondary(true);
				current.getAln2().setSecondary(true);
			}
			alns.add(current.getAln1());
			alns.add(current.getAln2());
		}
	}

	public List<ReadAlignmentPair> findPairs(List<ReadAlignment> alns1, List<ReadAlignment> alns2,boolean onlyProper){
		List<ReadAlignmentPair> pairEndAlns = new ArrayList<ReadAlignmentPair>();
		for (int i = 0; i < Math.min(alns1.size(),DEF_MAX_ALIGNMENTS); i++) {
			ReadAlignment aln1 = alns1.get(i);
			if(aln1.isPaired()) continue;
			ReadAlignmentPair alnPair = findPairForAlignment(aln1,alns2,onlyProper);
			if(alnPair!=null) {
				pairEndAlns.add(alnPair);
			}
		}
		return pairEndAlns;
	}

	public ReadAlignmentPair findPairForAlignment(ReadAlignment aln1, List<ReadAlignment> alns2,boolean onlyProper) {
		List<ReadAlignment> candidates = new ArrayList<ReadAlignment>();
		for (int i = 0; i < Math.min(alns2.size(), DEF_MAX_ALIGNMENTS); i++) {
			ReadAlignment current =alns2.get(i);
			if(!current.isPaired()) {
				if(isValidPair(aln1,current,onlyProper)) {
					candidates.add(current);
				}
			}
		}
		if(candidates.size()==1) {
			return buildPair(aln1, candidates.get(0), onlyProper);
		}
		else if (candidates.size()>1){
			return buildPair(aln1, getRandomReadAlignment(candidates), onlyProper);
		}
		return null;
	}

	public boolean isValidPair(ReadAlignment aln1, ReadAlignment aln2,boolean onlyProper) {
		if(!aln1.getSequenceName().equals(aln2.getSequenceName())) return false;
		int start1 = aln1.getFirst();
		int start2 = aln2.getFirst();
		int end1 = aln1.getLast();
		int end2 = aln2.getLast();

		int endMax = Math.max(end1, end2);
		int startMinimum = Math.min(start1, start2);
		if(endMax==end1 && aln1.isNegativeStrand() && startMinimum==start2&& aln2.isPositiveStrand()
				|| endMax==end2 && aln2.isNegativeStrand() && startMinimum==start1&& aln1.isPositiveStrand())
		{
			if(onlyProper) {
				int insertLength = endMax-startMinimum;
				return insertLength>=minInsertLength && insertLength<=maxInsertLength;
			}
			else{
				return true;
			}
		}
		return false;
	}

	private ReadAlignmentPair buildPair(ReadAlignment aln1, ReadAlignment aln2, boolean proper) {
		aln1.setPaired(true);
		
		aln2.setPaired(true);
		
		if(proper) {
			aln1.setProperPair(true);
			aln2.setProperPair(true);
		}
		int insertLength1 = aln1.getLast()-aln2.getFirst();
		int insertLength2 = aln2.getLast()-aln1.getFirst();
		if(insertLength1>insertLength2) {
			aln1.setInferredInsertSize(-insertLength1);
			aln2.setInferredInsertSize(insertLength1);
		} else {
			aln1.setInferredInsertSize(insertLength2);
			aln2.setInferredInsertSize(-insertLength2);
		}
		setMateInfo(aln1,aln2);
		setMateInfo(aln2,aln1);
		return new ReadAlignmentPair(aln1,aln2);
	}


	public GenomicRegion findTandemRepeat(ReadAlignment aln1) {
		if(knownSTRs==null)return null;
		List<GenomicRegion> l =knownSTRs.get(aln1.getSequenceName());
		return binaryContains(l, 0, l.size()-1, aln1);
	}

	private GenomicRegion binaryContains(List<GenomicRegion> l, int left, int rigth, GenomicRegion element) {
		if(rigth>=left) {
			int middle = left +(rigth-left)/2;
			GenomicRegion actual = l.get(middle);
			if(isOverlappgingSorted(actual, element)||isOverlappgingSorted(element, actual	)) {
				return actual;
			}
			if(l.get(middle).getFirst()>element.getFirst()) {
				return binaryContains(l, left, middle-1, element);
			}
			return binaryContains(l, middle+1,rigth, element);
		}
		return null;
	}

	private ArrayList<ReadAlignment> processUnMapped(RawRead read1, List<ReadAlignment> alns1, RawRead read2, List<ReadAlignment> alns2) {
		ArrayList<ReadAlignment> unMappedAlignments= new ArrayList<ReadAlignment>();
		ReadAlignment unmappedAln1 = null;
		ReadAlignment unmappedAln2 = null;
		if(alns1.size()==0) unmappedAln1 = createUnmappedAlignment(read1, true, true);
		if(alns2.size()==0) unmappedAln2 = createUnmappedAlignment(read2, true, false);
		if(unmappedAln1!= null && unmappedAln2!=null) {
			unmappedAln1.setMateUnmapped(true);
			unmappedAln2.setMateUnmapped(true);
			unMappedAlignments.add(unmappedAln1);
			unMappedAlignments.add(unmappedAln2);
		}
		else if(unmappedAln2!=null) readMappedmateUnMapped(unmappedAln2, alns1, unMappedAlignments);
		else if(unmappedAln1!=null) readMappedmateUnMapped(unmappedAln1, alns2, unMappedAlignments);
		return unMappedAlignments;
	}

	private void readMappedmateUnMapped(ReadAlignment unmappedAln, List<ReadAlignment> alns, List<ReadAlignment> unMappedAlignments) {
		unMappedAlignments.add(unmappedAln);
		
		for(int i=0;i<alns.size();i++) {
			ReadAlignment aln =  alns.get(i);
			aln.setPaired(true);
			aln.setMateUnmapped(true);
			if(i==0) {
				setMateInfo(unmappedAln, aln);
			} else {
				aln.setSecondary(true);
			}
			unMappedAlignments.add(aln);
		}
		
	}

	private ReadAlignment createUnmappedAlignment(RawRead read, boolean paired, boolean first) {
		ReadAlignment alnNoMap = new ReadAlignment(null, 0, 0, read.getLength(), ReadAlignment.FLAG_READ_UNMAPPED);
		alnNoMap.setReadName(read.getName());
		alnNoMap.setReadCharacters(read.getCharacters());
		alnNoMap.setQualityScores(read.getQualityScores());
		alnNoMap.setPaired(paired);
		if(paired) {
			alnNoMap.setFirstOfPair(first);
			alnNoMap.setSecondOfPair(!first);
		}
		return alnNoMap;
	}

	private ReadAlignment setMateInfo(ReadAlignment alignment,ReadAlignment mate) {
		alignment.setMateFirst(mate.getFirst());
		alignment.setMateSequenceName(mate.getSequenceName());
		alignment.setMateNegativeStrand(mate.isNegativeStrand());
		return alignment;
	}

	private static ReadAlignment getRandomReadAlignment(List<ReadAlignment> alns) {
		Random r = new Random();
		return alns.get(r.nextInt(alns.size())) ;
	}

	private ReadAlignment buildAln(CharSequence query, String qualityScores, String sequenceName, int first, int last,String cigar, double alnQual) {
		ReadAlignment aln = new ReadAlignment(sequenceName, first, last, query.length(), 0);
		aln.setReadCharacters(query);
		aln.setQualityScores(qualityScores);
		if(cigar!=null)aln.setCigarString(cigar);
		aln.setAlignmentQuality((short) Math.round(alnQual));
		//verify last exists
		if(!fMIndex.isValidPosition(sequenceName,aln.getLast())) return null;
		return aln;
	}

	public ReadAlignment verifyShortTandemRepeats(ReadAlignment aln, CharSequence read, String qualityScores, GenomicRegion region) {
		//Case 1
		if(aln.getFirst()<region.getFirst() && aln.getLast()<region.getLast()) {
			return alnInLeft(aln, read, qualityScores, region);
		}
		//Case 2
		if(aln.getFirst()>region.getFirst()&&aln.getLast()>region.getLast()) {
			return alnInRigth(aln, read, qualityScores, region);
		}
		//Case 3
		if(aln.getFirst()<region.getFirst()&&aln.getLast()>region.getLast()) {
			int overlapLength = region.length();

			int aln1First =Math.max(1,aln.getFirst());
			int aln1Last=region.getFirst()-1;
			String aln1SequenceName=aln.getSequenceName();
			if(aln1First>=aln1Last)return null;
			ReadAlignment aln1=alnInMid(read,qualityScores, aln1First, aln1Last,aln1SequenceName,overlapLength,false);

			int aln2First=Math.max(1,region.getLast())+1;
			int aln2Last=aln.getLast();
			if(aln2First>=aln2Last)return null;

			String aln2SequenceName=aln.getSequenceName();
			ReadAlignment aln2=alnInMid( read, qualityScores,  aln2First,  aln2Last, aln2SequenceName,overlapLength+aln1Last-aln1First+1,true);

			if(aln1==null||aln2==null)return null;
			int a = aln1.length()+aln2.length()+region.length();
			int b = aln.length();
			int c=a-b;
			//case 3a
			int first = aln1.getFirst();
			int last = aln2.getLast();
			String cigar = aln1.getCigarString()+""+region.length()+"M"+aln2.getCigarString();
			double alnQual = (aln1.getAlignmentQuality()+aln2.getAlignmentQuality())/2;

			//case 3b
			if(c>0) {
				System.out.println("c: "+c);
			}
			//case 3c
			else if(c<0) {
				int del1=aln1Last-aln1.getLast();
				int del2 = aln2.getFirst()-aln2First;
				cigar = aln1.getCigarString();
				cigar+=addIfnotEmpty(del1,"D");
				cigar+=addIfnotEmpty(region.length(),"M");
				cigar+=addIfnotEmpty(del2,"D");
				cigar+=aln2.getCigarString();
			}
			return buildAln(read, qualityScores, aln1.getSequenceName(), first, last, cigar, alnQual);
		}
		return null;
	}

	private String addIfnotEmpty(int del1, String string) {
		if(del1>0)return del1+string;
		return "";
	}

	private ReadAlignment alnInLeft(ReadAlignment aln, CharSequence read, String qualityScores, GenomicRegion region) {
		int overlapStart = region.getFirst();
		int overlapEnd = aln.getLast();
		int overlapLength = overlapEnd-overlapStart+1;
		String sequenceName=aln.getSequenceName();
		int first = Math.max(1, aln.getFirst());
		int last = overlapStart;
		CharSequence refSeq = fMIndex.getSequence(sequenceName, first, last);
		if(overlapLength>1 && overlapLength<90 && refSeq != null)
		{
			String readNoSTR=read.subSequence(0, read.length()-overlapLength).toString();
			AlignmentResult result = smithWatermanLocalAlignment(readNoSTR,refSeq);
			last = first+result.getSubjectLastIdx();
			first = first + result.getSubjectStartIdx();
			String cigar = result.getCigarString()+overlapLength+"S";
			double alnQual = 100.0* (read.length() - result.getDistance())/read.length();
			ReadAlignment newaln = buildAln(read, qualityScores, sequenceName, first, last, cigar, alnQual);
			return newaln;
		}
		return null;
	}


	private ReadAlignment alnInMid(CharSequence read,String qualityScores, int pFirst, int pLast,String sequenceName, int overlapLength, boolean rigth) {
		try {
			String readNoSTR=read.subSequence(0,pLast-pFirst+1).toString();
			if(rigth)readNoSTR=read.subSequence(overlapLength,read.length()).toString();
			CharSequence refSeq = fMIndex.getSequence(sequenceName, pFirst, pLast);
			AlignmentResult result = smithWatermanLocalAlignment(readNoSTR,refSeq);
			int last = pFirst+result.getSubjectLastIdx();
			int first = pFirst + result.getSubjectStartIdx();
			String cigar = result.getCigarString();
			double alnQual = 100.0* (read.length() - result.getDistance())/read.length();
			ReadAlignment newaln = buildAln(readNoSTR, qualityScores, sequenceName, first, last, cigar, alnQual);
			return newaln;
		}
		catch(Exception e) {
			return null;
		}
	}



	private ReadAlignment alnInRigth(ReadAlignment aln, CharSequence read, String qualityScores, GenomicRegion region) {
		int overlapStart=aln.getFirst();
		int overlapEnd=region.getLast();
		int overlapLength = overlapEnd-overlapStart+1;
		String sequenceName=aln.getSequenceName();
		int first =Math.max(1, region.getLast()+1);
		int last = aln.getLast();
		CharSequence refSeq = fMIndex.getSequence(sequenceName, first, last);
		if(overlapLength>0 && overlapLength<90 && refSeq != null)
		{
			String readNoSTR=read.subSequence(overlapLength, read.length()).toString();
			AlignmentResult result = smithWatermanLocalAlignment(readNoSTR,refSeq);
			last = first+result.getSubjectLastIdx();
			first = first + result.getSubjectStartIdx();
			String cigar = overlapLength+"S"+result.getCigarString();
			double alnQual = 100.0* (read.length() - result.getDistance())/read.length();
			ReadAlignment newaln = buildAln(read, qualityScores, sequenceName, first, last, cigar, alnQual);
			return newaln;
		}
		return null;
	}
	/**
	 * Old and stable
	 * @param read
	 * @return
	 */

	public List<ReadAlignment> alignRead(RawRead read, boolean asignSecondaryStatus) {
		List<ReadAlignment> alignments = kmerBasedInexactSearchAlgorithm(read);
		return filterAlignments(alignments, asignSecondaryStatus);
	}

	/**
	 * First approach to allow inexact search
	 * It iterates the Sequences of the genome and if there is at least MIN_ACCURACY percentage of the kmers
	 * it allow the alignment with the first an the last position of the kmers ocurrence
	 * @return 
	 */
	private List<ReadAlignment> kmerBasedInexactSearchAlgorithm (RawRead read) {
		List<ReadAlignment> alns = new ArrayList<>();
		DNAMaskedSequence readSeq = (DNAMaskedSequence)read.getCharacters();
		String qual = read.getQualityScores();
		alns.addAll(kmerBasedInexactSearchAlgorithm(readSeq, qual, read.getName()));

		if(!onlyPositiveStrand) {
			readSeq = readSeq.getReverseComplement();
			qual = new StringBuilder(qual).reverse().toString();
			List<ReadAlignment> alnsR = kmerBasedInexactSearchAlgorithm(readSeq, qual, read.getName());
			for (ReadAlignment aln:alnsR) {
				aln.setNegativeStrand(true);
			}
			alns.addAll(alnsR);
		}
		return alns;
	}
	/**
	 * First approach to allow inexact search
	 * It iterates the Sequences of the genome and if there is at least MIN_ACCURACY percentage of the kmers
	 * it allow the alignment with the first an the last position of the kmers ocurrence
	 * Only tries to align the given quey in the positive strand
	 * @return List<ReadAlignment>
	 */
	private List<ReadAlignment> kmerBasedInexactSearchAlgorithm (CharSequence query, String qualityScores, String readName) 
	{
		Map<Integer,CharSequence> kmersMap = KmersExtractor.extractKmersAsMap(query, kmerLength, kmerLength, true, true, true);
		List<ReadAlignment> finalAlignments =  new ArrayList<>();
		//System.out.println("Query: "+query.toString()+" kmers: "+kmers.size());
		int kmersCount=kmersMap.size();
		if(kmersCount==0) return finalAlignments;
		List<FMIndexUngappedSearchHit> initialKmerHits = searchKmers (kmersMap);
		List<KmerHitsCluster> clusteredKmerHits = clusterKmerHits(query, initialKmerHits); 
		//System.out.println("Initial kmer hits: "+initialKmerHits.size()+" Clusters: "+clusteredKmerHits.size());
		Collections.sort(clusteredKmerHits, (o1, o2) -> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());

		int kmersMaxCluster = 0;
		for (int i=0;i<clusteredKmerHits.size() && i<DEF_MAX_ALIGNMENTS;i++) {
			KmerHitsCluster cluster = clusteredKmerHits.get(i);
			//System.out.println("Processing cluster "+i+" spanning "+cluster.getSequenceName()+":"+cluster.getFirst()+"-"+cluster.getLast()+" Num kmers: "+cluster.getNumDifferentKmers()+" consistent: "+cluster.isAllConsistent());
			if(i==0) kmersMaxCluster = cluster.getNumDifferentKmers();
			else if (2*cluster.getNumDifferentKmers()<kmersMaxCluster) break;
			ReadAlignment readAln = createNewAlignmentFromConsistentKmers(cluster, kmersCount, query, qualityScores);
			if(readAln!=null) {
				readAln.setReadName(readName);
				finalAlignments.add(readAln);
			}
		}
		//System.out.println("Found "+finalAlignments.size()+" alignments for query: "+query.toString());
		return finalAlignments;
	}


	/**
	 * Searches the given kmers in the fmIndex 
	 * @param kmers to search
	 * @return List of alignments of each kmer. The read number of each alignment contains the kmer number.
	 */
	private List<FMIndexUngappedSearchHit> searchKmers(Map<Integer,CharSequence> kmersMap) {
		List<FMIndexUngappedSearchHit> answer = new ArrayList<>();
		for (int start:kmersMap.keySet()) {
			String kmer = kmersMap.get(start).toString();
			if(repetitiveKmers.contains(kmer)) continue;
			List<FMIndexUngappedSearchHit> kmerHits=fMIndex.exactSearch(kmer);
			if(kmerHits.size()>10) {
				repetitiveKmers.add(kmer);
				continue;
			}
			for(FMIndexUngappedSearchHit hit:kmerHits) {
				hit.setQueryIdx(start);
				answer.add(hit);
			}
		}
		return answer;
	}

	private List<KmerHitsCluster> clusterKmerHits(CharSequence query, List<FMIndexUngappedSearchHit> initialKmerHits) {
		List<KmerHitsCluster> clusters = new ArrayList<>();
		Map<String,List<FMIndexUngappedSearchHit>> hitsBySubjectName = new LinkedHashMap<String, List<FMIndexUngappedSearchHit>>();
		for(FMIndexUngappedSearchHit hit:initialKmerHits) {
			List<FMIndexUngappedSearchHit> hitsSeq = hitsBySubjectName.computeIfAbsent(hit.getSequenceName(), k -> new ArrayList<>());
			hitsSeq.add(hit);
		}
		for(List<FMIndexUngappedSearchHit> hitsSeq:hitsBySubjectName.values()) {
			Collections.sort(hitsSeq, new Comparator<FMIndexUngappedSearchHit>() {

				@Override
				public int compare(FMIndexUngappedSearchHit hit0, FMIndexUngappedSearchHit hit1) {
					return hit0.getStart()-hit1.getStart();
				}
			});
			clusters.addAll(clusterSequenceKmerAlns(query, hitsSeq));
		}
		return clusters;
	}
	private List<KmerHitsCluster> clusterSequenceKmerAlns(CharSequence query, List<FMIndexUngappedSearchHit> sequenceHits) {
		List<KmerHitsCluster> answer = new ArrayList<>();
		//System.out.println("Alns to cluster: "+sequenceAlns.size());
		KmerHitsCluster cluster=null;
		for(FMIndexUngappedSearchHit kmerHit:sequenceHits) {
			if(cluster==null || !cluster.addKmerHit(kmerHit, 0)) {
				cluster = new KmerHitsCluster(query, kmerHit);
				answer.add(cluster);
			}
		}
		return answer;
	}


	private ReadAlignment createNewAlignmentFromConsistentKmers(KmerHitsCluster cluster, int totalKmers, CharSequence query, String qualityScores) {
		int numDiffKmers = cluster.getNumDifferentKmers();
		double prop = (double) numDiffKmers/totalKmers;
		if(prop<minProportionKmers) return null;
		String sequenceName = cluster.getSequenceName();
		int first = cluster.getFirst();
		int last = cluster.getLast();
		int length = last-first+1;
		//System.out.println("Reference region from k-mers: "+first+"-"+last+" all consistent: "+cluster.isAllConsistent()+" firstaln: "+cluster.isFirstKmerPresent()+" last aln: "+cluster.isLastKmerPresent());
		String cigar = query.length()+"M";
		if(length>query.length()) {
			cigar+=(length-query.length())+"D";
		}
		double alnQual = 100.0;
		ReadAlignment aln = buildAln(query, qualityScores, sequenceName, first, last, cigar, alnQual);
		//System.out.println("Built alignment at "+sequenceName+":"+first+"-"+last+" CIGAR: "+cigar+" Aln: "+aln);
		if(aln!=null) {
			GenomicRegion region =findTandemRepeat(aln);
			if(region!=null) {
				ReadAlignment newaln=verifyShortTandemRepeats(aln,query,qualityScores,region);
				System.out.println("Found overlapping tandem repeat at "+region.getSequenceName()+":"+region.getFirst()+"-"+region.getLast()+" new aln: "+newaln);
				if(newaln!=null) return newaln;
			}
		}
		if (!cluster.isAllConsistent()) {
			if(!runFullAlignment) return null;
			//Perform smith waterman
			if(!cluster.isFirstKmerPresent()) first -=10;
			first = Math.max(1, first);
			if(!cluster.isLastKmerPresent()) last+=10;
			last = Math.min(fMIndex.getReferenceLength(sequenceName), cluster.getLast());
			CharSequence refSeq = fMIndex.getSequence(sequenceName, first, last);
			if(refSeq == null) return null;
			AlignmentResult result = smithWatermanLocalAlignment(query.toString(),refSeq);
			//TODO: Make better score
			if(result.getDistance()>0.5*query.length()) return null;
			//Last must be updated before first
			last = first+result.getSubjectLastIdx();
			first = first + result.getSubjectStartIdx();
			cigar = result.getCigarString();
			alnQual = 100.0* (query.length() - result.getDistance())/query.length();
			aln = buildAln(query, qualityScores, sequenceName, first, last, cigar, alnQual);
		}
		return aln;
	}

	private AlignmentResult smithWatermanLocalAlignment(CharSequence query, CharSequence subject) {
		//Matrix for dynamic programming saves the lowest weight from 0,0 to i,j
		int[][] scores = new int[query.length()+1][subject.length()+1]; 
		int minLastRow = scores.length*scores[0].length;
		int minJLastRow = -1;
		for (int i = 0; i < scores.length; i++) {	
			scores[i][0] = i==0?0:scores[i-1][0]+1;
			for (int j = 1; j < scores[0].length; j++) 
			{
				if(i==0)
				{
					scores[i][j]=0;
				}
				//This case can be reached from upper (i-1,j)
				//From left (i,j-1)
				//or from upper left diagonal (i-1,j-1)
				else
				{
					int d = query.charAt(i-1)!=subject.charAt(j-1)?1:0;
					int min = scores[i-1][j-1]+d;
					int a = 1 + scores[i-1][j];
					if(min>a) min= a;
					int b = 1 + scores[i][j-1];
					if(min>b) min=b;
					scores[i][j]=min;
				}
				if(i==scores.length-1 && (j==1 || minLastRow>=scores[i][j])) {
					minLastRow = scores[i][j];
					minJLastRow = j; 
				}
			}
		}
		//At this point we have the minimum cost to reach the corner (A.length-1,A[0].length-1)
		//Now we have to go back and remember the decisions from the minimum of the last row
		int i = scores.length-1;
		int j= minJLastRow;
		//System.out.println("Subject length: "+subject.length()+" minJ last row: "+j+" minscore: "+minLastRow);
		AlignmentResult result = new AlignmentResult();
		result.setSubjectLastIdx(minJLastRow-1);
		result.setDistance(minLastRow);
		//The algorithm keeps going back until reach origin
		while(i>0 && j>0) {
			//We want the path with lowest cost
			int d = query.charAt(i-1)!=subject.charAt(j-1)?1:0;
			int sD = scores[i-1][j-1]+d;
			int sA = 1 + scores[i-1][j];
			if( sD == scores[i][j]) {
				result.addBacktrack(ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH));
				i--;
				j--;
			} else if( sA == scores[i][j]) {
				result.addBacktrack(ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_INSERTION));
				i--;
			} else {
				result.addBacktrack(ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_DELETION));
				j--;
			}
		}
		result.setSubjectStartIdx(j);
		//System.out.println("Subject start: "+result.getSubjectStartIdx());
		while(i>0) {
			result.addBacktrack(ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_INSERTION));
			i--;
		}
		return result;
	}

	private List<ReadAlignment> filterAlignments(List<ReadAlignment> alignments, boolean assignSecondary) {
		if (alignments.size()==0) return alignments;
		Collections.sort(alignments, new Comparator<ReadAlignment>() {
			@Override
			public int compare(ReadAlignment aln1, ReadAlignment aln2) {
				return aln2.getAlignmentQuality() - aln1.getAlignmentQuality();
			}
		});
		short bestQual = alignments.get(0).getAlignmentQuality();
		//TODO. Investigate alignment score
		double threshold = 0.8*bestQual;
		List<ReadAlignment> filteredAlignments = new ArrayList<>();
		int n = alignments.size();
		for (int i=0;i<n;i++) {
			ReadAlignment aln = alignments.get(i);
			if(aln.getAlignmentQuality()<threshold) break;
			if(assignSecondary && i>0) aln.setSecondary(true);
			filteredAlignments.add(aln);
		}
		if(filteredAlignments.size()>1) {
			for(ReadAlignment aln:filteredAlignments) aln.setAlignmentQuality((short) Math.round(0.1*aln.getAlignmentQuality()));
		}
		return filteredAlignments;
	}
}





