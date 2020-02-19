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
import ngsep.sequences.LimitedSequence;
import ngsep.sequences.PairwiseAlignmentAffineGap;
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
	public static final int DEF_KMER_LENGTH = 25;
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
				int insertLength = endMax-startMinimum+1;
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
		int insertLength1 = aln1.getLast()-aln2.getFirst()+1;
		int insertLength2 = aln2.getLast()-aln1.getFirst()+1;
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
		if(first <=0) return null;
		ReadAlignment aln = new ReadAlignment(sequenceName, first, last, query.length(), 0);
		aln.setReadCharacters(query);
		aln.setQualityScores(qualityScores);
		if(cigar!=null)aln.setCigarString(cigar);
		aln.setAlignmentQuality((short) Math.round(alnQual));
		//verify last exists
		if(!fMIndex.isValidPosition(sequenceName,aln.getLast())) return null;
		return aln;
	}

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
		alns.addAll(kmerBasedInexactSearchAlgorithm(readSeq.toString(), qual, read.getName()));

		if(!onlyPositiveStrand) {
			readSeq = readSeq.getReverseComplement();
			qual = new StringBuilder(qual).reverse().toString();
			List<ReadAlignment> alnsR = kmerBasedInexactSearchAlgorithm(readSeq.toString(), qual, read.getName());
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
	private List<ReadAlignment> kmerBasedInexactSearchAlgorithm (String query, String qualityScores, String readName) 
	{
		Map<Integer,CharSequence> kmersMap = KmersExtractor.extractKmersAsMap(query, kmerLength, kmerLength, true, true, true);
		List<ReadAlignment> finalAlignments =  new ArrayList<>();
		//System.out.println("Query: "+query.toString()+" kmers: "+kmersMap.size());
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
		//System.out.println("Found "+finalAlignments.size()+" alignments for query: "+query);
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

	private List<KmerHitsCluster> clusterKmerHits(String query, List<FMIndexUngappedSearchHit> initialKmerHits) {
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
	
	private List<KmerHitsCluster> clusterSequenceKmerAlns(String query, List<FMIndexUngappedSearchHit> sequenceHits) {
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

	private ReadAlignment createNewAlignmentFromConsistentKmers(KmerHitsCluster cluster, int totalKmers, String query, String qualityScores) {
		int numDiffKmers = cluster.getNumDifferentKmers();
		double prop = (double) numDiffKmers/totalKmers;
		if(prop<minProportionKmers) return null;
		String sequenceName = cluster.getSequenceName();
		int first = cluster.getFirst();
		int last = cluster.getLast();
		int lastPerfect = first+query.length()-1;
		//System.out.println("Reference region from k-mers: "+first+"-"+last+" all consistent: "+cluster.isAllConsistent()+" lastPerfect: "+lastPerfect+" firstaln: "+cluster.isFirstKmerPresent()+" last aln: "+cluster.isLastKmerPresent());
		String cigar = query.length()+"M";
		double alnQual = 100.0;
		ReadAlignment aln = buildAln(query, qualityScores, sequenceName, first, lastPerfect, cigar, alnQual);
		if(aln!=null) {
			//System.out.println("Built alignment at "+sequenceName+":"+aln.getFirst()+"-"+aln.getLast()+" CIGAR: "+cigar);
			GenomicRegion region =findTandemRepeat(sequenceName,first,last);
			if(region!=null) {
				ReadAlignment newaln=verifyShortTandemRepeats(aln.getSequenceName(),aln.getFirst(), aln.getLast(),query,qualityScores,region);
				//System.out.println("Found overlapping tandem repeat at "+region.getSequenceName()+":"+region.getFirst()+"-"+region.getLast()+" new aln: "+newaln);
				if(newaln!=null) return newaln;
			}
			if(cluster.isAllConsistent()) {
				int mismatches = countMismatches (query, aln);
				//System.out.println("Mismatches alignment at "+aln.getSequenceName()+":"+aln.getFirst()+"-"+aln.getLast()+": "+mismatches);
				if(mismatches<3) return aln;
			}
		}
		if(!runFullAlignment) return null;
		//Perform smith waterman
		if(!cluster.isFirstKmerPresent()) first -=10;
		first = Math.max(1, first);
		if(!cluster.isLastKmerPresent()) last+=10;
		last = Math.min(fMIndex.getReferenceLength(sequenceName), cluster.getLast());
		CharSequence refSeq = fMIndex.getSequence(sequenceName, first, last);
		if(refSeq == null) return null;
		//System.out.println("Aligning reference from "+first+" to "+last+ " to query. length: "+refSeq.length());
		PairwiseAlignmentWithCigar pAln = new PairwiseAlignmentWithCigar(query, refSeq.toString(), false);
		//TODO: Make better score
		//System.out.println("Pairwise alignment found from relative : "+pAln.getSubjectStartIdx()+" to "+pAln.getSubjectLastIdx()+" CIGAR:" +pAln.getCigar()+" ");
		if(pAln.getMismatches()>0.05*query.length()) return null;
		//Last must be updated before first
		last = first + pAln.getSubjectLastIdx();
		first = first + pAln.getSubjectStartIdx();
		cigar = pAln.getCigar();
		//System.out.println("New genomic coordinates : "+first+"-"+last+" CIGAR:" +cigar);
		alnQual = 100.0* (query.length() - pAln.getMismatches())/query.length();
		return buildAln(query, qualityScores, sequenceName, first, last, cigar, alnQual);
	}

	private int countMismatches(CharSequence query, ReadAlignment aln) {
		int mismatches = 0;
		String refSeq = fMIndex.getSequence(aln.getSequenceName(), aln.getFirst(), aln.getLast()).toString();
		for (int i=0;i<query.length() && i<refSeq.length();i++ ) {
			if(query.charAt(i)!=refSeq.charAt(i)) {
				mismatches++;
			}
		}
		mismatches+=Math.abs(query.length()-refSeq.length());
		return mismatches;
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
	
	public GenomicRegion findTandemRepeat(String sequenceName, int first, int last) {
		GenomicRegionImpl region = new GenomicRegionImpl(sequenceName, first, last);
		if(knownSTRs==null) return null;
		List<GenomicRegion> l =knownSTRs.get(sequenceName);
		if(l==null) return null;
		return binaryContains(l, 0, l.size()-1, region);
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
	
	/**
	 * Creates an alignment taking into account that the region overlaps with a tandem repeat
	 * @param aln
	 * @param read
	 * @param qualityScores
	 * @param region
	 * @return
	 */
	public ReadAlignment verifyShortTandemRepeats(String sequenceName, int first, int last, String read, String qualityScores, GenomicRegion region) {
		PairwiseAlignmentWithCigar leftPart = null;
		PairwiseAlignmentWithCigar rightPart = null;
		int firstLeftPart = Math.max(first-10,1);
		int softClipLeft = 0;
		int softClipRight = 0;
		String cigarLeft=null;
		String cigarRight=null;
		char softClipChar = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_SKIPFROMREAD);
		if(first<region.getFirst()-5) {
			CharSequence refSeq = fMIndex.getSequence(sequenceName, firstLeftPart, region.getFirst()-1).toString();
			if(refSeq!=null) {
				int endReadSegment= Math.min(read.length(), region.getFirst()-first+5);
				String readSegment = read.substring(0,endReadSegment);
				//System.out.println(refSeq);
				//System.out.println(readSegment);
				leftPart = new PairwiseAlignmentWithCigar(readSegment, refSeq.toString(), false);
				softClipLeft = readSegment.length()-1 - leftPart.getQueryLastIdx();
				String softClipSubstr = "";
				if(softClipLeft>0) softClipSubstr = ""+softClipLeft+""+softClipChar;
				//System.out.println("Left part query start: "+leftPart.getQueryStartIdx()+" subject start: "+leftPart.getSubjectStartIdx()+" mismatches: "+leftPart.getMismatches()+" cigar: "+leftPart.getCigar()+" soft clip: "+softClipLeft);
				if(leftPart.getMismatches()>3 || firstLeftPart+leftPart.getSubjectLastIdx()!=region.getFirst()-1 || !leftPart.getCigar().endsWith(softClipSubstr)) leftPart = null;
				else {
					//Remove soft clip part and update soft clip total with removed segment of the read
					cigarLeft = leftPart.getCigar().substring(0,leftPart.getCigar().length()-softClipSubstr.length());
					softClipLeft+=(read.length()-readSegment.length());
				}
			}	
		}
		if(last>region.getLast()+5) {
			CharSequence refSeq = fMIndex.getSequence(sequenceName, region.getLast()+1, last+10);
			if(refSeq!=null) {
				int startReadSegment= Math.max(0, read.length()-(last-region.getLast())-5);
				String readSegment = read.substring(startReadSegment);
				//System.out.println(refSeq);
				//System.out.println(readSegment);
				rightPart = new PairwiseAlignmentWithCigar(readSegment, refSeq.toString(), false);
				softClipRight = rightPart.getQueryStartIdx();
				String softClipSubstr = "";
				if(softClipRight>0) softClipSubstr = ""+softClipRight+""+softClipChar;
				//System.out.println("Right part query start: "+rightPart.getQueryStartIdx()+" subject start: "+rightPart.getSubjectStartIdx()+" mismatches: "+rightPart.getMismatches()+" cigar: "+rightPart.getCigar()+" soft clip: "+softClipRight);
				if(rightPart.getMismatches()>3 || rightPart.getSubjectStartIdx()!=0 || !rightPart.getCigar().startsWith(softClipSubstr)) rightPart = null;
				else {
					cigarRight = rightPart.getCigar().substring(softClipSubstr.length());
					softClipRight+=startReadSegment;
				}
			}	
		}
		if(leftPart==null && rightPart ==null) {
			return null;
		}
		if(rightPart==null) {
			//Left alignment with right soft clip
			first = firstLeftPart + leftPart.getSubjectStartIdx();
			last = region.getFirst()-1;
			String cigar = cigarLeft+""+softClipLeft+""+softClipChar;
			double alnQual = 100.0* (read.length() - leftPart.getMismatches())/read.length();
			//System.out.println("Left alignment new genomic coordinates : "+first+"-"+last+" CIGAR:" +cigar+" quality: "+alnQual);
			return buildAln(read, qualityScores, sequenceName, first, last, cigar, alnQual);
		}
		if(leftPart==null) {
			//Right alignment with left soft clip
			first = region.getLast()+1;
			last = first + rightPart.getSubjectLastIdx();
			String cigar = ""+softClipRight+""+softClipChar+cigarRight;
			double alnQual = 100.0* (read.length() - rightPart.getMismatches())/read.length();
			//System.out.println("Right alignment new genomic coordinates : "+first+"-"+last+" CIGAR:" +cigar+" quality: "+alnQual);
			return buildAln(read, qualityScores, sequenceName, first, last, cigar, alnQual);
		}
		first = firstLeftPart + leftPart.getSubjectStartIdx();
		last = region.getLast()+1 + rightPart.getSubjectLastIdx();
		int alignedLeft = read.length()-softClipLeft;
		int alignedRight = read.length()-softClipRight;
		int middleLength = read.length()-alignedLeft-alignedRight;
		if(middleLength<0) return null;
		int difference = region.length()-middleLength;
		//System.out.println("Aligned left: "+alignedLeft+" aligned right: "+alignedRight+" middle length: "+middleLength+" difference: "+difference);
		String cigar = cigarLeft;
		if(difference>0) {
			// Region length > middle length. Add deletion
			cigar+=(difference)+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_DELETION);
			if(middleLength>0) cigar+= (middleLength)+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		} else if (difference<0) {
			// Region length < middle length. Add insertion
			cigar+=(-difference)+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_INSERTION);
			if(middleLength>0) cigar+= (region.length())+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		} else {
			if(middleLength>0) cigar+= (middleLength)+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		}
		
		cigar+=cigarRight;
		double alnQual = 100.0* (read.length() - leftPart.getMismatches()-rightPart.getMismatches())/read.length();
		//System.out.println("Building alignment from first "+first+" last: "+last+" softClipLeft: "+softClipLeft+" softClip right "+softClipRight+" cigar "+cigar);
		
		return buildAln(read, qualityScores, sequenceName, first, last, cigar, alnQual);
	}
}

class PairwiseAlignmentWithCigar {
	private static PairwiseAlignmentAffineGap aligner = new PairwiseAlignmentAffineGap(1, 2, 1, 1);
	private String query;
	private String subject;
	private int subjectStartIdx=-1;
	private int subjectLastIdx=-1;
	private int queryStartIdx=-1;
	private int queryLastIdx=-1;
	private String cigar;
	private int mismatches=0;
	
	public PairwiseAlignmentWithCigar (String query, String subject, boolean includeEndSubject) {
		//System.out.println("Query length: "+query.length()+" subject length: "+subject.length());
		this.query = query;
		this.subject = subject;
		String [] alignment = aligner.getAlignment(query, subject);
		String alnQuery = alignment[0];
		//System.out.println(alnQuery);
		String alnSubject = alignment[1];
		//System.out.println(alnSubject);
		StringBuilder cigarBuilder = new StringBuilder();
		int subjectIdx = 0;
		int queryIdx = 0;
		byte nextOpCode = -1;
		int nextOpLength = 0;
		for(int i=0;i<alnQuery.length();i++) {
			char q = alnQuery.charAt(i);
			char s = alnSubject.charAt(i);
			byte opCode = ReadAlignment.ALIGNMENT_MATCH;
			if(q!=LimitedSequence.GAP_CHARACTER) {
				if(s!=LimitedSequence.GAP_CHARACTER) {
					if(subjectStartIdx==-1) subjectStartIdx = subjectIdx;
					if(queryStartIdx==-1) queryStartIdx = queryIdx;
					queryLastIdx = queryIdx;
					if(q!=s) mismatches++;
					subjectIdx++;
				} else {
					if(subjectIdx==0 || subjectIdx>=subject.length()) {
						opCode = ReadAlignment.ALIGNMENT_SKIPFROMREAD;
					} else {
						opCode = ReadAlignment.ALIGNMENT_INSERTION;
					}
				}
				queryIdx++;
			} else if (s!=LimitedSequence.GAP_CHARACTER) {
				opCode = ReadAlignment.ALIGNMENT_DELETION;
				subjectIdx++;
			}
			if(opCode!=nextOpCode) {
				if(nextOpCode>=0) {
					if(nextOpCode!=ReadAlignment.ALIGNMENT_INSERTION && nextOpCode!=ReadAlignment.ALIGNMENT_SKIPFROMREAD) subjectLastIdx +=nextOpLength;
					if(cigarBuilder.length()>0 || nextOpCode!=ReadAlignment.ALIGNMENT_DELETION || includeEndSubject) {
						cigarBuilder.append(nextOpLength+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(nextOpCode));
						if(nextOpCode!=ReadAlignment.ALIGNMENT_MATCH && nextOpCode!=ReadAlignment.ALIGNMENT_SKIPFROMREAD) mismatches+=2;
						 
					}
				}
				nextOpCode = opCode;
				nextOpLength = 0;
			}
			nextOpLength++;
		}
		if(nextOpCode>=0) {
			if(nextOpCode!=ReadAlignment.ALIGNMENT_DELETION || includeEndSubject) {
				cigarBuilder.append(nextOpLength+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(nextOpCode));
				if(nextOpCode!=ReadAlignment.ALIGNMENT_MATCH && nextOpCode!=ReadAlignment.ALIGNMENT_SKIPFROMREAD) mismatches+=2;
				if(nextOpCode!=ReadAlignment.ALIGNMENT_INSERTION && nextOpCode!=ReadAlignment.ALIGNMENT_SKIPFROMREAD) subjectLastIdx += nextOpLength;
			}
		}
		cigar = cigarBuilder.toString();
	}

	public String getQuery() {
		return query;
	}

	public String getSubject() {
		return subject;
	}


	public int getSubjectStartIdx() {
		return subjectStartIdx;
	}
	
	public int getSubjectLastIdx() {
		return subjectLastIdx;
	}
	
	public int getQueryStartIdx() {
		return queryStartIdx;
	}

	public int getQueryLastIdx() {
		return queryLastIdx;
	}

	public String getCigar() {
		return cigar;
	}

	public int getMismatches() {
		return mismatches;
	}
}





