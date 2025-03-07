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
import java.util.List;
import java.util.logging.Logger;

import ngsep.alignments.PairwiseAlignerSimpleGap;
import ngsep.alignments.PairwiseAlignment;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.ReadsAligner;
import ngsep.alignments.ReadAlignment.Platform;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;

/**
 * 
 * @author Erick Duarte
 *
 */
public class CircularSequencesProcessor {
	private Logger log = Logger.getLogger(CircularSequencesProcessor.class.getName());
	private ProgressNotifier progressNotifier = null;
	//constant values for default
	public static final int DEF_MAX_LENGTH = 6000000;
	private static final int END_LENGTH = 1000;
	private static final int THRESHOLD=50;
	
	// Parameters
	private String inputFile;
	private String outputFile;
	private List<QualifiedSequence> starts; // Start sequences for circularization
	private int maxLength = DEF_MAX_LENGTH;
	
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

	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	
	public List<QualifiedSequence> getStarts() {
		return starts;
	}
	public void setStarts(List<QualifiedSequence> starts) {
		this.starts = starts;
	}
	public void setStarts(String filename) throws IOException {
		if(filename==null) return; 
		FastaSequencesHandler handler = new FastaSequencesHandler();
		starts = handler.loadSequences(filename);
	}
	
	public int getMaxLength() {
		return maxLength;
	}
	public void setMaxLength(int maxLength) {
		this.maxLength = maxLength;
	}
	public void setMaxLength(String value) {
		this.setMaxLength((int) OptionValuesDecoder.decode(value, Integer.class));
	}

	public static void main(String[] args) throws Exception {
		CircularSequencesProcessor instance = new CircularSequencesProcessor();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	public void run() throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> contigs = handler.loadSequences(inputFile);
		processContigs(contigs);
		try (PrintStream out = new PrintStream(outputFile)) {
			handler.saveSequences(contigs, out, 100);
		}
	}
	
	/**
	 * It receives a list of contigs to perform a processContig method in each of them
	 * @param List of contigs
	 */
	public void processContigs(List<QualifiedSequence> contigs) {
		for (QualifiedSequence contig:contigs) {
			processContig(contig);
			
		}
	}
	/**
	 * It process the contig in order to find a duplicated fragment at the beginning and end and then remove it
	 * @param contig
	 */
	public void processContig(QualifiedSequence contig) {
		if(contig.getLength()>maxLength) return;
		removeRepeatedEnds(contig);
		if(starts!=null) {
			findStart(contig);
		}
			
	}

	private void removeRepeatedEnds(QualifiedSequence contig) {
		String sequence = contig.getCharacters().toString();
		log.info("Next sequence: "+contig.getName()+" length: "+sequence.length());
		int maximumInterval = Math.min(10000, sequence.length()/2);
		PairwiseAlignerSimpleGap aligner = new PairwiseAlignerSimpleGap();
		//PairwiseAlignerDynamicKmers aligner = new PairwiseAlignerDynamicKmers();
		aligner.setForceStart1(false);
		aligner.setForceEnd2(false);
		for (int i=END_LENGTH; i<maximumInterval;i+=END_LENGTH) {
			String subseqini= sequence.substring(0, i);
			String subseqend= sequence.substring(sequence.length()-i);
			//System.err.println("subseq ini: "+subseqini);
			//System.err.println("subseq end: "+subseqend);
			PairwiseAlignment alignment = aligner.calculateAlignment(subseqini,subseqend);
			String[] seqs = {alignment.getAlignedSequence1(),alignment.getAlignedSequence2()}; 
			String [] matchingZone = fixAlignment(seqs);
			log.info("End length: "+i+" Next match: "+ matchingZone[0].length());
			if (passFilterMethod(matchingZone[0], matchingZone[1],i)) {
				int duplicatedLength=calculateDuplicatedLength(subseqini, matchingZone[0]);
				String modifiedSequence = removeDuplicatedSubSeq(sequence, duplicatedLength);
				log.info("Sequence "+contig.getName()+" trimmed basepairs "+duplicatedLength+" new length: "+modifiedSequence.length());
				System.out.println(modifiedSequence.length());
				contig.setCharacters(modifiedSequence);
				break;
			}
		}
	}

	/**
	 * it receives a matrix of the aligned fragments with gaps
	 * @param matrix of the alignment between both fragments with not forced start and end 
	 * @return a fixed matrix of the alignment between both fragments with no gaps
	 */
	private String [] fixAlignment(String[] algn) {
		String[] matrix = new String[2];
		String subseqalgn1= algn[0];
		String subseqalgn2=algn[1];
		char empty = '-';
		int p1=0;
		int p2=0;
		for (int i=subseqalgn1.length()-1;i>=0;i--) {
			char nucleotid1 = subseqalgn1.charAt(i);
			if(nucleotid1 != empty) {
				p2=i;
				//System.out.println(p2);
				break;
			}
				
		}	


			//System.out.println(matrix[0]);
		for (int j=0;j<subseqalgn2.length() && j<p2;j++) {
			char nucleotid2 = subseqalgn2.charAt(j);
			if(nucleotid2 != empty) {
				p1=j;
				//System.out.println(p1);
				break;
			}
		}
		// sacar substring de p1 a p2 que es el match que nos importa
		matrix[0] = subseqalgn1.substring(p1,p2);
		//System.out.println("MATRIZ 1 " +matrix[0]);
		matrix[1] = subseqalgn2.substring(p1, p2);
		//System.out.println("MATRIZ 2 " + matrix[1]);
		return matrix;
	}
	
	/**
	 * Determine which fragments have a good alignment according the THRESHOLD and HammingDistance
	 * @param subseq1 first fragment aligned
	 * @param subseq2 second fragment aligned
	 * @param totalLength length of the query sequence
	 * @return boolean indicating if the alignment meets the requirements
	 */
	public boolean passFilterMethod(String subseq1, String subseq2, int totalLength) {
		double lengthAlignment = subseq1.length();
		double calc = (lengthAlignment/totalLength)*100;
		//System.out.println("CALC " + calc);
		boolean passFilter = calc >THRESHOLD;
		double distance = (new HammingSequenceDistanceMeasure()).calculateNormalizedDistance(subseq1, subseq2);
		passFilter = passFilter && distance< 0.1;
		return passFilter; 
	}
	
	
	/**
	 * It calculates the length of the real sequence that it is duplicated according to the length of match
	 * @param searchingSeq sequence used as a query to look for a duplicated fragment at the end of the contig
	 * @param matchingSeq fragment of the query sequence that has made match with a fragment of the end-contig sequence
	 * @return the length of the sequence duplicated at the beginning and end of the contig
	 */
	public int  calculateDuplicatedLength (String searchingSeq, String matchingSeq) {
		int lengthSeachingSeq= searchingSeq.length();
		int lengthMatchingZone=matchingSeq.length();
		int duplicatedLength = lengthSeachingSeq + (lengthSeachingSeq-lengthMatchingZone);
		return duplicatedLength;
	}
	
	/**
	 * It removes the duplicated sequence of the contig
	 * @param String contig
	 * @param int length of the duplicated sequence
	 * @return the contig with removed duplication
	 */
	public String removeDuplicatedSubSeq (String contig, int duplicatedLength) {
		String substringRemove = contig.substring(0,contig.length()-duplicatedLength);
		return substringRemove;	
	}
	private void findStart(QualifiedSequence contig) {
		ReferenceGenome genome = new ReferenceGenome(contig);
		ReadsAligner aligner = new ReadsAligner(genome, Platform.PACBIO);
		aligner.setLog(log);
		List<ReadAlignment> alns = new ArrayList<>();
		for(QualifiedSequence startSeq:starts) {
			List<ReadAlignment> alnsRead = aligner.alignRead(startSeq);
			//System.err.println("Start seq: "+startSeq.getName()+" Alignments: "+alnsRead);
			alns.addAll(alnsRead);
		}
		if(alns.size()==0) return;
		Collections.sort(alns,(a1,a2)->a2.getAlignmentQuality()- a1.getAlignmentQuality());
		ReadAlignment aln = alns.get(0);
		int cutSite = Math.max(0, aln.getFirst()-aln.getSoftClipStart()-1);
		if(aln.isNegativeStrand()) cutSite = Math.min(contig.getLength(), aln.getLast()+aln.getSoftClipEnd());
		String sequence = contig.getCharacters().toString();
		String subseqLeft = sequence.substring(0,cutSite);
		String subseqRight = sequence.substring(cutSite);
		String reorganizedSequence = subseqRight+subseqLeft;
		if(aln.isNegativeStrand()) {
			reorganizedSequence = DNAMaskedSequence.getReverseComplement(reorganizedSequence).toString(); 
		}
		log.info("Sequence "+contig.getName()+" start "+cutSite+" negative strand: "+aln.isNegativeStrand());
		contig.setCharacters(reorganizedSequence);
	}		
}