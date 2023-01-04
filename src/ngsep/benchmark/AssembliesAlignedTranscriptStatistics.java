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
package ngsep.benchmark;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;

/**
 * Script that calculates statistics on transcripts aligned to different genome assemblies.
 * @author Jorge Duitama
 *
 */
public class AssembliesAlignedTranscriptStatistics {

	private Map<String,QualifiedSequence> transcripts = new HashMap<>();
	private Map<String,Integer> refTranscriptCompleteCounts = new HashMap<>();
	private Map<String,Integer> refTranscriptMaxAlnLength = new HashMap<>();
	private int [][] statistics;

	/**
	 * Calculates the following statistics on transcripts aligned to different genome assemblies:
	 * 1. Number of transcripts aligned
	 * 2. Single copy transcripts aligned one time
	 * 3. Single copy transcripts aligned two times
	 * 4. Single copy transcripts aligned more than two times
	 * 5. Total alignments (alleles) of single copy transcripts
	 * 6. Single copy transcripts partially aligned
	 * 7. Single copy transcripts not aligned
	 * 8. Multiple copy transcripts aligned one time
	 * 9. Multiple copy transcripts aligned two times
	 * 10. Multiple copy transcripts aligned more than two times
	 * 11. Total alignments (alleles) of multiple copy transcripts
	 * 12. Multiple copy transcripts partially aligned
	 * 13. Multiple copy transcripts not aligned
	 * The aligned transcripts should be provided in paf format and the file should be sorted by query name.
	 * @param args Files to process. The first is a fasta file with the transcripts.
	 * The remaining arguments should be paf files.
	 * To determine if a gene is unique, it is assumed that the first paf file corresponds to alignments against an unphased reference genome 
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		AssembliesAlignedTranscriptStatistics instance = new AssembliesAlignedTranscriptStatistics();
		String transcriptsFastaFile = args[0];
		String referenceAlnsPafFile = args[1];
		instance.initialize(transcriptsFastaFile,args.length-1);
		instance.processAlignments(referenceAlnsPafFile, 0);
		for(int i=2;i<args.length;i++) {
			instance.processAlignments(args[i],i-1);
		}
		instance.printResults(System.out);

	}
	private void initialize(String transcriptsFastaFile, int testAssemblies) throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> transcriptsList = handler.loadSequences(transcriptsFastaFile);
		for(QualifiedSequence transcript:transcriptsList) {
			refTranscriptCompleteCounts.put(transcript.getName(), 0);
			refTranscriptMaxAlnLength.put(transcript.getName(), 0);
			transcripts.put(transcript.getName(), transcript);
		}
		statistics = new int [13][testAssemblies];
	}
	private void processAlignments(String alignmentsFile, int datasetId) throws IOException {
		try (FileReader reader = new FileReader(alignmentsFile);
			 BufferedReader in = new BufferedReader(reader)) {
			String line = in.readLine();
			String nextTranscript = null;
			int nextCompleteCount = 0;
			int nextPartialCount = 0;
			Set<String> transcriptsMapped = new HashSet<>();
			for (int i = 0; line != null; i++) {
				String[] items = line.split("\t");
				String transcriptId = items[0];
				QualifiedSequence transcript = transcripts.get(transcriptId);
				if(transcript==null) {
					System.err.println("WARN: Transcript with id: "+transcriptId+" not found");
					line = in.readLine();
					continue;
				}
				int tlen = transcript.getLength();
				if(tlen != Integer.parseInt(items[1])) {
					System.err.println("WARN: Inconsistent length for transcript with id: "+transcriptId+" given: "+Integer.parseInt(items[1])+" real: "+tlen);
					line = in.readLine();
					continue;
				}
				int queryFirst = Integer.parseInt(items[2]);
				int queryLast = Integer.parseInt(items[3]);
				String subjectName = items[5];
				int subjectFirst = Integer.parseInt(items[7]);
				int subjectLast = Integer.parseInt(items[8]);
				int alnLen = queryLast-queryFirst;
				boolean complete = alnLen>0.9*tlen;
				boolean partial = !complete && alnLen>0.3*tlen;
				if(datasetId==0) {
					refTranscriptMaxAlnLength.compute(transcriptId,(k,v)->Math.max(alnLen, v));
					if (complete) refTranscriptCompleteCounts.compute(transcriptId,(k,v)->v+1);
				}
				if(!transcriptId.equals(nextTranscript)) {
					if(nextTranscript!=null) processCount(datasetId, nextTranscript,nextCompleteCount,nextPartialCount);
					nextTranscript = transcriptId;
					nextCompleteCount = 0;
					nextPartialCount = 0;	
				}
				if(complete) nextCompleteCount++;
				else if(partial) nextPartialCount++;
				if(complete || partial) transcriptsMapped.add(transcriptId);
				line = in.readLine();
			}
			if(nextTranscript!=null) processCount(datasetId, nextTranscript,nextCompleteCount,nextPartialCount);
			statistics[0][datasetId]=transcriptsMapped.size();
			//Counts for transcripts not present in the alignments
			for(String tid:transcripts.keySet()) {
				if(!transcriptsMapped.contains(tid)) {
					int refCount = calculateReferenceCount(tid);
					if(refCount<2) statistics[6][datasetId]++;
					else statistics[12][datasetId]++;
					
				}
			}
			statistics[0][datasetId]=transcriptsMapped.size();
		}
		
	}
	private void processCount(int datasetId, String tid, int completeCount, int partialCount) {
		int refCount = calculateReferenceCount(tid);
		if(refCount==0) System.err.println("WARN: Transcript "+tid+" found in dataset: "+datasetId+" and not found in reference alignments");
		int offsetRow = (refCount <2)?1:7;
		if(completeCount==0) {
			if (partialCount>0) statistics[offsetRow+4][datasetId]++;
			else statistics[offsetRow+5][datasetId]++;
		} else {
			statistics[offsetRow+3][datasetId]+=completeCount;
			if (completeCount == 1) {
				//if(refCount>1) System.out.println("Multigene: "+tid+" dataset: "+datasetId+" complete copies: "+completeCount+" partial copies: "+partialCount);
				statistics[offsetRow][datasetId]++;
			}
			else if (completeCount == 2) statistics[offsetRow+1][datasetId]++;
			else {
				//if(refCount>1) System.out.println("Multigene: "+tid+" complete copies: "+completeCount+" partial copies: "+partialCount);
				statistics[offsetRow+2][datasetId]++;
			}
		}
	}
	public int calculateReferenceCount(String tid) {
		int refCount = refTranscriptCompleteCounts.get(tid);
		if(refCount==0) {
			//If there was at least one partial alignment
			if(refTranscriptMaxAlnLength.get(tid)>0) refCount=1;
		}
		return refCount;
	}
	private void printResults(PrintStream out) {
		out.println("Total transcripts: "+transcripts.size());
		String [] titles = {"Aligned one time", "Aligned two times", "Aligned three or more times", "Total copies", "Partially aligned", "Not aligned"};
		int k=0;
		for(int i=0;i<statistics.length;i++) {
			if(i==1) out.println("Single copy");
			if(i==7) {
				k=0;
				out.println("Multiple copies");
			}
			if(i==0) out.print("Aligned");
			else {
				out.print(titles[k]);
				k++;
			}
			for(int j=0;j<statistics[i].length;j++) {
				out.print("\t"+statistics[i][j]);
			}
			out.println();
		}	
	}

}
