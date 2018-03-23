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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import ngsep.main.CommandsDescriptor;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.KmersCounter;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;

/**
 * Program to align reads to a refernece genome
 * @author German Andrade
 * @author Jorge Duitama
 */
public class ReadsAligner {
	
	static final int SEARCH_KMER_LENGTH = 15;
	static final double MIN_ACCURACY =0.5;

	public static void main(String[] args) throws Exception 
	{
		ReadsAligner instance = new ReadsAligner();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String fMIndexFile = args[i++];
		String readsFile = args[i++];
		String outFile = args[i++];
		try (PrintStream out = new PrintStream(outFile)){
			instance.alignReads(fMIndexFile, readsFile, out);
		}
	}
	/**
	 * Aligns readsFile with the fMIndexFile
	 * @param fMIndexFile Binary file with the serialization of an FMIndex
	 * @param readsFile Fastq file with the reads to align
	 * @param out
	 * @throws IOException
	 */
	public void alignReads( String fMIndexFile, String readsFile, PrintStream out) throws IOException {
		FMIndex fMIndex = FMIndex.loadFromBinaries(fMIndexFile);
		int totalReads = 0;
		int readsAligned = 0;
		int uniqueAlignments=0;
		long time = System.currentTimeMillis();
		try (FastqFileReader reader = new FastqFileReader(readsFile)) {
			Iterator<RawRead> it = reader.iterator();
			while(it.hasNext()) {
				RawRead read = it.next();
				int numAlns = alignRead(fMIndex, read, out);
				totalReads++;
				if(numAlns>0) readsAligned++;
				if(numAlns==1) uniqueAlignments++;
			}
		}
		System.out.println("Total reads: "+totalReads);
		System.out.println("Reads aligned: "+readsAligned);
		System.out.println("Unique alignments: "+uniqueAlignments);
		System.out.println("Overall alignment rate: "+(100.0*readsAligned/(double)totalReads)+"%");
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		System.out.println("Time: "+seconds+" seconds");
		
	}
	
	

	private int alignRead(FMIndex fMIndex, RawRead read, PrintStream out) {
		List<ReadAlignment> alignments = search(fMIndex, read);
		int i=0;
		for (ReadAlignment aln: alignments) {
			String readSeq = read.getSequenceString();
			String qual = read.getQualityScores();
			if(aln.isNegativeStrand()) {
				readSeq = DNAMaskedSequence.getReverseComplement(readSeq).toString();
				qual = new StringBuilder(qual).reverse().toString();
			}
			if(i>0) aln.setFlags(aln.getFlags()+ReadAlignment.FLAG_SECONDARY);
			out.println(
					//1.query name
					read.getName()+"\t"+
							
					//2.Flag
					aln.getFlags()+"\t"+
					
					//3.reference sequence name
					aln.getSequenceName()+"\t"+
					
					//4.POS
					aln.getFirst()+"\t"+
					
					//5.MAPQ
					"255\t"+
					
					//6.CIGAR
					read.getLength()+"M\t"+
					
					//7. RNEXT
					"*\t"+
					
					//8. PNEXT
					"0\t"+
					
					//9. TLEN
					"0\t"+
					
					//10. SEQ
					readSeq+"\t"+
					
					//11. QUAL
					qual
					
					);
			i++;
		}
		return alignments.size();
	}
	
	public List<ReadAlignment> search (FMIndex fMIndex, RawRead read) {
		return kmerBasedInexactSearchAlgorithm(fMIndex, read);
	}
	
	public List<ReadAlignment> exactSearch (FMIndex fMIndex, RawRead read) {
		return fMIndex.search(read.getSequenceString());
	}
	/**
	 * 
	 * @param fMIndex
	 * @param read 
	 * @return 
	 */
	private List<ReadAlignment> kmerBasedInexactSearchAlgorithm (FMIndex fMIndex, RawRead read) 
	{
		//Find Kmers
		CharSequence[] kmers = KmersCounter.extractKmers(read.getCharacters().toString(), SEARCH_KMER_LENGTH, true);
		
		//Stores the aligments of each sequence of the genome, the alignment store their kmerNumber in the 
		//MateFirst field
		HashMap<String,List<ReadAlignment>> seqHits =  new HashMap<String,List<ReadAlignment>>();
		
		//Avoid overlaps
		int kmersCount=0;
		for (int i = 0; i < kmers.length; i+=SEARCH_KMER_LENGTH) 
		{
			//Exit loop if kmers[i] is null
			if(kmers[i]==null)
				continue;
			
			String kmer =kmers[i].toString();
			
			//Where is located the kmer in exact way
			List<ReadAlignment> regions=fMIndex.search(kmer);
			
			for(ReadAlignment aln:regions)
			{
				//Use mate start to store the kmer start site producing the hit
				aln.setMateFirst(i);
				
				List<ReadAlignment> seqAlns = seqHits.get(aln.getSequenceName());

				if(seqAlns==null) {
					seqAlns = new ArrayList<>();
					seqHits.put(aln.getSequenceName(), seqAlns);
				}
				seqAlns.add(aln);
			}
			kmersCount++;
		}
		
		//Processing part
		
		List<ReadAlignment> finalAlignments =  new ArrayList<>();
		ReadAlignmentMateFirstComparator cmp = ReadAlignmentMateFirstComparator.getInstance();
		
		for (String sequenceName: seqHits.keySet())
		{
			Set<Integer> kmersInSequence= new HashSet<>();
			List<ReadAlignment> alns = seqHits.get(sequenceName);
			Collections.sort(alns,cmp);
			for (int i = 0; i < alns.size(); i++) {
				kmersInSequence.add(alns.get(i).getMateFirst());
			}
			if(kmersInSequence.size()/kmersCount>=MIN_ACCURACY)
			{
				ReadAlignment first = alns.get(0);
				ReadAlignment last = alns.get(alns.size()-1);
				ReadAlignment readAlignment =new ReadAlignment(first.getSequenceName(), first.getFirst(), 
						last.getLast(), last.getLast()-first.getFirst(), first.getFlags());
				finalAlignments.add(readAlignment);
			}
			
		}
		
		
		//TODO: Implement
		return finalAlignments;
	}
}
