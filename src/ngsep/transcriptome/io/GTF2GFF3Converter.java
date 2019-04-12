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
package ngsep.transcriptome.io;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.ProteinTranslator;
import ngsep.transcriptome.Transcript;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class GTF2GFF3Converter {

	private ReferenceGenome genome;
	public static void main(String[] args) throws Exception {
		GTF2GFF3Converter instance = new GTF2GFF3Converter();
		instance.genome = new ReferenceGenome(args[0]);
		instance.run(args[1], System.out);

	}

	public void run(String gtfFile, PrintStream out) throws IOException {
		QualifiedSequenceList seqNames = genome.getSequencesMetadata();
		Map<String,Gene> genes = new HashMap<>();
		try (FileReader reader = new FileReader(gtfFile);
			 BufferedReader in = new BufferedReader(reader)) {
			Transcript transcript = null;
			List<GenomicRegion> rawExons = new ArrayList<>();
			out.println("##gff-version	3");
			String line = in.readLine();
			while (line != null) {
				if(line.charAt(0)=='#') {
					//Skip comments
					line = in.readLine();
					continue;
				}
				String [] items = line.split("\t|;");
				QualifiedSequence seq = seqNames.get(items[0]);
				int first = Integer.parseInt(items[3]);
				int last = Integer.parseInt(items[4]);
				if(items[6].charAt(0)=='.') {
					//Ignoring feature without direction
					line = in.readLine();
					continue;
				}
				boolean negativeStrand = items[6].charAt(0)=='-';
				String type = items[2];
				Map<String,String> attributes = loadAttributes (items);
				String geneId = attributes.get("gene_id");
				if(geneId == null) {
					System.err.println("Gene id not found at line: "+line);
					line = in.readLine();
					continue;
				}
				String transcriptId = attributes.get("transcript_id");
				if(transcriptId == null) {
					System.err.println("Transcript id not found at line: "+line);
					line = in.readLine();
					continue;
				}
				if("transcript".equals(type)) {
					processTranscript (transcript,rawExons,out);
					Gene gene = genes.get(geneId);
					if(gene == null) {
						//Write gene
						gene = new Gene(geneId, geneId);
						genes.put(geneId, gene);
						out.println(items[0]+"\tSimpleGTF2GFF\tgene\t"+first+"\t"+last+"\t.\t"+items[6]+"\t.\tID="+geneId);
					}
					transcript = new Transcript(transcriptId, seq.getName(), first, last, negativeStrand);
					transcript.setGene(gene);
					rawExons.clear();
				} else if ("exon".equals(type)) {
					if(!transcript.getId().equals(transcriptId)) {
						System.err.println("Disorganized exon at line: "+line);
						line = in.readLine();
						continue;
					}
					if(transcript.getFirst()>first || transcript.getLast()<last) {
						System.err.println("Exon coordinates outside transcript coordinates. Line: "+line);
						line = in.readLine();
						continue;
					}
					rawExons.add(new GenomicRegionImpl(seq.getName(), first, last));
				}
				
				line = in.readLine();
			}
			processTranscript (transcript,rawExons,out);
		}
		
	}

	private Map<String, String> loadAttributes(String[] items) {
		Map<String, String> attributes = new HashMap<>();
		for(int i=8;i<items.length;i++) {
			String nextPair = items[i].trim();
			int j = nextPair.indexOf(' ');
			//System.out.println("Next pair: "+nextPair+". space chr: "+j);
			if(j>0 && j<nextPair.length()-1 && nextPair.charAt(j+1)=='"') {
				attributes.put(nextPair.substring(0,j), nextPair.substring(j+2,nextPair.length()-1));
			}
		}
		return attributes;
	}
	
	private void processTranscript(Transcript transcript, List<GenomicRegion> rawExons, PrintStream out) {
		if(transcript==null || rawExons.size()==0) return; 
		StringBuilder mRNAB = new StringBuilder();
		Collections.sort(rawExons,GenomicRegionPositionComparator.getInstance());
		List<GenomicRegion> exonsByTranscript = new ArrayList<>(rawExons);
		for(GenomicRegion exon:rawExons) {
			CharSequence nextSeq = genome.getReference(exon);
			if(nextSeq==null) {
				System.err.println("Sequence for exon at "+exon.getSequenceName()+":"+exon.getFirst()+"-"+exon.getLast()+" could not be retrieved");
				return;
			}
			mRNAB.append(nextSeq);
		}
		String mRNA = mRNAB.toString();
		if(transcript.isNegativeStrand()) {
			mRNA = DNAMaskedSequence.getReverseComplement(mRNA);
			Collections.reverse(exonsByTranscript);
		}
		int maxL = -1;
		int startLongest = -1;
		ProteinTranslator translator = new ProteinTranslator();
		char [] mRNACharArray = mRNA.toCharArray();
		for(int i=0;i<mRNA.length();i++) {
			char [] nextP = translator.getProteinSequence(mRNACharArray,i,false);
			if(nextP.length>0 && nextP[0]=='M') {
				if(maxL < nextP.length) {
					maxL = nextP.length;
					startLongest = i;
				}
			}
		}
		if(startLongest>=0) {
			String negStr = transcript.isNegativeStrand()?"-":"+";
			//Write mRNA feature
			out.println(transcript.getSequenceName()+"\tSimpleGTF2GFF\tmRNA\t"+transcript.getFirst()+"\t"+transcript.getLast()+"\t.\t"+negStr+"\t.\tID="+transcript.getId()+";Parent="+transcript.getGeneId());
			int exId = 1;
			for(GenomicRegion exon:rawExons) {
				//Write exon feature
				out.println(exon.getSequenceName()+"\tSimpleGTF2GFF\texon\t"+exon.getFirst()+"\t"+exon.getLast()+"\t.\t"+negStr+"\t.\tID="+transcript.getId()+"_e"+exId+";Parent="+transcript.getId());
				exId++;
			}
			
			
			char [] protein = translator.getProteinSequence(mRNACharArray,startLongest,false);
			int start3PUTR = startLongest+3*(protein.length+1);
			//System.out.println("mRNA len: "+mRNA.length()+" Start coding: "+startLongest+" start 3PUTR: "+start3PUTR+" length: "+protein.length+" protein: "+new String(protein));
			//System.out.println("mRNA: "+mRNA);
			int exonRelativeFirst=0;
			int module = 0;
			for(GenomicRegion exon:exonsByTranscript) {
				int exonRelativeLast = exonRelativeFirst+exon.length()-1;
				if(exonRelativeLast<startLongest) {
					//Five prime UTR
					out.println(exon.getSequenceName()+"\tSimpleGTF2GFF\tfive_prime_UTR\t"+exon.getFirst()+"\t"+exon.getLast()+"\t.\t"+negStr+"\t.\tParent="+transcript.getId());
				} else if (exonRelativeFirst>=start3PUTR) {
					//Three prime UTR
					out.println(exon.getSequenceName()+"\tSimpleGTF2GFF\tthree_prime_UTR\t"+exon.getFirst()+"\t"+exon.getLast()+"\t.\t"+negStr+"\t.\tParent="+transcript.getId());
				} else if (startLongest<=exonRelativeFirst && exonRelativeLast<start3PUTR) {
					//CDS
					out.println(exon.getSequenceName()+"\tSimpleGTF2GFF\tCDS\t"+exon.getFirst()+"\t"+exon.getLast()+"\t.\t"+negStr+"\t"+getPhase(module)+"\tParent="+transcript.getId());
					int cdsLenPlusModule = module + exon.getLast() - exon.getFirst() + 1;
					module = cdsLenPlusModule % 3;
				} else {
					int cdsFirst = exon.getFirst();
					int cdsLast = exon.getLast();
					if (exonRelativeFirst<startLongest) {
						int diff = startLongest-exonRelativeFirst;
						if(transcript.isPositiveStrand()) {
							cdsFirst +=diff;
							out.println(exon.getSequenceName()+"\tSimpleGTF2GFF\tfive_prime_UTR\t"+exon.getFirst()+"\t"+(cdsFirst-1)+"\t.\t+\t.\tParent="+transcript.getId());
						} else {
							cdsLast -=diff;
							out.println(exon.getSequenceName()+"\tSimpleGTF2GFF\tfive_prime_UTR\t"+(cdsLast+1)+"\t"+exon.getLast()+"\t.\t-\t.\tParent="+transcript.getId());
						}
					}
					int threePUTRFirst = 0;
					int threePUTRLast = 0;
					if (exonRelativeLast>=start3PUTR) {
						int diff = exonRelativeLast-start3PUTR+1;
						if(transcript.isPositiveStrand()) {
							cdsLast -=diff;
							threePUTRFirst = cdsLast+1;
							threePUTRLast = exon.getLast();
						} else {
							cdsFirst +=diff;
							threePUTRFirst = exon.getFirst();
							threePUTRLast = cdsFirst-1;
						}
					}
					out.println(exon.getSequenceName()+"\tSimpleGTF2GFF\tCDS\t"+cdsFirst+"\t"+cdsLast+"\t.\t"+negStr+"\t"+getPhase(module)+"\tParent="+transcript.getId());
					int cdsLenPlusModule = module + cdsLast - cdsFirst + 1;
					module = cdsLenPlusModule % 3;
					if(threePUTRFirst>0) {
						out.println(exon.getSequenceName()+"\tSimpleGTF2GFF\tthree_prime_UTR\t"+threePUTRFirst+"\t"+threePUTRLast+"\t.\t"+negStr+"\t.\tParent="+transcript.getId());
					}
				}
				exonRelativeFirst=exonRelativeLast+1;
			}
		} else {
			System.err.println("No protein found for "+transcript.getSequenceName()+":"+transcript.getFirst()+"-"+transcript.getLast()+" "+transcript.getId());
			return;
		}
		
	}

	private int getPhase(int module) {
		int phase = 0;
		if(module == 1) phase =2;
		if(module == 2) phase =1;
		return phase;
	}

	

}
