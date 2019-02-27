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
package ngsep.transcriptome;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.math.Distribution;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;

/**
 * @author Tatiana Garcia
 * @author Jorge Duitama
 */
public class TranscriptomeAnalyzer {

	private ReferenceGenome genome;
	
	private boolean selectCompleteProteins = false;
	
	private Logger log = Logger.getLogger(TranscriptomeAnalyzer.class.getName());
	
	private ProteinTranslator translator=  new ProteinTranslator();
	
	public static void main(String[] args) throws Exception {
		TranscriptomeAnalyzer instance = new TranscriptomeAnalyzer();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String transcriptomeFile = args[i++];
		String outPrefix = args[i++];
		instance.processTranscriptome(transcriptomeFile, outPrefix);
	}

	/**
	 * @return the log
	 */
	public Logger getLog() {
		return log;
	}

	/**
	 * @param log the log to set
	 */
	public void setLog(Logger log) {
		this.log = log;
	}


	/**
	 * @return the genome
	 */
	public ReferenceGenome getGenome() {
		return genome;
	}

	/**
	 * @param genome the genome to set
	 */
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	
	

	/**
	 * @return the selectCompleteProteins
	 */
	public boolean isSelectCompleteProteins() {
		return selectCompleteProteins;
	}

	/**
	 * @param selectCompleteProteins the selectCompleteProteins to set
	 */
	public void setSelectCompleteProteins(boolean selectCompleteProteins) {
		this.selectCompleteProteins = selectCompleteProteins;
	}
	public void setSelectCompleteProteins(Boolean selectCompleteProteins) {
		this.setSelectCompleteProteins(selectCompleteProteins.booleanValue());
	}

	public void processTranscriptome(String transcriptomeFile, String outPrefix) throws IOException {
		Distribution geneLengthDist = new Distribution (0,10000,200);
		Distribution transcriptLengthDist = new Distribution (0,5000,200); 
		Distribution numberOfExonsDist = new Distribution (0,50,1);
		Distribution proteinLengthDist = new Distribution (0,2500,200);
		Distribution transcriptsPerGeneDist = new Distribution (0,20,1);
		
		//Load transcriptome
		GFF3TranscriptomeHandler gff3Handler = new GFF3TranscriptomeHandler();
		gff3Handler.setLog(log);
		QualifiedSequenceList sequenceNames;
		if(genome!=null) {
			sequenceNames = genome.getSequencesMetadata();
			gff3Handler.setSequenceNames(sequenceNames);
		}
		Transcriptome transcriptome = gff3Handler.loadMap(transcriptomeFile); 			
		if(genome!=null) transcriptome.fillSequenceTranscripts(genome);
		else sequenceNames = gff3Handler.getSequenceNames();
		
		List<Transcript> transcriptsList = transcriptome.getAllTranscripts();
		Set<String> visitedGeneIDs = new HashSet<>();
		QualifiedSequenceList proteome = new QualifiedSequenceList();
		
		for(Transcript t: transcriptsList) {
			//Collect transcript structure statistics
			if (!t.isCoding()) continue;
			int start = t.getCodingRelativeStart();
			int end = t.getCodingRelativeEnd();
			if(t.getCodingRelativeEnd()<4) continue;
			transcriptLengthDist.processDatapoint(t.length());
			List<TranscriptSegment> exons=t.getTranscriptSegments();
			numberOfExonsDist.processDatapoint(exons.size());		
			String geneId= t.getGeneId(); 

			if(!visitedGeneIDs.contains(geneId)) {
				//Collect gene statistics
				List<Transcript> geneTranscripts = transcriptome.getTranscriptsByGene(geneId);
				transcriptsPerGeneDist.processDatapoint(geneTranscripts.size());
				int minFirst = -1;
				int maxLast = -1; 
				for (Transcript tr: geneTranscripts ) {
					int first = tr.getFirst();
					int last = tr.getLast();
					if(minFirst==-1 || minFirst>first) minFirst = first;
					if(maxLast==-1 || maxLast<last) maxLast = last;	
				}
				if(minFirst>=0 && maxLast>=0) geneLengthDist.processDatapoint(maxLast-minFirst+1);
				
				visitedGeneIDs.add(geneId);
			}
			String protein= t.getProteinSequence(translator);
			if(protein==null) continue;
			
			//Collect sequence statistics
			proteinLengthDist.processDatapoint(protein.length());
			String cdnaSequence = t.getCDNASequence().toString();
			Codon startCodon = translator.getCodon(cdnaSequence.charAt(start), cdnaSequence.charAt(start+1), cdnaSequence.charAt(start+2));
			if(startCodon==null) {
				log.info("Transcript "+t.getId()+" has an invalid start codon. Sequence: "+cdnaSequence.substring(start, start+3));
				continue;
			}
			else if(!startCodon.isStart()) log.info("Transcript "+t.getId()+" does not have a standard start codon. Codon: "+startCodon.getRnaSequence());
			
			Codon stopCodon = translator.getCodon(cdnaSequence.charAt(end-2), cdnaSequence.charAt(end-1), cdnaSequence.charAt(end));
			if(stopCodon==null) log.info("Transcript "+t.getId()+" has an invalid stop codon. Sequence: "+cdnaSequence.substring(end-2, end+1));
			else if(!stopCodon.isStop()) log.info("Transcript "+t.getId()+" does not have a standard stop codon. Codon: "+stopCodon.getRnaSequence());
			
			if(!selectCompleteProteins || (startCodon.isStart() && stopCodon!=null && stopCodon.isStop())) {
				QualifiedSequence qp = new QualifiedSequence(t.getId(),protein);
				qp.setComments(t.getGeneId()+" "+t.getGeneName());
				proteome.add(qp);
			}
		}
		try (PrintStream out = new PrintStream(outPrefix+"_stats.txt")) {
			out.println("Gene length"); 
			geneLengthDist.printDistributionInt(out);
			out.println();
			out.println("Transcript length");
			transcriptLengthDist.printDistributionInt(out);
			out.println();
			out.println("Number of exons");
			numberOfExonsDist.printDistributionInt(out);
			out.println();
			out.println("Transcripts per gene");
			transcriptsPerGeneDist.printDistributionInt(out);
			if(proteome.size()>0) {
				out.println();
				out.println("Protein length");
				proteinLengthDist.printDistributionInt(out);
			}
		}
		if(proteome.size()>0) {
			FastaSequencesHandler handler = new FastaSequencesHandler();
			try (PrintStream out = new PrintStream(outPrefix+"_proteins.fa")) {
				handler.saveSequences(proteome, out, 100);
			}
		}
	}
}
