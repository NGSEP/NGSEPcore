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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionCollectionDensityCalculator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
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

	// Constants for default values
	
	// Logging and progress
	private Logger log = Logger.getLogger(TranscriptomeAnalyzer.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private String outputPrefix = null;
	private ReferenceGenome genome;
	
	// Model attributes
	
	private ProteinTranslator translator=  new ProteinTranslator();
	private GenomicRegionCollectionDensityCalculator densityCalculator = new GenomicRegionCollectionDensityCalculator();
	
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
	
	public String getOutputPrefix() {
		return outputPrefix;
	}
	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}

	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
		densityCalculator.setGenome(genome);
	}
	public void setGenome(String genomeFile) throws IOException {
		setGenome(OptionValuesDecoder.loadGenome(genomeFile,log));
	}
	
	public static void main(String[] args) throws Exception {
		TranscriptomeAnalyzer instance = new TranscriptomeAnalyzer();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	
	public void run () throws IOException {
		if(inputFile==null) throw new IOException("The input transcriptome file is a required parameter");
		if (genome == null) throw new IOException("The file with the reference genome is a required parameter");
		if(outputPrefix==null) throw new IOException("The prefix of the output files is a required parameter");
		processTranscriptome(inputFile, outputPrefix);
		log.info("Process finished");
	}

	public void processTranscriptome(String transcriptomeFile, String outPrefix) throws IOException {
		Distribution geneLengthDist = new Distribution (0,9999,200);
		Distribution transcriptLengthDist = new Distribution (0,4999,200);
		Distribution transcriptsPerGeneDist = new Distribution (0,20,1);
		Distribution numberOfExonsDist = new Distribution (0,50,1);
		Distribution cdsLengthDist = new Distribution (0,4999,200);
		Distribution proteinLengthDist = new Distribution (0,1999,100);
		
		QualifiedSequenceList sequenceNames= genome.getSequencesMetadata();
		//Load transcriptome
		GFF3TranscriptomeHandler gff3Handler = new GFF3TranscriptomeHandler(sequenceNames);
		gff3Handler.setLog(log);
		Transcriptome transcriptome = gff3Handler.loadMap(transcriptomeFile); 			
		transcriptome.fillSequenceTranscripts(genome, log);

		Map<String,Gene> visitedGenes = new HashMap<>();
		List<GenomicRegion> allRawExons = new ArrayList<GenomicRegion>();
		QualifiedSequenceList proteome = new QualifiedSequenceList();
		QualifiedSequenceList codingSequences = new QualifiedSequenceList();
		QualifiedSequenceList cdnaSequences = new QualifiedSequenceList();
		
		int n = 0;
		for(QualifiedSequence seq: sequenceNames) {
			List<Transcript> transcriptsList = transcriptome.getTranscripts(seq.getName()).asList();
			for(Transcript t: transcriptsList) {
				//Collect transcript structure statistics
				int start = t.getCodingRelativeStart();
				int end = t.getCodingRelativeEnd();
				String geneId= t.getGeneId(); 
				transcriptLengthDist.processDatapoint(t.length());
				List<GenomicRegion> exons=t.getRawExons();
				numberOfExonsDist.processDatapoint(exons.size());
				allRawExons.addAll(exons);
				
				if(!visitedGenes.containsKey(geneId)) {
					//Collect gene statistics
					List<Transcript> geneTranscripts = transcriptome.getTranscriptsByGene(geneId);
					if(geneTranscripts==null) {
						log.info("Gene not found for transcript with id"+t.getId()+" gene id: "+geneId);
						continue;
					}
					transcriptsPerGeneDist.processDatapoint(geneTranscripts.size());
					Gene g = t.getGene();
					geneLengthDist.processDatapoint(g.length());
					visitedGenes.put(geneId,g);
				}
				String comments = t.getGeneId();
				if(t.getGeneName()!=null) comments+= " "+t.getGeneName();
				CharSequence cdnaSequence = t.getCDNASequence();
				if(cdnaSequence==null) {
					continue;
				}
				QualifiedSequence qCDNA = new QualifiedSequence(t.getId(),cdnaSequence);
				qCDNA.setComments(comments);
				cdnaSequences.add(qCDNA);
				
				String protein = null;
				Codon startCodon = t.getStartCodon();
				if(startCodon==null) {
					log.info("Transcript "+t.getId()+" has an invalid start codon. Length: "+cdnaSequence.length()+" start: "+start);
					continue;
				}
				else if(!startCodon.isStart()) log.info("Transcript "+t.getId()+" does not have a standard start codon. Codon: "+startCodon.getRnaSequence()+" Length: "+cdnaSequence.length()+"start: "+start);
				Codon stopCodon = t.getStopCodon();
				if(stopCodon==null) log.info("Transcript "+t.getId()+" has an invalid stop codon. Length: "+cdnaSequence.length()+" end: "+end);
				else if(!stopCodon.isStop()) log.info("Transcript "+t.getId()+" does not have a standard stop codon. Codon: "+stopCodon.getRnaSequence()+" Length: "+cdnaSequence.length()+" end: "+end);
				
				// Process cds sequence
				CharSequence cds = t.getCDSSequence();
				if(cds == null) {
					log.info("Invalid coding sequence for transcript "+t.getId()+" Length: "+cdnaSequence.length()+"start: "+start+" end: "+end);
					continue;
				}
				QualifiedSequence qCDS = new QualifiedSequence(t.getId(),cds);
				qCDS.setComments(comments);
				codingSequences.add(qCDS);
				cdsLengthDist.processDatapoint(cds.length());
				
				// Process protein sequence
				protein = t.getProteinSequence(translator);
				if(protein == null) {
					log.info("Invalid protein for transcript "+t.getId()+" Length: "+cdnaSequence.length()+"start: "+start+" end: "+end);
					continue;
				}
				proteinLengthDist.processDatapoint(protein.length());
				QualifiedSequence qp = new QualifiedSequence(t.getId(),protein);
				qp.setComments(comments);
				proteome.add(qp);
				n++;
				if (progressNotifier!=null && n%100==0) {
					int progress = n/100;
					if (!progressNotifier.keepRunning(progress)) {
						return;
					}
				}
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
			if(codingSequences.size()>0) {
				out.println();
				out.println("CDS length");
				cdsLengthDist.printDistributionInt(out);
			}
			if(proteome.size()>0) {
				out.println();
				out.println("Protein length");
				proteinLengthDist.printDistributionInt(out);
			}
		}
		if(cdnaSequences.size()>0) {
			FastaSequencesHandler handler = new FastaSequencesHandler();
			try (PrintStream out = new PrintStream(outPrefix+"_cdna.fa")) {
				handler.saveSequences(cdnaSequences, out, 100);
			}
		}
		if(codingSequences.size()>0) {
			FastaSequencesHandler handler = new FastaSequencesHandler();
			try (PrintStream out = new PrintStream(outPrefix+"_cds.fa")) {
				handler.saveSequences(codingSequences, out, 100);
			}
		}
		if(proteome.size()>0) {
			FastaSequencesHandler handler = new FastaSequencesHandler();
			try (PrintStream out = new PrintStream(outPrefix+"_proteins.fa")) {
				handler.saveSequences(proteome, out, 100);
			}
		}
		if(visitedGenes.size()>0) {
			try (PrintStream out = new PrintStream(outPrefix+"_geneDensity.txt")) {
				densityCalculator.calculateDensity(visitedGenes.values(), out);	
			}
		}
		if(allRawExons.size()>0) {
			try (PrintStream out = new PrintStream(outPrefix+"_exonDensity.txt")) {
				densityCalculator.calculateDensity(allRawExons, out);	
			}
		}
		
	}
}
