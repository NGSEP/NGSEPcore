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
package ngsep.vcf;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;

/**
 * Builds a genome in fasta format from homozygous variants in a VCF file
 * @author Jorge Duitama
 *
 */
public class VCFIndividualGenomeBuilder {

	// Logging and progress
	private Logger log = Logger.getLogger(VCFIndividualGenomeBuilder.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private ReferenceGenome genome;
	private String variantsFile = null;
	private String outputFile = null;
	
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
    
	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	public void setGenome(String genomeFile) throws IOException {
		setGenome(OptionValuesDecoder.loadGenome(genomeFile,log));
	}
	
	public String getVariantsFile() {
		return variantsFile;
	}
	public void setVariantsFile(String variantsFile) {
		this.variantsFile = variantsFile;
	}
	
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	

	public static void main(String[] args) throws Exception {
		VCFIndividualGenomeBuilder instance = new VCFIndividualGenomeBuilder();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	
	public void run() throws IOException {
		logParameters();
		if(genome==null) throw new IOException("The input file with a genome assembly is a required parameter");
		if(variantsFile==null) throw new IOException("A file with variants in VCF format is required");
		if(outputFile==null) throw new IOException("An output file path is required");
		
		try (PrintStream out = new PrintStream(outputFile)) {
			makeGenomeFromVCF(variantsFile, out);
		}
		log.info("Process finished");
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		if (genome!=null) out.println("Input genome loaded from: "+genome.getFilename());
		out.println("Variants file:"+ variantsFile);
		out.println("Output file:"+ outputFile);
		log.info(os.toString());
	}
	public void makeGenomeFromVCF(String vcfFile, PrintStream out) throws IOException {
		QualifiedSequenceList seqMetadata = genome.getSequencesMetadata();
		QualifiedSequenceList individualGenome = new QualifiedSequenceList();
		
		GenomicRegionSortedCollection<GenomicVariant> allVariants = loadVariants(vcfFile,seqMetadata);
		if(allVariants.size()==0) {
			throw new IOException("No variants were loaded from "+vcfFile);
		}
		for(QualifiedSequence sequence:seqMetadata) {
			String seqName = sequence.getName();
			log.info("Building assembly for sequence "+seqName);
			int l = sequence.getLength();
			int nextPos = 1;
			StringBuilder outSequence = new StringBuilder();
			List<GenomicVariant> sequenceVarsList = allVariants.getSequenceRegions(seqName).asList();
			for(GenomicVariant var:sequenceVarsList) {
				String [] alleles = var.getAlleles();
				if(nextPos<var.getFirst()) {
					//Fill haplotypes with non variant segment
					CharSequence segment = genome.getReference(seqName, nextPos, var.getFirst()-1);
					if(segment!=null) {
						String nonVariantSegment = segment.toString().toUpperCase();
						outSequence.append(nonVariantSegment);
					} else {
						log.warning("Error loading segment "+seqName+":"+nextPos+"-"+(var.getFirst()-1));
					}
				}
				outSequence.append(alleles[1]);
				nextPos = var.getLast()+1;
			}
			if(nextPos<l) {
				//End of a chromosome
				CharSequence nonVarLast = genome.getReference(seqName, nextPos, l);
				if(nonVarLast!=null) {
					String nonVariantSegment = nonVarLast.toString().toUpperCase();
					outSequence.append(nonVariantSegment);
				} else {
					log.warning("Error loading segment "+seqName+":"+nextPos+"-"+l);
				}
			}
			individualGenome.add(new QualifiedSequence(seqName, outSequence));
		}
		FastaSequencesHandler outHandler = new FastaSequencesHandler();
		outHandler.saveSequences(individualGenome, out, 100);
	}
	
	private GenomicRegionSortedCollection<GenomicVariant> loadVariants(String filename, QualifiedSequenceList seqMetadata) throws IOException {
		GenomicRegionSortedCollection<GenomicVariant> variants = new GenomicRegionSortedCollection<>(seqMetadata);
		try (VCFFileReader reader = new VCFFileReader(filename)){		
			Iterator<VCFRecord> it = reader.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				GenomicVariant variant = record.getVariant();
				
				List<CalledGenomicVariant> genotypeCalls = record.getCalls();
				if(genotypeCalls.size()<1) continue;
				CalledGenomicVariant call = genotypeCalls.get(0);
				if(!call.isBiallelic()) continue;
				if(call.isHeterozygous()) continue;
				if(call.isHomozygousReference()) continue;
				if(call.getGenotypeQuality()<20) continue;
				variants.add(variant);
			}
		}
		return variants;
	}
}