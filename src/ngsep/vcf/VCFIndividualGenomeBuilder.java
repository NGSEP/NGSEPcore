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
	
	//Constants
	public static final short DEF_PLOIDY = GenomicVariant.DEFAULT_PLOIDY; 
	
	// Parameters
	private ReferenceGenome genome;
	private String variantsFile = null;
	private String outputFile = null;
	private short ploidy = DEF_PLOIDY;
	
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
	
	public short getPloidy() {
		return ploidy;
	}
	public void setPloidy(short ploidy) {
		this.ploidy = ploidy;
	}
	
	public void setPloidy(String value) {
		this.setPloidy((short) OptionValuesDecoder.decode(value, Short.class));
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
		out.println("Sample ploidy:"+ ploidy);
		log.info(os.toString());
	}
	public void makeGenomeFromVCF(String vcfFile, PrintStream out) throws IOException {
		QualifiedSequenceList seqMetadata = genome.getSequencesMetadata();
		QualifiedSequenceList individualGenome = new QualifiedSequenceList();
		
		GenomicRegionSortedCollection<CalledGenomicVariant> allVariants = loadCalls(vcfFile,seqMetadata);
		if(allVariants.size()==0) {
			throw new IOException("No variants were loaded from "+vcfFile);
		}
		for(QualifiedSequence sequence:seqMetadata) {
			String seqName = sequence.getName();
			log.info("Building assembly for sequence "+seqName);
			int l = sequence.getLength();
			int nextPos = 1;
			StringBuilder [] outSequences = new StringBuilder[ploidy];
			for(int i=0;i<ploidy;i++) outSequences[i]=new StringBuilder(l);
			List<CalledGenomicVariant> sequenceCallsList = allVariants.getSequenceRegions(seqName).asList();
			for(CalledGenomicVariant call:sequenceCallsList) {
				String[] phasedAlleles = calculatePhasedAlleles(call);
				
				if(nextPos<call.getFirst()) {
					//Fill haplotypes with non variant segment
					CharSequence segment = genome.getReference(seqName, nextPos, call.getFirst()-1);
					if(segment!=null) {
						String nonVariantSegment = segment.toString().toUpperCase();
						for(int i=0;i<ploidy;i++) outSequences[i].append(nonVariantSegment);
					} else {
						log.warning("Error loading segment "+seqName+":"+nextPos+"-"+(call.getFirst()-1));
					}
				}
				for(int i=0;i<ploidy;i++) outSequences[i].append(phasedAlleles[i]);
				nextPos = call.getLast()+1;
			}
			if(nextPos<l) {
				//End of a chromosome
				CharSequence nonVarLast = genome.getReference(seqName, nextPos, l);
				if(nonVarLast!=null) {
					String nonVariantSegment = nonVarLast.toString().toUpperCase();
					for(int i=0;i<ploidy;i++) outSequences[i].append(nonVariantSegment);
				} else {
					log.warning("Error loading segment "+seqName+":"+nextPos+"-"+l);
				}
			}
			for(int i=0;i<ploidy;i++) {
				String outSeqName = seqName;
				if(ploidy>1) outSeqName+="_H"+(i+1);
				individualGenome.add(new QualifiedSequence(outSeqName, outSequences[i]));
			}
		}
		FastaSequencesHandler outHandler = new FastaSequencesHandler();
		outHandler.saveSequences(individualGenome, out, 100);
	}
	private String[] calculatePhasedAlleles(CalledGenomicVariant call) {
		String [] answer = new String [ploidy] ;
		if(ploidy==1) {
			//Polishing mode. Ensure that the reference allele is copied unless the call is homozygous alternative
			if(call.isHomozygous()) answer[0] = call.getCalledAlleles()[0];
			else answer[0] = call.getAlleles()[0];
			return answer;
		}
		if(call.isPhased()) {
			String [] phasedAlleles = call.getPhasedAlleles();
			if(phasedAlleles.length==ploidy) {
				for(int i=0;i<ploidy;i++) answer[i] = phasedAlleles[i];
				return answer;
			}
		}
		if(call.isHomozygous()) {
			String calledAllele = call.getCalledAlleles()[0];
			for(int i=0;i<ploidy;i++) answer[i] = calledAllele;
		} else {
			//TODO: Improve managenent of acn for calls with different copy number than ploidy
			String [] alleles = call.getAlleles();
			short [] acn = call.getAllelesCopyNumber();
			if (call.getCopyNumber()!=ploidy) {
				log.warning("Heterozygous call at "+call.getSequenceName()+": "+call.getFirst()+" has a copy number different than the ploidy.");	
			}
			if(acn.length!=alleles.length) {
				log.warning("Heterozygous call at "+call.getSequenceName()+": "+call.getFirst()+" has an alleles dosage array size "+acn.length+" different than the alleles array "+alleles.length);
			}
			
			int i = 0;
			for(int a=0;i<ploidy && a<alleles.length;a++) {
				for(int j=0;i<ploidy && j<acn[a];j++) {
					answer[i] = alleles[a];
					i++;
				}
			}
			for(;i<ploidy;i++) answer[i] = alleles[0];
		}
		return answer;
	}
	
	private GenomicRegionSortedCollection<CalledGenomicVariant> loadCalls(String filename, QualifiedSequenceList seqMetadata) throws IOException {
		GenomicRegionSortedCollection<CalledGenomicVariant> calls = new GenomicRegionSortedCollection<>(seqMetadata);
		try (VCFFileReader reader = new VCFFileReader(filename)){		
			Iterator<VCFRecord> it = reader.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				if(record.getVariant().isStructural()) continue;
				List<CalledGenomicVariant> genotypeCalls = record.getCalls();
				if(genotypeCalls.size()<1) continue;
				CalledGenomicVariant call = genotypeCalls.get(0);
				if(call.isHomozygousReference()) continue;
				calls.add(call);
			}
		}
		return calls;
	}
}