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
package ngsep.simulation;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class SingleIndividualSimulator {

	// Constants for default values
	public static final double DEF_SNV_RATE=0.001;
	public static final double DEF_INDEL_RATE=0.0001;
	public static final double DEF_MUTATED_STR_FRACTION=0.1;
	public static final int DEF_STR_UNIT_INDEX=14;
	public static final short DEF_PLOIDY=2;
	public static final String DEF_SAMPLE_ID="Simulated";
	
	// Logging and progress
	private Logger log = Logger.getLogger(SingleIndividualSimulator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private ReferenceGenome genome;
	private String outputPrefix = null;
	private double snvRate = DEF_SNV_RATE;
	private double indelRate = DEF_INDEL_RATE;
	private double mutatedSTRFraction = DEF_MUTATED_STR_FRACTION;
	private String strsFile = null;
	private int strUnitIndex = DEF_STR_UNIT_INDEX;
	private String sampleId = DEF_SAMPLE_ID;
	private byte ploidy = DEF_PLOIDY;
	
	// Model attributes
	private GenomicRegionSortedCollection<SNV> snvs;
	private GenomicRegionSortedCollection<GenomicVariant> indels;
	private GenomicRegionSortedCollection<STR> strs;
	
	//Calls including reference allele regions and simulated variants
	private List<CalledGenomicVariant> genomicCalls;
	private List<QualifiedSequence> individualGenome;
	
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
	
	public String getOutputPrefix() {
		return outputPrefix;
	}
	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}
	
	public double getSnvRate() {
		return snvRate;
	}
	public void setSnvRate(double snvRate) {
		this.snvRate = snvRate;
	}
	public void setSnvRate(String value) {
		this.setSnvRate((double)OptionValuesDecoder.decode(value, Double.class));
	}

	public double getIndelRate() {
		return indelRate;
	}
	public void setIndelRate(double indelRate) {
		this.indelRate = indelRate;
	}
	public void setIndelRate(String value) {
		this.setIndelRate((double)OptionValuesDecoder.decode(value, Double.class));
	}
	
	public double getMutatedSTRFraction() {
		return mutatedSTRFraction;
	}
	public void setMutatedSTRFraction(double mutatedSTRFraction) {
		this.mutatedSTRFraction = mutatedSTRFraction;
	}
	public void setMutatedSTRFraction(String value) {
		this.setMutatedSTRFraction((double)OptionValuesDecoder.decode(value, Double.class));
	}

	public byte getPloidy() {
		return ploidy;
	}
	public void setPloidy(byte ploidy) {
		this.ploidy = ploidy;
	}
	public void setPloidy(String value) {
		this.setPloidy((byte)OptionValuesDecoder.decode(value, Byte.class));
	}

	public int getStrUnitIndex() {
		return strUnitIndex;
	}
	public void setStrUnitIndex(int strUnitIndex) {
		this.strUnitIndex = strUnitIndex;
	}
	public void setStrUnitIndex(String value) {
		this.setStrUnitIndex((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public String getSampleId() {
		return sampleId;
	}
	public void setSampleId(String sampleId) {
		this.sampleId = sampleId;
	}
	
	public String getStrsFile() {
		return strsFile;
	}
	public void setStrsFile(String strsFile) {
		this.strsFile = strsFile;
	}
	
	public static void main(String[] args) throws Exception {
		SingleIndividualSimulator instance = new SingleIndividualSimulator();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	
	public void run() throws IOException {
		logParameters();
		if(genome==null) throw new IOException("A file with the reference genome is a required parameter");
		if(outputPrefix==null) throw new IOException("A prefix for the output files is a required parameter");
		loadSTRs();
		simulateVariants();
		buildAssembly();
		try (PrintStream outGenome = new PrintStream(outputPrefix+".fa")) {
			saveIndividualGenome(outGenome);
		}
		try (PrintStream outVariants = new PrintStream(outputPrefix+".vcf")) {
			saveVariants(outVariants);
		}
	}

	public void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		if (genome!=null) out.println("Genome for simulation loaded from file: "+genome.getFilename());
		out.println("Prefix of the output files:"+ outputPrefix);
		out.println("SNV rate: "+snvRate);
		out.println("Indel rate: "+indelRate);
		if(strsFile!=null) {
			out.println("Mutated STR fraction: "+mutatedSTRFraction);
			out.println("Load STRs from: "+strsFile);
			out.println("Column within the STRs file with the unit sequence (zero-based): "+strUnitIndex);
		}
		out.println("Sample id: "+sampleId);
		out.println("Ploidy: "+ploidy);
		log.info(os.toString());
	}
	public void loadSTRs() throws IOException {
		if(strsFile==null) return; 
		log.info("Loading STRs from "+strsFile);
		QualifiedSequenceList seqMetadata = genome.getSequencesMetadata();
		GenomicRegionSortedCollection<STR> strsC = new GenomicRegionSortedCollection<>(seqMetadata);
		try (FileReader reader = new FileReader(strsFile);
			BufferedReader in = new BufferedReader(reader)) {
			String line = in.readLine();
			while (line != null) {
				String [] items = line.split(" |\t");
				QualifiedSequence seq = seqMetadata.get(items[0]);
				if(seq==null) {
					log.warning("Can not load line "+line+". Unrecognized sequence name");
					line = in.readLine();
					continue;
				}
				String seqName = seq.getName();
				int first = Math.max(1, Integer.parseInt(items[1])-1);
				int last = Math.min(Integer.parseInt(items[2])+1,seq.getLength()-1);
				String unit = items[strUnitIndex];
				CharSequence segment = genome.getReference(seqName, first, last);
				if(segment == null) {
					log.warning("Can not load line "+line+". Sequence segment not found");
					line = in.readLine();
					continue;
				}
				String reference = segment.toString().toUpperCase();
				List<String> alleles = new ArrayList<>();
				alleles.add(reference);
				STR str = new STR(seqName, first, last, alleles, unit);
				str.setVariantQS((short) 255);
				strsC.add(str);
				line = in.readLine();
			}
		}
		strs = new GenomicRegionSortedCollection<>(seqMetadata);
		for(QualifiedSequence seq:seqMetadata) {
			List<STR> seqSTRs = strsC.getSequenceRegions(seq.getName()).asList();
			STR lastStr = null;
			for(STR str:seqSTRs) {
				if(lastStr==null || (lastStr.getLast()<str.getFirst())) {
					strs.add(str);
					lastStr = str;
				}
			}
		}
		
		log.info("Loaded "+strs.size()+" non overlapping STRs");
	}
	
	public void simulateVariants() {
		simulateIndels();
		simulateSNVs();
		simulateVariantSTRs();
	}
	
	private void simulateIndels() {
		QualifiedSequenceList seqMetadata = genome.getSequencesMetadata();
		long length = genome.getTotalLength();
		indels = new GenomicRegionSortedCollection<>(seqMetadata);
		int numIndels = (int) Math.round(length*indelRate);
		Random random = new Random();
		Set<Long> selected = new HashSet<>();
		for(int i=0;i<numIndels;i++) {
			while (true) {
				long absoluteFirst = random.nextLong()%length;
				if(selected.contains(absoluteFirst)) continue;
				selected.add(absoluteFirst);
				long currentFirst = 0;
				String seqName = null; 
				int first = 0;
				for(QualifiedSequence seq: seqMetadata) {
					if(currentFirst+seq.getLength()>absoluteFirst) {
						seqName = seq.getName();
						first = (int)(absoluteFirst - currentFirst) + 1;
						break;
					}
					currentFirst+=seq.getLength();
				}
				int eventLength = random.nextInt(10)+1;
				boolean deletion = random.nextBoolean();
				int last = first+1;
				if(deletion) last = first + eventLength + 1;
				
				if (strs!=null && strs.findSpanningRegions(seqName, first, last).size()>0) continue;
				List<String> alleles = new ArrayList<>();
				CharSequence segment = genome.getReference(seqName, first, last);
				if(segment == null) continue;
				String refAllele = segment.toString().toUpperCase();
				if(refAllele == null) continue;
				alleles.add(refAllele);
				String altAllele;
				if(deletion) {
					altAllele = ""+refAllele.charAt(0)+""+refAllele.charAt(refAllele.length()-1);
				} else {
					altAllele = ""+refAllele.charAt(0)+makeRandomDNA(eventLength)+refAllele.charAt(refAllele.length()-1);
				}
				alleles.add(altAllele);
				GenomicVariantImpl indel = new GenomicVariantImpl(seqName, first, last, alleles);
				indel.setType(GenomicVariant.TYPE_INDEL);
				indel.setVariantQS((short) 255);
				indels.add(indel);
				break;
			}
		}
		log.info("Simulated "+indels.size()+" indels");
	}

	private String makeRandomDNA(int eventLength) {
		Random r = new Random();
		StringBuilder randomSequence = new StringBuilder();
		for(int j=0;j<eventLength;j++) {
			int bpI = r.nextInt(4);
			randomSequence.append(DNASequence.BASES_STRING.charAt(bpI));
		}
		return randomSequence.toString();
	}

	private void simulateSNVs() {
		QualifiedSequenceList seqMetadata = genome.getSequencesMetadata();
		long length = genome.getTotalLength();
		int numSNVs = (int) Math.round(length*snvRate);
		snvs = new GenomicRegionSortedCollection<>(seqMetadata);
		Random random = new Random();
		Set<Long> selected = new HashSet<>();
		for(int i=0;i<numSNVs;i++) {
			while (true) {
				long absoluteFirst = random.nextLong()%length;
				if(selected.contains(absoluteFirst)) continue;
				selected.add(absoluteFirst);
				long currentFirst = 0;
				String seqName = null; 
				int pos = 0;
				for(QualifiedSequence seq: seqMetadata) {
					if(currentFirst+seq.getLength()>absoluteFirst) {
						seqName = seq.getName();
						pos = (int)(absoluteFirst - currentFirst) + 1;
						break;
					}
					currentFirst+=seq.getLength();
				}
				if (strs!=null && strs.findSpanningRegions(seqName, pos).size()>0) continue;
				if (indels!=null && indels.findSpanningRegions(seqName, pos).size()>0) continue;
				char refBase = genome.getReferenceBase(seqName, pos);
				int refIdx = DNASequence.BASES_STRING.indexOf(refBase);
				if(refIdx==-1) continue;
				int altIdx = refIdx;
				while(refIdx==altIdx) altIdx = random.nextInt(DNASequence.BASES_STRING.length());
				char altBase = DNASequence.BASES_STRING.charAt(altIdx);
				SNV snv = new SNV(seqName, pos, refBase, altBase);
				snv.setVariantQS((short) 255);
				snvs.add(snv);
				break;
			}
		}
		log.info("Simulated "+snvs.size()+" SNVs");
	}

	private void simulateVariantSTRs() {
		if(strs==null) return;
		List<STR> strsList = strs.asList();
		Random random = new Random();
		int numSTRs = (int) Math.round(strsList.size()*mutatedSTRFraction);
		log.info("Simulating mutations in "+numSTRs+" STRs from a total of "+strsList.size()+" loaded STRs");
		Set<Integer> selected = new HashSet<>();
		for(int i=0;i<numSTRs;i++) {
			while (true) {
				int j = random.nextInt(strsList.size());
				if(selected.contains(j)) continue;
				selected.add(j);
				STR str = strsList.get(j);
				String reference = str.getReference();
				String unit = str.getUnitSequence();
				int copies = (reference.length()-2)/unit.length();
				int affectedCopies = random.nextInt(copies)+1; 
				int eventLength = affectedCopies*unit.length();
				boolean deletion = random.nextBoolean();
				String altAllele;
				if(deletion) {
					int end1 = reference.length() - eventLength - 1; 
					altAllele = reference.substring(0, end1);
				} else {
					altAllele = reference.substring(0, reference.length()-1);
					for(int k=0;k<affectedCopies;k++) altAllele+=unit;
				}
				altAllele+=reference.charAt(reference.length()-1);
				str.addAllele(altAllele);
				//System.out.println("Alleles str: "+str.getAlleles().length);
				break;
			}
		}
		log.info("Simulated mutations in "+numSTRs+" STRs");
	}
	
	public void buildAssembly() {
		
		QualifiedSequenceList seqMetadata = genome.getSequencesMetadata();
		genomicCalls = new ArrayList<>();
		individualGenome = new QualifiedSequenceList();
		GenomicRegionSortedCollection<GenomicVariant> allVariants = new GenomicRegionSortedCollection<>(seqMetadata);
		if(snvs!=null) allVariants.addAll(snvs);
		if (indels!=null) allVariants.addAll(indels);
		if (strs!=null) {
			for(STR str: strs) {
				if(str.isBiallelic()) allVariants.add(str);
			}
		}
		if(allVariants.size()==0) {
			throw new RuntimeException("No variants were simulated");
		}
		for(QualifiedSequence sequence:seqMetadata) {
			String seqName = sequence.getName();
			log.info("Building assembly for sequence "+seqName);
			int l = sequence.getLength();
			int nextPos = 1;
			StringBuilder [] haplotypes = new StringBuilder [ploidy];
			for(int i=0;i<haplotypes.length;i++) haplotypes[i] = new StringBuilder();
			
			int numHet = 0;
			Random random = new Random();
			List<GenomicVariant> sequenceVarsList = allVariants.getSequenceRegions(seqName).asList();
			for(GenomicVariant var:sequenceVarsList) {
				String [] alleles = var.getAlleles();
				if(nextPos<var.getFirst()) {
					//Fill haplotypes with non variant segment
					CharSequence segment = genome.getReference(seqName, nextPos, var.getFirst()-1);
					if(segment==null) {
						log.warning("Error loading segment "+seqName+":"+nextPos+"-"+(var.getFirst()-1));
					}
					String nonVariantSegment = segment.toString().toUpperCase();
					//Create reference call
					CalledGenomicVariant refCall = createReferenceCall(seqName, nextPos, nonVariantSegment);
					genomicCalls.add(refCall);
					//Update haplotype sequences
					for(int i=0;i<haplotypes.length;i++) (haplotypes[i]).append(nonVariantSegment);
				}
				//Simulate genotype as alternative allele count (always homozygous alternative for haploids
				byte altAlleleCount = 1;
				if(ploidy>1) altAlleleCount = (byte) (random.nextInt(ploidy)+1);
				boolean homozygousAlt = (altAlleleCount == ploidy);
				//Defaults for homozygous alternative
				byte [] indexesCalledAlleles= {(byte)1};
				short [] allelesCopyNumber= {(short)0,(short)ploidy};
				byte [] indexesPhasedAlleles= new byte [ploidy];
				Arrays.fill(indexesPhasedAlleles, (byte)1);
				if(!homozygousAlt) {
					//Values for heterozygous
					numHet++;
					indexesCalledAlleles= new byte [2];
					indexesCalledAlleles[0]=0;
					indexesCalledAlleles[1]=1;
					//Determine alleles copy number from simulated alternative count
					allelesCopyNumber[0] = (short) (ploidy-altAlleleCount);
					allelesCopyNumber[1] = altAlleleCount;
					
					//Simulate random assignment of alleles in haplotypes
					if(allelesCopyNumber[0] <= allelesCopyNumber[1] ) {
						randomDistribute(indexesPhasedAlleles,(byte)0,allelesCopyNumber[0],random);
					} else {
						Arrays.fill(indexesPhasedAlleles, (byte)0);
						randomDistribute(indexesPhasedAlleles,(byte)1,allelesCopyNumber[1],random);
					}
				}
				
				//Create variant call
				CalledGenomicVariantImpl call = new CalledGenomicVariantImpl(var, indexesCalledAlleles);
				call.setAllelesCopyNumber(allelesCopyNumber);
				call.setIndexesPhasedAlleles(indexesPhasedAlleles);
				genomicCalls.add(call);
				
				//Update haplotype sequences
				for(int i=0;i<indexesPhasedAlleles.length;i++) {
					byte nextAlleleIdx = indexesPhasedAlleles[i];
					(haplotypes[i]).append(alleles[nextAlleleIdx]);
				}
				nextPos = var.getLast()+1;
			}
			if(nextPos<l) {
				//End of a chromosome
				CharSequence nonVarLast = genome.getReference(seqName, nextPos, l-1);
				if(nonVarLast!=null) {
					String nonVariantSegment = nonVarLast.toString().toUpperCase();
					//Create reference call
					CalledGenomicVariant refCall = createReferenceCall(seqName, nextPos, nonVariantSegment);
					genomicCalls.add(refCall);
					//Update haplotype sequences
					for(int i=0;i<haplotypes.length;i++) (haplotypes[i]).append(nonVariantSegment);
				}
			}
			for(int i=0;i<haplotypes.length;i++) {
				String haplotype = haplotypes[i].toString();
				individualGenome.add(new QualifiedSequence(sampleId+"_"+seqName+"_Hap_"+i, haplotype));
			}
			log.info("Simulated "+numHet+" heterozygous calls for sequence "+seqName);
		}
	}

	public CalledGenomicVariant createReferenceCall(String seqName, int nextPos, String nonVariantSegment) {
		List<String> refCallAlleles = new ArrayList<>();
		refCallAlleles.add(nonVariantSegment);
		byte [] idxCalledAlleles = {0};
		CalledGenomicVariant refCall = new CalledGenomicVariantImpl(seqName, nextPos, refCallAlleles, idxCalledAlleles);
		return refCall;
	}

	private void randomDistribute(byte[] array, byte value, int number, Random random) {
		for(int i=0;i<number;i++) {
			while(true) {
				int j = random.nextInt(array.length);
				if(array[j]!=value) {
					array[j] = value;
					break;
				}
			}
		}
		
	}

	public void saveIndividualGenome(PrintStream out) {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		handler.saveSequences(individualGenome, out, 100);		
	}
	
	public void saveVariants(PrintStream out) {
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		Sample s = new Sample(sampleId);
		s.setNormalPloidy(ploidy);
		header.addSample(s, true);
		VCFFileWriter writer = new VCFFileWriter();
		writer.printHeader(header, out);
		for(CalledGenomicVariant call:genomicCalls) {
			VCFRecord record = new VCFRecord(call, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL, call, header);
			writer.printVCFRecord(record, out);
		}
	}

}
class STR extends GenomicVariantImpl implements GenomicVariant {
	private String unitSequence;
	public STR(String sequenceName, int first, int last, List<String> alleles, String unitSequence) {
		super(sequenceName, first, last, alleles);
		setType(GenomicVariant.TYPE_STR);
		this.unitSequence = unitSequence;
	}
	/**
	 * @return the unitSequence
	 */
	public String getUnitSequence() {
		return unitSequence;
	}
	/**
	 * @param unitSequence the unitSequence to set
	 */
	public void setUnitSequence(String unitSequence) {
		this.unitSequence = unitSequence;
	}
	

}
