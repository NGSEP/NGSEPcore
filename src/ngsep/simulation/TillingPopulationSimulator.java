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

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.SNV;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class TillingPopulationSimulator {
	private Logger log = Logger.getLogger(TillingPopulationSimulator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	public static final int DEF_MUTATIONS=1000;
	public static final int DEF_INDIVIDUALS=500;
	public static final int DEF_NUM_FRAGMENTS_POOL=1000000;
	public static final int DEF_READ_LENGTH=250;
	public static final double DEF_ERROR_RATE=0.01;
	
	private ReferenceGenome genome;
	private int numIndividuals = DEF_INDIVIDUALS;
	private int numFragments = DEF_NUM_FRAGMENTS_POOL;
	private int readLength = DEF_READ_LENGTH;
	private double errorRate = DEF_ERROR_RATE;
	private List<GenomicRegion> sequencedRegions;
	
	//Variants indexed by individual
	private List<SimulatedDiploidIndividual> individuals;
	private List<List<SimulatedDiploidIndividual>> pools;
	
	
	public static void main(String[] args) throws Exception {
		TillingPopulationSimulator instance = new TillingPopulationSimulator();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String referenceFile = args[i++];
		String sequencedRegionsFile = args[i++];
		String outPrefix = args[i++];
		instance.genome = new ReferenceGenome(referenceFile);
		instance.runSimulation (sequencedRegionsFile, outPrefix);
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
	 * @return the progressNotifier
	 */
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}

	/**
	 * @param progressNotifier the progressNotifier to set
	 */
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
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
	 * @return the numIndividuals
	 */
	public int getNumIndividuals() {
		return numIndividuals;
	}
	
	/**
	 * @param numIndividuals the numIndividuals to set
	 */
	public void setNumIndividuals(int numIndividuals) {
		this.numIndividuals = numIndividuals;
	}

	public void setNumIndividuals(String value) {
		this.setNumIndividuals((int) OptionValuesDecoder.decode(value, Integer.class));
	}

	/**
	 * @return the numFragments
	 */
	public int getNumFragments() {
		return numFragments;
	}

	/**
	 * @param numFragments the numFragments to set
	 */
	public void setNumFragments(int numFragments) {
		this.numFragments = numFragments;
	}
	
	public void setNumFragments(String value) {
		this.setNumFragments((int) OptionValuesDecoder.decode(value, Integer.class));
	}

	/**
	 * @return the readLength
	 */
	public int getReadLength() {
		return readLength;
	}

	/**
	 * @param readLength the readLength to set
	 */
	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}

	public void setReadLength(String value) {
		this.setReadLength((int) OptionValuesDecoder.decode(value, Integer.class));
	}

	/**
	 * @return the errorRate
	 */
	public double getErrorRate() {
		return errorRate;
	}

	/**
	 * @param errorRate the errorRate to set
	 */
	public void setErrorRate(double errorRate) {
		this.errorRate = errorRate;
	}
	
	public void setErrorRate(String value) {
		this.setErrorRate((double)OptionValuesDecoder.decode(value, Double.class));
	}

	public void runSimulation(String sequencedRegionsFile, String outPrefix) throws IOException {
		loadSequencedRegions(sequencedRegionsFile);
		simulatePopulation();
		printMutations(outPrefix+".vcf");
		simulatePools();
		for(int i=0;i<pools.size();i++) {
			List<SimulatedDiploidIndividual> pool = pools.get(i);
			simulatePoolReads(pool, outPrefix+"P"+i+"_1.fastq", outPrefix+"P"+i+"_2.fastq");
		}
		
		
	}


	public void loadSequencedRegions(String sequencedRegionsFile) throws IOException {
		SimpleGenomicRegionFileHandler regionsFH = new SimpleGenomicRegionFileHandler();
		sequencedRegions = regionsFH.loadRegions(sequencedRegionsFile);
		QualifiedSequenceList seqsMetadata = genome.getSequencesMetadata();
		Collections.sort(sequencedRegions,new GenomicRegionComparator(seqsMetadata));
	}

	/**
	 * Simulates a tilling population
	 * Creates individuals, simulates random mutations and builds the mutated allele sequences per individual
	 */
	public void simulatePopulation() {
		Random random = new Random();
		String alphabet = DNASequence.BASES_STRING;
		//Step 1: create objects for individuals
		individuals = new ArrayList<>(numIndividuals);
		for(int i=0;i<numIndividuals;i++) {
			individuals.add(new SimulatedDiploidIndividual(i));
		}
		//TODO: Step 2: create random mutations within the given regions and assign each mutation to a random individual
		//Use the reference genome to derive reference alleles. Create objects of the class SNV as random mutations
		for(int j=0; j < DEF_MUTATIONS; j++) {
			SimulatedDiploidIndividual targetInd = individuals.get(random.nextInt(individuals.size()));
			GenomicRegion targetGR = sequencedRegions.get(random.nextInt(sequencedRegions.size()));
			int location = random.nextInt(targetGR.getLast())+targetGR.getFirst();
			char refBase = genome.getReferenceBase(targetGR.getSequenceName(), location);
			String mutated = alphabet.replaceAll(Character.toString(refBase) , "");
			GenomicVariant variant = new SNV(targetGR.getSequenceName(), location, refBase, mutated.charAt(random.nextInt(3)));
			targetInd.addMutation(variant);
		}
			
		//Step 3: build individual allele sequences from the mutations
		List<DNAMaskedSequence> referenceSequences = getReferenceSequencesRegions();
		for(SimulatedDiploidIndividual ind: individuals) {
			ind.buildAlleleSequences(sequencedRegions, referenceSequences);
		}
	}

	/**
	 * Extracts from the genome the reference sequences of the regions to sequence in the experiment
	 * @return List<DNAMaskedSequence> Reference allele of each sequenced region
	 */
	private List<DNAMaskedSequence> getReferenceSequencesRegions() {
		// TODO Auto-generated method stub
		return null;
	}

	/**
	 * Writes a VCF file with the simulated mutations
	 * @param filename Name of the file to write
	 * @throws IOException 
	 */
	public void printMutations(String filename) throws IOException {
		List<CalledGenomicVariant> allCalls = new ArrayList<>();
		for(SimulatedDiploidIndividual individual:individuals) {
			List<CalledGenomicVariant> callsIndv = individual.getCalls();
			allCalls.addAll(callsIndv);
		}
		QualifiedSequenceList seqsMetadata = genome.getSequencesMetadata();
		GenomicRegionComparator cmp = new GenomicRegionComparator(seqsMetadata);
		Collections.sort(allCalls,cmp);
		
		VCFFileWriter vcfWriter = new VCFFileWriter();
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		for(SimulatedDiploidIndividual indv:individuals) header.addDefaultSample(""+indv.getId());
		
		try (PrintStream out = new PrintStream(filename)) {
			vcfWriter.printHeader(header, out);
			List<CalledGenomicVariant> altCallsVar = new ArrayList<>();
			for(CalledGenomicVariant call: allCalls) {
				if(altCallsVar.size()>0 && cmp.compare(altCallsVar.get(0),call)!=0) {
					VCFRecord record = buildRecord (altCallsVar,header);
					vcfWriter.printVCFRecord(record, out);
					altCallsVar.clear();
				}
				altCallsVar.add(call);
			}
			if(altCallsVar.size()>0) {
				VCFRecord record = buildRecord (altCallsVar,header);
				vcfWriter.printVCFRecord(record, out);
			}
		}
	}
	private VCFRecord buildRecord(List<CalledGenomicVariant> altCallsVar, VCFFileHeader header) {
		Set<String> allelesSet = new TreeSet<>();
		CalledGenomicVariant [] callsArray = new CalledGenomicVariant [numIndividuals];
		String sequenceName=null;
		int first=0;
		String reference=null;
		for(CalledGenomicVariant call:altCallsVar) {
			int n = allelesSet.size();
			if(n==0) {
				sequenceName = call.getSequenceName();
				first = call.getFirst();
				reference = call.getReference();
			}
			allelesSet.addAll(Arrays.asList(call.getAlleles()));
		}
		List<String> alleles = new ArrayList<>(allelesSet.size());
		alleles.add(reference);
		for(String allele:allelesSet) {
			if(!allele.equals(reference)) alleles.add(allele);
		}
		GenomicVariant finalVariant = new GenomicVariantImpl(sequenceName, first, alleles);
		for(CalledGenomicVariant call:altCallsVar) {
			int i = Integer.parseInt(call.getSampleId()); 
			callsArray[i] = new CalledGenomicVariantImpl(finalVariant, call.getCalledAlleles());
		}
		List<CalledGenomicVariant> finalCalls = new ArrayList<>();
		for(int i=0;i<callsArray.length;i++) {
			if(callsArray[i]==null) {
				CalledGenomicVariant refCall = new CalledGenomicVariantImpl(finalVariant, 0);
				finalCalls.add(refCall);
			} else {
				finalCalls.add(callsArray[i]);
			}
		}
		return new VCFRecord(finalVariant, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL, finalCalls, header);
	}




	/**
	 * Simulate the pools to sequence the simulated individuals
	 */
	public void simulatePools() {
		//TODO: Implement
		
	}

	/**
	 * Simulates sequencing reads for the given pool of individuals
	 * @param pool Individuals with allele sequences to simulate reads
	 * @param file1 Output file for first end of paired end reads
	 * @param file2 Output file for second end of paired end reads
	 */
	public void simulatePoolReads(List<SimulatedDiploidIndividual> pool, String file1, String file2) {
		Random random = new Random();
		
		//For each fragment select a random individual, then select an allele sequence at random and build the reads from the two ends of the sequence 
		
	}	
}
class SimulatedDiploidIndividual {
	private int id;
	private Map<String,GenomicVariant> mutationsAllele1 = new HashMap<>();
	private Map<String,GenomicVariant> mutationsAllele2 = new HashMap<>();
	private List<DNAMaskedSequence> alleleSequences = new ArrayList<>();
	private static Random random = new Random();
	public SimulatedDiploidIndividual(int id) {
		super();
		this.id = id;
	}
	public int getId() {
		return id;
	}
	public void addMutation (GenomicVariant mutation) {
		String key = buildKey(mutation);
		if(random.nextBoolean()) {
			if(!mutationsAllele1.containsKey(key)) mutationsAllele1.put(key,mutation);
		} else {
			if(!mutationsAllele1.containsKey(key)) mutationsAllele2.put(key,mutation);
		}
	}
	public List<CalledGenomicVariant> getCalls() {
		List<CalledGenomicVariant> answer = new ArrayList<>();
		String sampleId = ""+id;
		for(String key:mutationsAllele1.keySet()) {
			GenomicVariant mutation1 = mutationsAllele1.get(key);
			GenomicVariant mutation2 = mutationsAllele2.get(key);
			// TODO: Conciliate better overlapping mutations
			CalledGenomicVariant call;
			if(mutation2 == null) {
				call = new CalledGenomicVariantImpl(mutation1, 1);
			} else if (mutation1.isCompatible(mutation2)) {
				call = new CalledGenomicVariantImpl(mutation1, 2);
			} else {
				List<String> alleles = new ArrayList<>();
				alleles.add(mutation1.getReference());
				alleles.add(mutation1.getAlleles()[1]);
				alleles.add(mutation2.getAlleles()[1]);
				GenomicVariantImpl triallelicVariant = new GenomicVariantImpl(mutation1.getSequenceName(), mutation1.getFirst(), alleles);
				byte [] calledAlleles = {1,2};
				call = new CalledGenomicVariantImpl (triallelicVariant,calledAlleles);
			}
			call.setSampleId(sampleId);
			answer.add(call);
		}
		for(String key:mutationsAllele2.keySet()) {
			GenomicVariant mutation1 = mutationsAllele1.get(key);
			GenomicVariant mutation2 = mutationsAllele2.get(key);
			CalledGenomicVariant call;
			if(mutation1 == null) {
				call = new CalledGenomicVariantImpl(mutation2, 1);
				call.setSampleId(sampleId);
				answer.add(call);
			}
		}
		return answer;
	}
	private String buildKey(GenomicVariant mutation) {
		return mutation.getSequenceName()+"-"+mutation.getFirst()+"-"+mutation.getLast();
	}
	
	/**
	 * Calculates the sequences within an individual from the given mutations
	 * @param regions
	 * @param referenceSequences
	 */
	public void buildAlleleSequences (List<GenomicRegion> regions, List<DNAMaskedSequence> referenceSequences) {
		//TODO: implement
	}
	public List<DNAMaskedSequence> getAlleleSequences () {
		return Collections.unmodifiableList(alleleSequences);
	}
	public DNAMaskedSequence getRandomSequence() {
		int idx = random.nextInt(alleleSequences.size());
		return alleleSequences.get(idx);
	}
}
