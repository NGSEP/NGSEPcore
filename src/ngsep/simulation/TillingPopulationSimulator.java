package ngsep.simulation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;

public class TillingPopulationSimulator {
	private Logger log = Logger.getLogger(SingleIndividualSimulator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	public static final int DEF_MUTATIONS=1000;
	public static final int DEF_INDIVIDUALS=500;
	
	private ReferenceGenome genome;
	private int numIndividuals = DEF_INDIVIDUALS;
	private List<GenomicRegion> sequencedRegions;
	
	//Variants indexed by individual
	private List<SimulatedDiploidIndividual> individuals;
	private List<List<SimulatedDiploidIndividual>> pools;
	
	
	public static void main(String[] args) throws Exception {
		TillingPopulationSimulator instance = new TillingPopulationSimulator();
		String referenceFile = args[0];
		String sequencedRegionsFile = args[1];
		String outPrefix = args[2];
		instance.genome = new ReferenceGenome(referenceFile);
		instance.runSimulation (sequencedRegionsFile, outPrefix);
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
		//Step 2: create random mutations within the given regions and assign each mutation to a random individual
		//Use the reference genome to derive reference alleles. Crete objects of the class SNV as random mutations
		
		
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


	public void printMutations(String filename) {
		// TODO Auto-generated method stub
		
	}
	/**
	 * Simulate the pools to sequence the simulated individuals
	 */
	public void simulatePools() {
		
		
	}

	/**
	 * Simulates sequencing reads for the given pool of individuals
	 * @param pool Individuals with allele sequences to simulate reads
	 * @param file1 Output file for first end of paired end reads
	 * @param file2 Output file for second end of paired end reads
	 */
	public void simulatePoolReads(List<SimulatedDiploidIndividual> pool, String file1, String file2) {
		
		
	}	
}
class SimulatedDiploidIndividual {
	private int id;
	private List<GenomicVariant> mutationsAllele1 = new ArrayList<>();
	private List<GenomicVariant> mutationsAllele2 = new ArrayList<>();
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
		if(random.nextBoolean()) {
			mutationsAllele1.add(mutation);
		} else {
			mutationsAllele2.add(mutation);
		}
		
	}
	public List<GenomicVariant> getMutations() {
		List<GenomicVariant> answer = new ArrayList<>();
		answer.addAll(mutationsAllele1);
		answer.addAll(mutationsAllele2);
		return answer;
	}
	/**
	 * Calculates the sequences within an individual from the given mutations
	 * @param regions
	 * @param referenceSequences
	 */
	public void buildAlleleSequences (List<GenomicRegion> regions, List<DNAMaskedSequence> referenceSequences) {
		
	}
	public List<DNAMaskedSequence> getAlleleSequences () {
		return Collections.unmodifiableList(alleleSequences);
	}
	public DNAMaskedSequence getRandomSequence() {
		int idx = random.nextInt(alleleSequences.size());
		return alleleSequences.get(idx);
	}
}
