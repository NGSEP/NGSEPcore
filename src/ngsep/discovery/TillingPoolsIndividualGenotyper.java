package ngsep.discovery;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.simulation.TillingPopulationSimulator;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class TillingPoolsIndividualGenotyper {

	private Logger log = Logger.getLogger(TillingPopulationSimulator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private static HashMap<String,ArrayList<Integer>> poolConfiguration;
	private static HashMap<String,List<CalledGenomicVariant>> individualVariants;
	private static List<List<CalledGenomicVariant>> poolVariants;
	private static List<CalledGenomicVariant>[][] variantIntersection;
	private static ArrayList<Integer> numPools;
	private static int totPools;
	private ReferenceGenome genome;
	
	public static void main(String[] args) throws Exception {
		TillingPoolsIndividualGenotyper instance = new TillingPoolsIndividualGenotyper();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String poolsVCFDir = args[i++];
		String poolsDescriptor = args[i++];
		String referenceFile = args[i++];
		String outFilename = args[i++];
		instance.run(poolsVCFDir, poolsDescriptor, referenceFile, outFilename);
		
	}
	
	/**
	 * Runs the detection process itself.
	 * @param poolsVCFDir The directory with the vcf files with the variants called per pool.
	 * @param poolsDescriptor A csv file with one row per individual in the population. Each row should contain the identifier of the individual, and then its row, column,
	 * and plate pool.
	 * @param referenceFile The file with the reference sequences for each of the regions analyzed in the experiment.
	 * @param outFilename The prefix of the output vcf file.
	 * @throws IOException
	 */
	public void run(String poolsVCFDir, String poolsDescriptor, String referenceFile, String outFilename) throws IOException {
		loadPools(poolsDescriptor);
		loadPoolVariants(poolsVCFDir);
		for(int k=0;k<numPools.get(0);k++) {
			for (int j=0;j<totPools;j++)	{
				variantIntersection[k][j]= new ArrayList<CalledGenomicVariant>();
			}
		}
		
		genome = new ReferenceGenome(referenceFile);
		intersectPairs();
		callIndVariants();
		printMutations(outFilename);
	}
	
	/**
	 * Loads the csv file which indicates the pools assigned to each individual. 
	 * Also (for now) calculates the number of row, column, and plate pools (although the user could indicate this...).	
	 * @param poolsDescriptor
	 * @throws IOException
	 */

	public void loadPools(String poolsDescriptor) throws IOException {
		// TODO Auto-generated method stub
		numPools.add(0);
		numPools.add(0);
		numPools.add(0);
		
		BufferedReader reader;
		reader = new BufferedReader(new FileReader(poolsDescriptor));
		String line = reader.readLine();
		
		while (line != null) {
			line = reader.readLine();
			String[] indInfo = line.split(";");
			ArrayList<Integer> pools = new ArrayList();
			pools.add(Integer.parseInt(indInfo[1]));
			pools.add(Integer.parseInt(indInfo[2]));
			pools.add(Integer.parseInt(indInfo[3]));
			poolConfiguration.put(indInfo[0], pools);
			
			if (pools.get(0)>numPools.get(0)){
				numPools.set(0, pools.get(0));
			}
			
			if (pools.get(1)>numPools.get(1)){
				numPools.set(1, pools.get(1));
			}
			
			if (pools.get(2)>numPools.get(2)){
				numPools.set(2, pools.get(2));
			}
			
		}
		reader.close();
		
		totPools=numPools.get(0)+numPools.get(1)+numPools.get(2);
		
	}
	
	/**
	 * Loads the variants found in the vcf of each pool, all located in the directory poolsVCFDir.
	 * For now it is assumed that the vcf files are named with only the pool number (e.g.: 1.vcf), starting with the row pools, then the columns ones, and then the 
	 * plate ones.
	 * @param poolsVCFDir
	 * @throws IOException
	 */
	
	public void loadPoolVariants(String poolsVCFDir) throws IOException {
		File dir = new File(poolsVCFDir);
		int i = 0;
		File[] directoryListing = dir.listFiles();
		 if (directoryListing != null) {
		   for (File poolVCF : directoryListing) {
			    List<CalledGenomicVariant> indPool = VCFFileReader.loadCalledVariantsSingleIndividualVCF(poolVCF.getName());
			    poolVariants.add(indPool);
			    i+=1;
	    	}
	 	}
	}
	
	/**
	 * Determines which variants are common in every possible pair of pools were one (and only one) member of the pair is a column pool. 
	 * This is done to speed things up when calling the variants per individual.
	 */
	
	public void intersectPairs() {
		for(int i=0;i<numPools.get(0);i++) {
			for (int j=numPools.get(0);j<totPools;j++)	{
				List<CalledGenomicVariant> colQuery = poolVariants.get(i);
				for(int k=0;k<colQuery.size();k++) {
					CalledGenomicVariant varQuery=colQuery.get(k);
					List<CalledGenomicVariant> otherQuery = poolVariants.get(j);
					for(int l=0;l<otherQuery.size();l++) {
						if(varQuery.isCompatible(otherQuery.get(l))) {
							variantIntersection[i][j].add(varQuery);
						}
					}
					
				}
			}
		}
	}

	/**
	 * The run command itself takes the result of intersectPairs and does a the third, final, intersection based on its assigned pools.
	 * @param outFilename
	 */
	public void callIndVariants() {
		// TODO Auto-generated method stub
		for (String ind : poolConfiguration.keySet()) {
		    List<Integer> indPools = poolConfiguration.get(ind);
		    List<CalledGenomicVariant> varDef= new  ArrayList<CalledGenomicVariant>();
		    List<CalledGenomicVariant> varInt1=variantIntersection[indPools.get(0)][indPools.get(1)];
			for(int i=0;i<varInt1.size();i++) {
				CalledGenomicVariant varQuery=varInt1.get(i);
				List<CalledGenomicVariant> varInt2=variantIntersection[indPools.get(0)][indPools.get(2)];
				for(int j=0;j<varInt2.size();j++) {
					if(varQuery.isCompatible(varInt2.get(j))) {
						varDef.add(varQuery);
					}
				}
					
			}
			individualVariants.put(ind, varDef);
		}
	}
	
	/**
	 * Writes a VCF file with the simulated mutations.
	 * This is taken directly from TillingPopulationSimulator. Adjustments are done to match the variables and objects in this class.
	 * @param filename Name of the file to write
	 * @throws IOException 
	 */
	public void printMutations(String outFilename) throws IOException {
		List<CalledGenomicVariant> allCalls = new ArrayList<>();
		for(String individual:individualVariants.keySet()) {
			List<CalledGenomicVariant> callsIndv = individualVariants.get(individual);
			allCalls.addAll(callsIndv);
		}
		QualifiedSequenceList seqsMetadata = genome.getSequencesMetadata();
		GenomicRegionComparator cmp = new GenomicRegionComparator(seqsMetadata);
		Collections.sort(allCalls,cmp);
		
		VCFFileWriter vcfWriter = new VCFFileWriter();
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		for(String individual:individualVariants.keySet()) header.addDefaultSample(""+individual);
		
		try (PrintStream out = new PrintStream(outFilename)) {
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
		CalledGenomicVariant [] callsArray = new CalledGenomicVariant [individualVariants.size()];
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

}

