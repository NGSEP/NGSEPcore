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
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class TillingPoolsIndividualGenotyper {

	private Logger log = Logger.getLogger(TillingPoolsIndividualGenotyper.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private static HashMap<String,ArrayList<Integer>> poolConfiguration;
	private static HashMap<String,List<CalledGenomicVariant>> individualVariants;
	private static ArrayList<List<CalledGenomicVariant>> poolVariants;
	private static ArrayList<ArrayList<ArrayList<CalledGenomicVariant>>> variantIntersection;
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
		int num_columnPools = Integer.valueOf(args[i++]);
		int num_rowPools = Integer.valueOf(args[i++]);
		int num_platePools = Integer.valueOf(args[i++]);
		instance.run(poolsVCFDir, poolsDescriptor, referenceFile, outFilename, num_columnPools, num_rowPools, num_platePools);
		
	}
	
	/**
	 * Runs the detection process itself and yields the variants per individual.
	 * @param poolsVCFDir The directory with the vcf files with the variants called per pool.
	 * @param poolsDescriptor A csv file with one row per individual in the population. Each row should contain the identifier of the individual, and then its row, column,
	 * and plate pool.
	 * @param referenceFile The file with the reference sequences for each of the regions analyzed in the experiment.
	 * @param outFilename The prefix of the output vcf file.
	 * @param num_columnPools Number of column pools.
	 * @param num_rowPools Number of row pools.
	 * @param num_platePools Number of plate pools.
	 * @throws IOException
	 */
	public void run(String poolsVCFDir, String poolsDescriptor, String referenceFile, String outFilename, int num_columnPools, int num_rowPools, int num_platePools) throws IOException {
		loadPools(poolsDescriptor, num_columnPools, num_rowPools, num_platePools);
		loadPoolVariants(poolsVCFDir);
		
		variantIntersection = new ArrayList<ArrayList<ArrayList<CalledGenomicVariant>>>();
		
		for(int k=0;k<numPools.get(0);k++) {
			ArrayList<ArrayList<CalledGenomicVariant>> temp = new ArrayList<ArrayList<CalledGenomicVariant>>();
			for (int j=0;j<totPools;j++)	{
				temp.add(new ArrayList<CalledGenomicVariant>());
			}
			variantIntersection.add(temp);
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

	public void loadPools(String poolsDescriptor, int num_columnPools, int num_rowPools, int num_platePools) throws IOException {
		// TODO Auto-generated method stub
		numPools = new ArrayList<Integer>();
		numPools.add(num_columnPools);
		numPools.add(num_rowPools);
		numPools.add(num_platePools);
		
		poolConfiguration= new HashMap<String,ArrayList<Integer>>(); 
		
		BufferedReader reader;
		reader = new BufferedReader(new FileReader(poolsDescriptor));
		String line = reader.readLine();
		
		while (line != null) {

			String[] indInfo = line.split(";");
			ArrayList<Integer> pools = new ArrayList<>();
			pools.add(Integer.parseInt(indInfo[1]));
			pools.add(Integer.parseInt(indInfo[2]));
			pools.add(Integer.parseInt(indInfo[3]));
			poolConfiguration.put(indInfo[0], pools);
			
			line = reader.readLine();
			
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
		poolVariants=new ArrayList<List<CalledGenomicVariant>>();
		
		File dir = new File(poolsVCFDir);
		int i = 0;
		File[] directoryListing = dir.listFiles();
		Arrays.sort(directoryListing);	
		 if (directoryListing != null) {
		   for (File poolVCF : directoryListing) {
			   System.out.println(poolVCF.getName());
			    List<CalledGenomicVariant> indPool = VCFFileReader.loadCalledVariantsSingleIndividualVCF(poolsVCFDir+poolVCF.getName());
			    for(int j=0;j<indPool.size();j++) {
			    	indPool.get(j).setSampleId(String.valueOf(i)+String.valueOf(j));
			    }
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
		for(int i=0;i<numPools.get(0)-1;i++) {
			for (int j=numPools.get(0)-1;j<totPools-1;j++)	{
				List<CalledGenomicVariant> colQuery = poolVariants.get(i);
				for(int k=0;k<colQuery.size();k++) {
					CalledGenomicVariant varQuery=colQuery.get(k);
					List<CalledGenomicVariant> otherQuery = poolVariants.get(j);
					for(int l=0;l<otherQuery.size();l++) {
						if(varQuery.isCompatible(otherQuery.get(l))) {
							variantIntersection.get(i).get(j).add(varQuery);
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
		
		individualVariants = new HashMap<String,List<CalledGenomicVariant>>(); 
		// TODO Auto-generated method stub
		for (String ind : poolConfiguration.keySet()) {
			int z=0;
		    List<Integer> indPools = poolConfiguration.get(ind);
		    List<CalledGenomicVariant> varDef= new  ArrayList<CalledGenomicVariant>();
		    List<CalledGenomicVariant> varInt1=variantIntersection.get(indPools.get(0)).get(indPools.get(1));
			for(int i=0;i<varInt1.size();i++) {
				CalledGenomicVariant varQuery=varInt1.get(i);
				List<CalledGenomicVariant> varInt2=variantIntersection.get(indPools.get(0)).get(indPools.get(2));
				for(int j=0;j<varInt2.size();j++) {
					if(varQuery.isCompatible(varInt2.get(j))) {
						varQuery.setSampleId(ind.substring(3)+"_"+String.valueOf(z));
						varDef.add(varQuery);
						z+=1;
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
		
		try (PrintStream out = new PrintStream(outFilename+".vcf")) {
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
		/**System.out.println(individualVariants.size());**/
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
		/** Made change here to correct error, but we need to verify whether it affects the process of calling variants (pending further revision with simulated data)**/
			String[] ident = call.getSampleId().split("_");
			callsArray[Integer.parseInt(ident[0])] = new CalledGenomicVariantImpl(finalVariant, call.getCalledAlleles());
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

