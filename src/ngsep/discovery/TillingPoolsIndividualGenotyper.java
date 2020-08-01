package ngsep.discovery;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
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
	
	private HashMap<String,ArrayList<Integer>> poolConfiguration;
	private Map<Integer,List<CalledGenomicVariant>> poolVariants;
	private HashMap<String,List<CalledGenomicVariant>> individualVariants;
	
	private String [] inputFiles;
	private String poolsDescriptor = null;
	private String outputFile = null;
	private ReferenceGenome genome;
	
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

	public String getPoolsDescriptor() {
		return poolsDescriptor;
	}
	public void setPoolsDescriptor(String poolsDescriptor) {
		this.poolsDescriptor = poolsDescriptor;
	}

	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
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

	public static void main(String[] args) throws Exception {
		TillingPoolsIndividualGenotyper instance = new TillingPoolsIndividualGenotyper();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.inputFiles = Arrays.copyOfRange(args, i, args.length);
		instance.run();
		
	}
	
	/**
	 * Runs the detection process itself and yields the variants per individual.
	 * @throws IOException If files can not be read or written
	 */
	public void run() throws IOException {
		loadPools();
		loadPoolVariants();
		callIndVariants();
		saveVCF();
	}
	
	/**
	 * Loads the file which indicates the pools assigned to each individual.
	 * @throws IOException
	 */

	public void loadPools() throws IOException {
		poolConfiguration= new HashMap<String,ArrayList<Integer>>(); 
		try (BufferedReader reader = new BufferedReader(new FileReader(poolsDescriptor))) {
			String line = reader.readLine();
			if(line.startsWith("Individual")) line = reader.readLine();
			while (line != null) {
				String[] indInfo = line.split(";");
				ArrayList<Integer> poolNumbers = new ArrayList<>();
				for(int i=1;i<indInfo.length;i++) {
					poolNumbers.add(Integer.parseInt(indInfo[i]));
				}
				poolConfiguration.put(indInfo[0], poolNumbers);
				line = reader.readLine();
			}
		}	
	}
	
	/**
	 * Loads the variants found in the vcf of each pool
	 * @throws IOException
	 */
	
	public void loadPoolVariants() throws IOException {
		poolVariants=new HashMap<Integer,List<CalledGenomicVariant>>();
		for (String poolVCF : inputFiles) {
			log.info("Loading variants from "+poolVCF);
			List<CalledGenomicVariant> callsPool = VCFFileReader.loadCalledVariantsSingleIndividualVCF(poolVCF);
			if (callsPool.size()>0) {
				//TODO: use variant files descriptor
				int poolId = Integer.parseInt(callsPool.get(0).getSampleId());
				poolVariants.put(poolId,callsPool);
			}
			//TODO: add sample ids for VCFs created with other tools
			/*for(int j=0;j<callsPool.size();j++) {
				callsPool.get(j).setSampleId(String.valueOf(i));
			}*/
			
		}
	}

	/**
	 * The run command itself takes the result of intersectPairs and does a the third, final, intersection based on its assigned pools.
	 * @param outFilename
	 */
	public void callIndVariants() {
		individualVariants = new HashMap<String,List<CalledGenomicVariant>>();
		QualifiedSequenceList seqsMetadata = genome.getSequencesMetadata();
		GenomicRegionComparator cmp = new GenomicRegionComparator(seqsMetadata);
		for (String ind : poolConfiguration.keySet()) {
			List<Integer> indPools = poolConfiguration.get(ind);
			int firstPoolId = indPools.get(0);
			List<CalledGenomicVariant> callsFirstPool = poolVariants.get(firstPoolId);
			if(callsFirstPool==null) continue;
			List<CalledGenomicVariant> callsPerIndividual = new  ArrayList<CalledGenomicVariant>();
			for(CalledGenomicVariant call:callsFirstPool) {
				boolean consistent = true;
				for(int i=1;i<indPools.size() && consistent;i++) {
					List<CalledGenomicVariant> callsNextPool = poolVariants.get(indPools.get(i));
					consistent = callsNextPool!=null && Collections.binarySearch(callsNextPool, call, cmp)>=0;
				}
				if(consistent) callsPerIndividual.add(call);
		    }
			System.out.println("Calls for individual "+ind+" : "+callsPerIndividual.size());
			individualVariants.put(ind, callsPerIndividual);
		}
	}
	
	/**
	 * Writes a VCF file with the identified mutations.
	 * @throws IOException 
	 */
	public void saveVCF() throws IOException {
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
		
		try (PrintStream out = new PrintStream(outputFile)) {
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

