package ngsep.discovery;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
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
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class TillingPoolsIndividualGenotyper {

	private Logger log = Logger.getLogger(TillingPoolsIndividualGenotyper.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private List<String> individualIds;
	private Map<String,List<String>> poolConfiguration;
	private Map<String,List<CalledGenomicVariant>> poolVariants;
	private HashMap<String,List<CalledGenomicVariant>> individualVariants;
	
	private String [] inputFiles;
	private String poolsDescriptor = null;
	private String outputFile = null;
	private ReferenceGenome genome;
	private int maxPools = 0;
	private boolean onlyBiallelic = false; 
	
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
	
	public String[] getInputFiles() {
		return inputFiles;
	}
	public void setInputFiles(String[] inputFiles) {
		this.inputFiles = inputFiles;
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
	
	public int getMaxPools() {
		return maxPools;
	}
	public void setMaxPools(int maxPools) {
		this.maxPools = maxPools;
	}
	
	public void setMaxPools(String value) {
		setMaxPools((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public boolean isOnlyBiallelic() {
		return onlyBiallelic;
	}
	public void setOnlyBiallelic(boolean onlyBiallelic) {
		this.onlyBiallelic = onlyBiallelic;
	}
	public void setOnlyBiallelic(Boolean onlyBiallelic) {
		this.setOnlyBiallelic(onlyBiallelic.booleanValue());
	}
	
	public Map<String, List<String>> getPoolConfiguration() {
		return poolConfiguration;
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
		individualIds = new ArrayList<String>();
		poolConfiguration = new HashMap<String, List<String>>();	 
		try (BufferedReader reader = new BufferedReader(new FileReader(poolsDescriptor))) {
			String line = reader.readLine();
			if(line.startsWith("Individual")) line = reader.readLine();
			while (line != null) {
				String[] indInfo = line.split(";");
				individualIds.add(indInfo[0]);
				ArrayList<String> poolIds = new ArrayList<>();
				for(int i=1;i<indInfo.length;i++) {
					poolIds.add(indInfo[i]);
				}
				poolConfiguration.put(indInfo[0], poolIds);
				line = reader.readLine();
			}
		}
	}
	
	/**
	 * Loads the variants found in the vcf of each pool
	 * @throws IOException
	 */
	
	public void loadPoolVariants() throws IOException {
		poolVariants=new HashMap<String,List<CalledGenomicVariant>>();
		
		for (String poolVCF : inputFiles) {
			log.info("Loading variants from "+poolVCF);
			try (VCFFileReader reader = new VCFFileReader(poolVCF)) {
				Iterator<VCFRecord> it = reader.iterator();
				while(it.hasNext()) {
					VCFRecord record = it.next();
					List<CalledGenomicVariant> callsRecord = record.getCalls();
					for(CalledGenomicVariant call:callsRecord) {
						if(call.isUndecided() ) continue;
						if(call.isHomozygousReference()) continue;
						String poolId = call.getSampleId();
						List<CalledGenomicVariant> callsPool = poolVariants.computeIfAbsent(poolId, v->new ArrayList<CalledGenomicVariant>());
						callsPool.add(call);
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
		
		QualifiedSequenceList seqsMetadata = genome.getSequencesMetadata();
		GenomicRegionComparator cmp = new GenomicRegionComparator(seqsMetadata);
		for (String ind : poolConfiguration.keySet()) {
			List<String> indPools = poolConfiguration.get(ind);
			if(indPools==null || indPools.size()==0) System.err.println("Individual "+ind+" inconsistent in configuration. found 0 pools");
			String firstPoolId = indPools.get(0);
			List<CalledGenomicVariant> callsFirstPool = poolVariants.get(firstPoolId);
			if(callsFirstPool==null) continue;
			List<CalledGenomicVariant> callsPerIndividual = new  ArrayList<CalledGenomicVariant>();
			for(CalledGenomicVariant call:callsFirstPool) {
				if(!call.isBiallelic()) continue;
				
				boolean consistent = true;
				for(int i=1;i<indPools.size() && consistent;i++) {
					List<CalledGenomicVariant> callsNextPool = poolVariants.get(indPools.get(i));
					consistent = callsNextPool!=null && Collections.binarySearch(callsNextPool, call, cmp)>=0;
				}
				if(consistent) {
					CalledGenomicVariantImpl callIndv = new CalledGenomicVariantImpl(call, CalledGenomicVariant.GENOTYPE_HETERO);
					callIndv.setSampleId(ind);
					callsPerIndividual.add(callIndv);
				}
		    }
			System.out.println("Calls for individual "+ind+" : "+callsPerIndividual.size());
			individualVariants.put(ind, callsPerIndividual);
		}
	}
	
	private Map<String, Integer> calculateVariantPoolCounts(GenomicRegionComparator cmp) {
		Map<String, Integer> answer = new HashMap<String, Integer>();
		for(List<CalledGenomicVariant> callsPool:poolVariants.values()) {
			for(CalledGenomicVariant call: callsPool) {
				String key = getKey(call);
				answer.compute(key, (k,v)->(v==null?1:v+1));
			}
		}
		return answer;
	}
	private String getKey(CalledGenomicVariant call) {
		return call.getSequenceName()+":"+call.getFirst()+"-"+call.getLast();
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
		Map<String,Integer> variantPoolCounts = calculateVariantPoolCounts(cmp);
		
		VCFFileWriter vcfWriter = new VCFFileWriter();
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		for(String individual:individualIds) header.addDefaultSample(individual);
		
		try (PrintStream out = new PrintStream(outputFile)) {
			vcfWriter.printHeader(header, out);
			Map<String,CalledGenomicVariant> altCallsVar = new HashMap<String, CalledGenomicVariant>();
			CalledGenomicVariant firstCallGroup = null;
			for(CalledGenomicVariant call: allCalls) {
				if(firstCallGroup!=null && cmp.compare(firstCallGroup,call)!=0) {
					int poolsCount = variantPoolCounts.get(getKey(firstCallGroup));
					VCFRecord record = buildRecord (altCallsVar,poolsCount, header);
					if(passFilters(record, poolsCount)) vcfWriter.printVCFRecord(record, out);
					altCallsVar.clear();
					firstCallGroup = null; 
				}
				if(firstCallGroup==null) firstCallGroup = call;
				altCallsVar.put(call.getSampleId(),call);
			}
			if(firstCallGroup!=null) {
				int poolsCount = variantPoolCounts.get(getKey(firstCallGroup));
				VCFRecord record = buildRecord (altCallsVar, poolsCount, header);
				if(passFilters(record, poolsCount)) vcfWriter.printVCFRecord(record, out);
			}
		}
	}
	private boolean passFilters(VCFRecord record, int poolsCount) {
		if(maxPools>0 && poolsCount >maxPools) return false;
		if(onlyBiallelic && !record.getVariant().isBiallelic()) return false;
		return true;
	}
	private VCFRecord buildRecord(Map<String,CalledGenomicVariant> altCallsVar, int poolsCount, VCFFileHeader header) {
		Set<String> allelesSet = new TreeSet<>();
		String sequenceName=null;
		int first=0;
		String reference=null;
		for(CalledGenomicVariant call:altCallsVar.values()) {
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
		List<CalledGenomicVariant> finalCalls = new ArrayList<>();
		for(int i=0;i<individualIds.size();i++) {
			String sampleId = individualIds.get(i);
			CalledGenomicVariant altCall = altCallsVar.get(sampleId);
			if(altCall==null) {
				CalledGenomicVariant refCall = new CalledGenomicVariantImpl(finalVariant, CalledGenomicVariant.GENOTYPE_HOMOREF);
				refCall.setSampleId(individualIds.get(i));
				finalCalls.add(refCall);
			} else {
				CalledGenomicVariant altFinalCall = new CalledGenomicVariantImpl(finalVariant, altCall.getCalledAlleles());
				altFinalCall.setSampleId(sampleId);
				finalCalls.add(altFinalCall);
			}
		}
		VCFRecord record = new VCFRecord(finalVariant, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL, finalCalls, header);
		record.addAnnotation(new GenomicVariantAnnotation(finalVariant, "COUNTPOOLS", poolsCount));
		return record;
	}

}

