package ngsep.haplotyping;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import java.util.Iterator;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.math.NumberArrays;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class SingleIndividualHaplotyper {

	// Constants for default values
	public static final String ALGORITHM_NAME_DGS="DGS";
	public static final String ALGORITHM_NAME_REFHAP="Refhap";
	public static final String DEF_ALGORITHM_NAME=ALGORITHM_NAME_REFHAP;
	
	public static final int DEF_MIN_MQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	// Logging and progress
	private Logger log = Logger.getLogger(SingleIndividualHaplotyper.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private String alignmentsFile = null;
	private String outputFile = null;
	private ReferenceGenome genome = null;
	private String algorithmName = DEF_ALGORITHM_NAME;
	private SIHAlgorithm algorithm;
	private int minMQ = DEF_MIN_MQ;
	
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

	public String getAlignmentsFile() {
		return alignmentsFile;
	}
	public void setAlignmentsFile(String alignmentsFile) {
		this.alignmentsFile = alignmentsFile;
	}

	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}

	
	public String getAlgorithmName() {
		return algorithmName;
	}
	public void setAlgorithmName(String name) {
		algorithmName = name;
	}
	
	public int getMinMQ() {
		return minMQ;
	}
	public void setMinMQ(int minMQ) {
		this.minMQ = minMQ;
	}
	public void setMinMQ(String value) {
		this.setMinMQ((int)OptionValuesDecoder.decode(value, Integer.class));
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
		SingleIndividualHaplotyper instance = new SingleIndividualHaplotyper();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	
	public void run () throws IOException {
		logParameters();
		if (inputFile == null) throw new IOException("The input VCF file is a required parameter");
		if (alignmentsFile == null) throw new IOException("The file with read alignments is a required parameter");
		if(outputFile==null) process (inputFile, alignmentsFile, System.out);
		else {
			try (PrintStream out=new PrintStream(outputFile)) {
				process (inputFile, alignmentsFile, out);
			}
		}
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input VCF file: "+inputFile);
		out.println("Alignments file: "+alignmentsFile);
		if(outputFile != null) out.println("Output VCF file: "+outputFile);
		else out.println("VCF file written to standard output");
		out.println("Minimum mapping quality: "+minMQ);
		out.println("Algorithm: "+algorithmName);
		log.info(""+os.toString());
	}
	
	/**
	 * 
	 * @param vcfFilename Input VCF
	 * @param bamFilename Input BAM file
	 * @param out Output VCF
	 */
	public void process(String vcfFilename, String bamFilename, PrintStream out) throws IOException {
		loadAlgorithm();
		List<VCFRecord> records = new ArrayList<>();
		List<CalledGenomicVariant> hetCalls = new ArrayList<>();
		
		
		VCFFileWriter vcfWriter = new VCFFileWriter();
		try (VCFFileReader inputVCF = new VCFFileReader(vcfFilename);
			 ReadAlignmentFileReader alnReader = new ReadAlignmentFileReader(bamFilename,genome)) {
			VCFFileHeader header = inputVCF.getHeader();
			vcfWriter = new VCFFileWriter();
			vcfWriter.printHeader(header, out);
			alnReader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_ALIGNMENT_SEQUENCE);
			alnReader.setMinMQ(minMQ);
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			filterFlags+=ReadAlignment.FLAG_MULTIPLE_ALN;
			alnReader.setFilterFlags(filterFlags);
			Iterator<ReadAlignment> alnIt = alnReader.iterator();
			ReadAlignment nextAln = alnIt.next();
			String lastSeqName = null;
			Iterator<VCFRecord> iter = inputVCF.iterator();
			while(iter.hasNext())
			{
				VCFRecord record = iter.next();
				if(!record.getSequenceName().equals(lastSeqName)) {
					if(records.size()>0) {
						log.info("Phasing "+records.size()+" variants from VCF file for sequence "+lastSeqName+ " heterozygous calls: "+hetCalls.size());
						nextAln = phaseSequenceVariants(lastSeqName, hetCalls, nextAln, alnIt);
						vcfWriter.printVCFRecords(records, out);
					}
					records.clear();
					hetCalls.clear();
					lastSeqName = record.getSequenceName();
				}
				records.add(record);
				CalledGenomicVariant call = record.getCalls().get(0);
				//TODO: Multiallelic and polyploid
				if(call.isBiallelic() && call.getCopyNumber()==2 && !call.isUndecided()) {
					if(call.isHeterozygous()) hetCalls.add(call);
					else if (call instanceof CalledSNV) ((CalledSNV)call).setPhasingCN2(!call.isHomozygousReference());
					else if (call instanceof CalledGenomicVariantImpl) ((CalledGenomicVariantImpl)call).setPhasedHomozygous();
				}
				
			}
			if(records.size()>0) {
				log.info("Phasing "+records.size()+" variants from VCF file for sequence "+lastSeqName);
				phaseSequenceVariants(lastSeqName, hetCalls, nextAln, alnIt);
				vcfWriter.printVCFRecords(records, out);
			}
		}
	}
	private void loadAlgorithm() throws IOException {
		try {
			String algorithmClassName = "ngsep.haplotyping."+algorithmName+"SIHAlgorithm";
			Class<?> algClass = Class.forName(algorithmClassName);
			Constructor<?> constructor = algClass.getDeclaredConstructors()[0];
			algorithm = (SIHAlgorithm) constructor.newInstance();
			algorithm.setLog(log);
		} catch (Exception e) {
			throw new IOException("Can not load algorithm: "+algorithmName,e);
		}
		
	}
	private ReadAlignment phaseSequenceVariants(String seqName, List<CalledGenomicVariant> hetCalls, ReadAlignment nextAln, Iterator<ReadAlignment> alnIt) throws IOException {
		log.info("Sequence: "+seqName+" Phasing "+hetCalls.size()+" het calls");
		if(nextAln!=null) log.info("First alignment. "+nextAln.getSequenceName()+":"+nextAln.getFirst());
		List<ReadAlignment> sequenceAlignments = new ArrayList<ReadAlignment>();
		while(nextAln!=null && nextAln.getSequenceName().equals(seqName)) {
			sequenceAlignments.add(nextAln);
			if(alnIt.hasNext()) nextAln = alnIt.next();
			else nextAln = null;
		}
		phaseSequenceVariants(seqName, hetCalls, sequenceAlignments);
		if(nextAln!=null) System.err.println("First alignment for next sequence. "+nextAln.getSequenceName()+":"+nextAln.getFirst());
		return nextAln;
	}
	
	public List<HaplotypeBlock> phaseSequenceVariants (String seqName, List<CalledGenomicVariant> hetCalls, List<ReadAlignment> alignments) throws IOException {
		//HaplotypeBlock block = new HaplotypeBlock(hetCalls);
		List<HaplotypeFragment> fragments = new ArrayList<>();
		List<HaplotypeBlock> answer = new ArrayList<HaplotypeBlock>();
		int i=0;
		int firstNextBlock = 0;
		int lastNextBlock = -1;
		for(ReadAlignment aln:alignments) {
			//Advance i
			GenomicVariant firstHetVar = null;
			while(i<hetCalls.size()) {
				firstHetVar = hetCalls.get(i);
				if(firstHetVar.getFirst()==aln.getFirst() || GenomicRegionPositionComparator.getInstance().compare(firstHetVar, aln)>=0) {
					break;
				}
				i++;
			}
			if(i==hetCalls.size()) break;
			//Extract relevant calls from alignment
			int lastAln = aln.getLast();
			List<Byte> calls = new ArrayList<>(50);
			int realCalls = 0;
			int first = i;
			int firstAln = i;
			for(int j=i;j<hetCalls.size();j++) {
				GenomicVariant var = hetCalls.get(j);
				if(var.getFirst()>lastAln) {
					break;
				}
				String [] alleles = var.getAlleles();
				CharSequence callS = aln.getAlleleCall(var.getFirst(), var.getLast());
				String call = null;
				if(callS!=null) {
					call = callS.toString();
				}
				if(alleles[0].equals(call)) {
					calls.add(CalledGenomicVariant.ALLELE_REFERENCE);
					realCalls++;
				} else if(alleles[1].equals(call)) {
					calls.add(CalledGenomicVariant.ALLELE_ALTERNATIVE);
					realCalls++;
				} else if (calls.size()==0) {
					first=j+1;
				} else {
					calls.add(CalledGenomicVariant.ALLELE_UNDECIDED);
				}
			}
			//Trim last undecided calls
			for(int j=calls.size()-1;j>=0;j--) {
				Byte call = calls.get(j);
				if(call!=CalledGenomicVariant.ALLELE_UNDECIDED) {
					break;
				}
				calls.remove(j);
			}
			
			if(realCalls==0) continue;
			if(lastNextBlock>=0 && firstAln>lastNextBlock) {
				CalledGenomicVariant lastCall = hetCalls.get(lastNextBlock);
				CalledGenomicVariant nextCall = hetCalls.get(firstAln);
				log.info("Discontiguity in haplotype block for sequence: "+seqName+". Last SNP with information "+lastNextBlock +" "+lastCall.getFirst()+" next SNP: "+firstAln+" "+nextCall.getFirst()+" next alignment: "+aln);
				if(fragments.size()>0) {	
					if(algorithm==null) loadAlgorithm();
					List<CalledGenomicVariant> blockCalls = selectBlockCalls(hetCalls,firstNextBlock,lastNextBlock);
					HaplotypeBlock block = new HaplotypeBlock(blockCalls,fragments);
					algorithm.buildHaplotype(block);
					block.phaseCallsWithHaplotype(firstNextBlock, lastNextBlock);
					log.info("Phased block of "+seqName+" between "+firstNextBlock+" and "+lastNextBlock+" with "+blockCalls.size()+" variants and "+block.getNumFragments()+" fragments. MEC: "+block.calculateMECCurrentHaplotypes()+" calls proportion: "+block.calculateRelativeCallsProportion());
					answer.add(block);
				}
				fragments = new ArrayList<>();
				firstNextBlock = firstAln;
				lastNextBlock = first;
			}
			if(first<firstNextBlock) System.err.println("Error calculating start of fragment. First: "+first+" firstNextBlock: "+firstNextBlock+" aln : "+aln);
			fragments.add(new HaplotypeFragment(aln.getReadNumber(), first-firstNextBlock,NumberArrays.toByteArray(calls)));
			if(fragments.size()%1000==0) log.info("Added "+fragments.size()+" fragments"+" calls last fragment: "+realCalls);
			lastNextBlock=Math.max(lastNextBlock, first+calls.size()-1);
			
		}
		
		if(fragments.size()>0) {	
			if(algorithm==null) loadAlgorithm();
			List<CalledGenomicVariant> blockCalls = selectBlockCalls(hetCalls,firstNextBlock,lastNextBlock);
			HaplotypeBlock block = new HaplotypeBlock(blockCalls,fragments);
			algorithm.buildHaplotype(block);
			block.phaseCallsWithHaplotype(firstNextBlock, lastNextBlock);
			log.info("Phased block of "+seqName+" between "+firstNextBlock+" and "+lastNextBlock+" with "+blockCalls.size()+" variants and "+block.getNumFragments()+" fragments. MEC: "+block.calculateMECCurrentHaplotypes()+" calls proportion: "+block.calculateRelativeCallsProportion());
			answer.add(block);
		}
		return answer;
	}
	private List<CalledGenomicVariant> selectBlockCalls(List<CalledGenomicVariant> hetCalls, int firstNextBlock, int lastNextBlock) {
		List<CalledGenomicVariant> answer = new ArrayList<>(lastNextBlock-firstNextBlock+1);
		for(int i=firstNextBlock;i<=lastNextBlock;i++) answer.add(hetCalls.get(i));
		return answer;
	}
}
