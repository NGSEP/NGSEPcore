package ngsep.haplotyping;

import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import java.util.Iterator;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.math.NumberArrays;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class SingleIndividualHaplotyper {

	// Constants for default values
	public static final String DEF_ALGORITHM_NAME="Refhap";
	public static final int DEF_MIN_MQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	// Logging and progress
	private Logger log = Logger.getLogger(SingleIndividualHaplotyper.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private String alignmentsFile = null;
	private String outputFile = null;
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
	
	
	
	public static void main(String[] args) throws Exception {
		SingleIndividualHaplotyper instance = new SingleIndividualHaplotyper();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		
	}
	
	public void run () throws IOException {
		if (inputFile == null) throw new IOException("The input VCF file is a required parameter");
		if (alignmentsFile == null) throw new IOException("The file with read alignments is a required parameter");
		if(outputFile==null) process (inputFile, alignmentsFile, System.out);
		else {
			try (PrintStream out=new PrintStream(outputFile)) {
				process (inputFile, alignmentsFile, out);
			}
		}
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
			 ReadAlignmentFileReader alnReader = new ReadAlignmentFileReader(bamFilename)) {
			VCFFileHeader header = inputVCF.getHeader();
			vcfWriter = new VCFFileWriter();
			vcfWriter.printHeader(header, out);
			alnReader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_SEQUENCE);
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
						nextAln = phaseSequenceVariants(lastSeqName, hetCalls, nextAln, alnIt);
						vcfWriter.printVCFRecords(records, out);
					}
					records.clear();
					hetCalls.clear();
					lastSeqName = record.getSequenceName();
				}
				records.add(record);
				CalledGenomicVariant call = record.getCalls().get(0);
				if(call.isBiallelic() && call.getCopyNumber()==2 && !call.isUndecided()) {
					if(call.isHeterozygous()) hetCalls.add(call);
					else if (call instanceof CalledSNV) ((CalledSNV)call).setPhasingCN2(!call.isHomozygousReference());
					//else if (call instanceof CalledGenomicVariantImpl) ((CalledGenomicVariantImpl)call).setPhasedAlleles(phasedAlleles);
				}
				
			}
			if(records.size()>0) {
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
		} catch (Exception e) {
			throw new IOException("Can not load algorithm: "+algorithmName,e);
		}
		
	}
	private ReadAlignment phaseSequenceVariants(String seqName, List<CalledGenomicVariant> hetCalls, ReadAlignment nextAln, Iterator<ReadAlignment> alnIt) {
		System.err.println("Sequence: "+seqName+" Phasing "+hetCalls.size()+" het calls");
		if(nextAln!=null) System.err.println("First alignment. "+nextAln.getSequenceName()+":"+nextAln.getFirst());
		HaplotypeBlock block = new HaplotypeBlock(hetCalls);
		int i=0;
		while(nextAln!=null && nextAln.getSequenceName().equals(seqName)) {
			//Advance i
			GenomicVariant firstHetVar = null;
			while(i<hetCalls.size()) {
				firstHetVar = hetCalls.get(i);
				if(GenomicRegionPositionComparator.getInstance().compare(firstHetVar, nextAln)>=0) {
					break;
				}
				i++;
			}
			if(i==hetCalls.size()) {
				//Try to go to next alignment
				if(alnIt.hasNext()) nextAln = alnIt.next();
				else nextAln = null;
				continue;
			}
			//Extract relevant calls from alignment
			int lastAln = nextAln.getLast();
			List<Byte> calls = new ArrayList<>(50);
			int realCalls = 0;
			int first = i;
			for(int j=i;j<hetCalls.size();j++) {
				GenomicVariant var = hetCalls.get(j);
				if(var.getFirst()>lastAln) {
					break;
				}
				String [] alleles = var.getAlleles();
				CharSequence callS = nextAln.getAlleleCall(var.getFirst(), var.getLast());
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
			
			if(realCalls>1) {
				block.addFragment (first,NumberArrays.toByteArray(calls));
				if(block.getNumFragments()%1000==0) System.err.println("Added "+block.getNumFragments()+" fragments"+" calls last fragment: "+realCalls);
			}
			//Try to go to next alignment
			if(alnIt.hasNext()) nextAln = alnIt.next();
			else nextAln = null;
		}
		System.err.println("Phasing sequence "+seqName+" with "+block.getNumFragments()+" fragments");
		if(nextAln!=null) System.err.println("First alignment for next sequence. "+nextAln.getSequenceName()+":"+nextAln.getFirst());
		if(block.getNumFragments()>0) {
			phaseBlockVariants (seqName, block,hetCalls);
		}
		
		return nextAln;
	}
	private void phaseBlockVariants(String seqName, HaplotypeBlock block, List<CalledGenomicVariant> hetCalls) {
		algorithm.buildHaplotype(block);
		System.err.println("Phased sequence "+seqName+" with "+block.getNumFragments()+" fragments");
		block.phaseCallsWithHaplotype();
	}
}
