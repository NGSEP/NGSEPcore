package ngsep.haplotyping;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.main.CommandsDescriptor;
import ngsep.math.NumberArrays;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class SingleIndividualHaplotyper {

	private String algorithmClassName = "ngsep.haplotyping.RefhapSIHAlgorithm";
	private SIHAlgorithm algorithm;
	
	private int minMQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	public static void main(String[] args) throws Exception {
		SingleIndividualHaplotyper instance = new SingleIndividualHaplotyper();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		
		String vcfFilename = args[i++];
		String bamFilename = args[i++];
		instance.process (vcfFilename,bamFilename,System.out);
	}
	
	/**
	 * @return the minMQ
	 */
	public int getMinMQ() {
		return minMQ;
	}

	/**
	 * @param minMQ the minMQ to set
	 */
	public void setMinMQ(int minMQ) {
		this.minMQ = minMQ;
	}
	
	/**
	 * @param minMQ the minMQ to set
	 */
	public void setMinMQ(Integer minMQ) {
		this.setMinMQ(minMQ.intValue());
	}
	
	public void setAlgorithmName(String name) {
		algorithmClassName = "ngsep.haplotyping."+name+"SIHAlgorithm";
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
		VCFFileReader inputVCF = null;
		ReadAlignmentFileReader alnReader = null;
		VCFFileWriter vcfWriter = new VCFFileWriter();
		try
		{
			inputVCF = new VCFFileReader(vcfFilename);
			VCFFileHeader header = inputVCF.getHeader();
			vcfWriter = new VCFFileWriter();
			vcfWriter.printHeader(header, out);
			alnReader = new ReadAlignmentFileReader(bamFilename);
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
		finally{
			if(inputVCF!=null) inputVCF.close();
			if(alnReader!=null) alnReader.close();
		}
	}
	private void loadAlgorithm() throws IOException {
		try {
			algorithm = (SIHAlgorithm)Class.forName(algorithmClassName).newInstance();
		} catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
			throw new IOException("Can not load algorithm with classname: "+algorithmClassName,e);
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
