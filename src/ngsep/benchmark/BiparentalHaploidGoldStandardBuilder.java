package ngsep.benchmark;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import ngsep.genome.ReferenceGenome;
import ngsep.variants.AlleleCompatibilityGenomicVariantComparator;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class BiparentalHaploidGoldStandardBuilder {
	private List<String> filesParent1=new ArrayList<>();
	private List<String> filesParent2=new ArrayList<>();
	private ReferenceGenome genome;
	public static void main(String[] args) throws Exception {
		BiparentalHaploidGoldStandardBuilder instance = new BiparentalHaploidGoldStandardBuilder();
		instance.genome = new ReferenceGenome(args[0]);
		instance.loadDescriptor(args[1]);
		try (PrintStream out = new PrintStream(args[2])) {
			instance.buildGoldStandardF1(out);
		}
	}
	public void loadDescriptor(String filename) throws IOException {
		try (FileReader reader = new FileReader(filename); 
			 BufferedReader in = new BufferedReader(reader)) {
			String line = in.readLine();
			while (line != null) {
				String [] items = line.split("\t");
				if("1".equals(items[0])) filesParent1.add(items[1]);
				if("2".equals(items[0])) filesParent2.add(items[1]);
				line = in.readLine();
			}
		}
	}
	public void buildGoldStandardF1(PrintStream out) throws IOException{
		List<String> allFiles = new ArrayList<>();
		allFiles.addAll(filesParent1);
		allFiles.addAll(filesParent2);
		int n = allFiles.size();
		VCFFileWriter writer = new VCFFileWriter();
		List<VCFFileReader> readers = new ArrayList<VCFFileReader>();
		List<Iterator<VCFRecord>> iterators = new ArrayList<Iterator<VCFRecord>>();
		VCFFileHeader outHeader = VCFFileHeader.makeDefaultEmptyHeader();
		outHeader.addDefaultSample("GoldStandardF1");
		VCFRecord [] calls = new VCFRecord [n];
		for(int i=0;i<n;i++) {
			VCFFileReader reader = new VCFFileReader(allFiles.get(i));
			readers.add(reader);
			Iterator<VCFRecord> it = reader.iterator();
			if(it.hasNext()) calls[i] = it.next();
			else calls[i] = null;
			iterators.add(it);
		}
		AlleleCompatibilityGenomicVariantComparator comparator = new AlleleCompatibilityGenomicVariantComparator(genome);
		writer.printHeader(outHeader, out);
		while (true){
			VCFRecord record = createNextRecord(iterators,calls,filesParent1.size(), comparator,outHeader);
			if(record == null) break;
			writer.printVCFRecord(record, out);
		}
		for(VCFFileReader reader:readers) reader.close();
	}
	private VCFRecord createNextRecord(List<Iterator<VCFRecord>> iterators, VCFRecord[] fileRecords, int firstIdxP2, AlleleCompatibilityGenomicVariantComparator comparator, VCFFileHeader outHeader) {
		GenomicVariant nextVariant = null;
		int n = fileRecords.length;
		//Look for next variant
		for(int i=0;i<n;i++) {
			VCFRecord nextRecord = fileRecords[i];
			if(nextRecord!=null) {
				if(nextVariant==null || comparator.compare(nextRecord.getVariant(), nextVariant)<0) {
					nextVariant = nextRecord.getVariant();
				}
			}
		}
		if(nextVariant == null) return null;
		List<CalledGenomicVariant> outCalls = new ArrayList<CalledGenomicVariant>();
		boolean [] nextVarSupport = new boolean[n];
		List<GenomicVariant> consistentVars = new ArrayList<>();
		Arrays.fill(nextVarSupport, false);
		for(int i=0;i<n;i++) {
			VCFRecord record = fileRecords[i];
			if(record != null) {
				GenomicVariant var2 = record.getVariant();
				if(comparator.compare(nextVariant, var2)==0) {
					nextVarSupport[i] = true;
					consistentVars.add(var2);
					Iterator<VCFRecord> it = iterators.get(i);
					if(it.hasNext()) fileRecords[i] = it.next();
					else fileRecords[i] = null;
				}
			}
		}
		GenomicVariant outVariant = buildMergedVariant(consistentVars, comparator);
		CalledGenomicVariant outCall = decideCall(outVariant, nextVarSupport,firstIdxP2);
		outCalls.add(outCall);
		VCFRecord answer = new VCFRecord(nextVariant, new ArrayList<String>(), new ArrayList<GenomicVariantAnnotation>(), VCFRecord.DEF_FORMAT_ARRAY_MINIMAL, outCalls, outHeader);
		return answer;
	}
	/**
	 * Merging of consistent variants in the variant that will make it to the 
	 * @param consistentVars
	 * @return GenomicVariant Merged variant
	 */
	private GenomicVariant buildMergedVariant(List<GenomicVariant> consistentVars, AlleleCompatibilityGenomicVariantComparator comparator) {
		int first = -1;
		int last = -1;
		short quality = 0;
		String sequenceName = null;
		byte type= GenomicVariant.TYPE_UNDETERMINED;
		for (GenomicVariant gv: consistentVars) {
			sequenceName = gv.getSequenceName();
			if(first == -1 || first > gv.getFirst()) {
				first = gv.getFirst();
			}
			if(last == -1 || last < gv.getLast()) {
				last = gv.getLast();
			}
			if(quality<gv.getVariantQS()) quality = gv.getVariantQS();
			type = GenomicVariantImpl.mergeType(type, gv.getType());
		}
		List<String> alleles;
		if(first == last) {
			alleles = Arrays.asList(consistentVars.get(0).getAlleles());
		} else {
			alleles = buildAlleleStrings(consistentVars.get(0), first, last);
		}
		GenomicVariant gv = new GenomicVariantImpl(sequenceName, first, alleles);
		gv.setType(type);
		gv.setVariantQS(quality);
		return gv;
	}
	
	private List<String> buildAlleleStrings(GenomicVariant v, int first, int last) {
		String seqName = v.getSequenceName();
		CharSequence left = null;
		if(first <=v.getFirst()-1) left = genome.getReference(seqName, first, v.getFirst()-1);
		CharSequence right = null;
		if(last>=v.getLast()+1) right = genome.getReference(seqName, v.getLast()+1, last);
		List<String> allelesList = new ArrayList<>();
		String [] alleles = v.getAlleles();
		for(int i=0;i<alleles.length;i++) {
			String allele = "";
			if(left!=null)allele+=left;
			allele+=alleles[i];
			if(right!=null)allele+=right;
			allelesList.add(allele.toUpperCase());
		}
		return allelesList;
	}
	private CalledGenomicVariant decideCall(GenomicVariant variant, boolean[] nextVarSupport, int firstIdxP2) {
		int supportP1 = 0;
		int supportP2 = 0;
		int nP1 = firstIdxP2;
		int nP2 = nextVarSupport.length-nP1;
		
		for(int i=0;i<nextVarSupport.length;i++) {
			if(nextVarSupport[i]) {
				if(i<firstIdxP2) supportP1++;
				else supportP2++;
			}
		}
		boolean aP1 = (supportP1>=nP1-1);
		boolean aP2 = (supportP2>=nP2-1);
		byte [] indexesCalledAlleles = new byte[0];
		if ( aP1 && aP2 ) {
			indexesCalledAlleles = new byte[1];
			indexesCalledAlleles[0]=1;
		} else if ((aP1 && supportP2==0) || (aP2 && supportP1==0) ) {
			indexesCalledAlleles = new byte[2];
			indexesCalledAlleles[0]=0;
			indexesCalledAlleles[1]=1;
		}
		return new CalledGenomicVariantImpl(variant, indexesCalledAlleles);
	}
}
