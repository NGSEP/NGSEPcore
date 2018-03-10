package ngsep.benchmark;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;


public class VCFGoldStandardComparator {
	private ReferenceGenome genome;
	private String outVCF = null;
	private int mode = 0;
	private short minQuality = 0;
	private int posPrint = 5068;

	VCFFileWriter writer = new VCFFileWriter();
	public static void main(String[] args) throws Exception {
		VCFGoldStandardComparator instance = new VCFGoldStandardComparator();
		String referenceFile = args[0];
		String file1 = args[1];
		String file2 = args[2];
		
		if(args.length>5) {
			instance.outVCF = args[3];
			instance.mode = Integer.parseInt(args[4]);
			instance.minQuality = Short.parseShort(args[5]);
		}
		
		instance.genome = new ReferenceGenome(referenceFile);
		instance.compareFiles(file1,file2);
				
	}

	public void compareFiles(String vcfGS, String vcfTest) throws IOException {
		QualifiedSequenceList sequenceNames = genome.getSequencesMetadata();
		GenomicRegionComparator compRegion = new GenomicRegionComparator(sequenceNames);
		Map<Byte, GoldStandardComparisonCounts> countsPerType = new HashMap<>();
		countsPerType.put(GenomicVariant.TYPE_BIALLELIC_SNV, new GoldStandardComparisonCounts());
		countsPerType.put(GenomicVariant.TYPE_INDEL, new GoldStandardComparisonCounts());
		countsPerType.put(GenomicVariant.TYPE_STR, new GoldStandardComparisonCounts());
		//countsPerType.put(GenomicVariant.TYPE_UNDETERMINED, new GoldStandardComparisonCounts());
		List<CalledGenomicVariant> callsTest = VCFFileReader.loadCalledVariantsSingleIndividualVCF(vcfTest);
		int lastRowCounts = GoldStandardComparisonCounts.NUM_ROWS_COUNTS-1;
		PrintStream out = null;
		try (VCFFileReader inGS = new VCFFileReader(vcfGS)){ 
			QualifiedSequenceList seqNames = genome.getSequencesMetadata();
			inGS.setSequences(seqNames);
			inGS.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			if(outVCF!=null) {
				out = new PrintStream(outVCF);
				writer.printHeader(inGS.getHeader(), out);
			}
			Iterator<VCFRecord> itGS = inGS.iterator();
			int idxTest = 0;
			while (itGS.hasNext()) {
				VCFRecord rGS = itGS.next();
				CalledGenomicVariant callGS = rGS.getCalls().get(0);
				boolean referenceRegion = callGS.isHomozygousReference();
				byte typeGS = loadType(callGS);
				List<CalledGenomicVariant> callsIntersect = new ArrayList<>();
				while(idxTest<callsTest.size()) {
					CalledGenomicVariant callTest = callsTest.get(idxTest);
					if(callTest.isUndecided()) {
						idxTest++;
						continue;
					}
					
					byte type = loadType(callTest);
					short qualTest = callTest.getGenotypeQuality();
					int column = getGenotypeNumber(callTest);
					int cmp = compRegion.compare(callGS, callTest);
					if(posPrint==callGS.getFirst()) System.out.println("Call GS: "+callGS.getFirst()+" call test: "+callTest.getFirst()+" type: "+type+" column: "+column+" comparison: "+cmp);
					if(cmp>1) {
						column+=12;
					}
					if(cmp >1 || (cmp>=0 && referenceRegion)) {
						countsPerType.get(type).update(0,Math.min(qualTest/10, lastRowCounts),column);
						VCFRecord record = new VCFRecord(callTest, VCFRecord.DEF_FORMAT_ARRAY_QUALITY, callTest, null);
						if(mode == 2 && qualTest>=minQuality && !callTest.isHomozygousReference()) writer.printVCFRecord(record, out);
						idxTest++;
						continue;
					} else if (cmp<-1) break;
					else if (cmp == -1){
						if((callGS.isUndecided() || referenceRegion) && !callTest.isHomozygousReference()){
							int overlap = callGS.getLast()-callTest.getFirst()+1;
							int remainder = callTest.getLast() - callGS.getLast();
							if(overlap < remainder) break;
						}
					}
					if(!referenceRegion) {
						callsIntersect.add(callTest);
					} else {
						countsPerType.get(type).update(0,Math.min(qualTest/10, lastRowCounts),column);
					}
					idxTest++;
				}
				if(posPrint==callGS.getFirst()) System.out.println("GS call at "+callGS.getSequenceName()+":"+callGS.getFirst()+" type: "+typeGS+" undecided: "+callGS.isUndecided()+" Calls intersect: "+callsIntersect.size()+" next test call: "+callsTest.get(idxTest).getFirst());
				if(callGS.isUndecided() || referenceRegion ) {
					//Ignore variants that intersect with unknown gold standard sites
					continue;
				}
				//Process non-reference call in the gold standard
				int n1 = getGenotypeNumber(callGS);
				if(callsIntersect.size()==0) {
					//Gold standard variant not found
					countsPerType.get(typeGS).update(0,lastRowCounts,9+n1);
					if(mode == 1) writer.printVCFRecord(rGS, out);
					continue;
				}
				boolean covered = false;
				for(CalledGenomicVariant callTest:callsIntersect) {
					Set<String> allelesGS = buildAlleleStrings(callGS, callTest.getFirst(), callTest.getLast());
					Set<String> allelesTest = buildAlleleStrings(callTest, callGS.getFirst(), callGS.getLast());
					boolean consistent = allelesGS.equals(allelesTest); 
					byte typeTest = loadType(callTest);
					short qualTest = callTest.getGenotypeQuality();
					int n2 = getGenotypeNumber(callTest);
					VCFRecord rTest = new VCFRecord(callTest, VCFRecord.DEF_FORMAT_ARRAY_QUALITY, callTest, null);
					if(consistent) {
						GoldStandardComparisonCounts counts = countsPerType.get(typeGS);
						int k = 3*n1+n2;
						int row = Math.min(qualTest/10, lastRowCounts); 
						counts.update(0,row,k);
						counts.update(row+1,lastRowCounts,9+n1);
						if(mode == 1 &&  qualTest<minQuality) {
							writer.printVCFRecord(rGS, out);
							writer.printVCFRecord(rTest, out);
						}
						if(mode==3 && n1!=n2 &&  qualTest>=minQuality) {
							writer.printVCFRecord(rGS, out);
							writer.printVCFRecord(rTest, out);
						}
					} else {
						//False positive for inconsistent variant
						countsPerType.get(typeTest).update(0,Math.min(qualTest/10, lastRowCounts),12+n2);
						if(mode==3 && qualTest>=minQuality) {
							writer.printVCFRecord(rGS, out);
							writer.printVCFRecord(rTest, out);
						}
					}
				}
				if(!covered) {
					GoldStandardComparisonCounts counts = countsPerType.get(typeGS);
					//False negative for GS not properly called
					counts.update(0,lastRowCounts,9+n1);
					
					
				}
			}
			while(idxTest<callsTest.size()) {
				//Process test variants at the end
				CalledGenomicVariant callTest = callsTest.get(idxTest);
				if(!callTest.isUndecided()) {
					byte typeTest = loadType(callTest);
					short qualTest = callTest.getGenotypeQuality();
					int n2 = getGenotypeNumber(callTest);
					countsPerType.get(typeTest).update(0,Math.min(qualTest/10, lastRowCounts),12+n2);
					VCFRecord rTest = new VCFRecord(callTest, VCFRecord.DEF_FORMAT_ARRAY_QUALITY, callTest, null);
					if(mode == 2 && qualTest>=minQuality && !callTest.isHomozygousReference()) writer.printVCFRecord(rTest, out);
				}
			}
		} finally {
			if(out!=null) {
				out.flush();
				out.close();
			}
		}
		System.out.println("SNVs");
		countsPerType.get(GenomicVariant.TYPE_BIALLELIC_SNV).print(System.out);
		System.out.println("Indels");
		countsPerType.get(GenomicVariant.TYPE_INDEL).print(System.out);
		System.out.println("STRs/OTHER");
		countsPerType.get(GenomicVariant.TYPE_STR).print(System.out);
	}
	
	private Set<String> buildAlleleStrings(GenomicVariant v, int firstRegion, int lastRegion) {
		String seqName = v.getSequenceName();
		String left = null;
		if(firstRegion<v.getFirst()) left = genome.getReference(seqName, firstRegion, v.getFirst()-1).toString();
		String right = null;
		if(v.getLast()<lastRegion) right = genome.getReference(seqName, v.getLast()+1, lastRegion).toString();
		Set<String> allelesSet = new TreeSet<>();
		String [] alleles = v.getAlleles();
		for(int i=0;i<alleles.length;i++) {
			String allele = alleles[i];
			int cut = v.getLast()-lastRegion-(alleles[0].length()-allele.length());
			if(cut>0 && cut < allele.length()) allele = allele.substring(0, allele.length()-cut);
			String allO = "";
			if(left!=null) allO += left;
			allO+=allele;
			if(right!=null) allO += right;
			allelesSet.add(allO.toUpperCase());
		}
		return allelesSet;
	}
	
	private byte loadType(CalledGenomicVariant call) {
		byte type = call.getType();
		if(type == GenomicVariant.TYPE_BIALLELIC_SNV) return type;
		if(type == GenomicVariant.TYPE_INDEL) return type;
		if(type == GenomicVariant.TYPE_STR) return type;
		String [] alleles = call.getAlleles();
		if(alleles.length==2 && alleles[0].length()!=alleles[1].length()) return GenomicVariant.TYPE_INDEL;
		else if(alleles.length==2 && alleles[0].length()==1) return GenomicVariant.TYPE_BIALLELIC_SNV;
		return GenomicVariant.TYPE_STR;
	}

	private static int getGenotypeNumber(CalledGenomicVariant call) {
		if(call.isHeterozygous()) return 1;
		else if (!call.isHomozygousReference()) return 2;
		return 0;
	}
}
class GoldStandardComparisonCounts {
	public static final int NUM_ROWS_COUNTS = 10;
	private int [][] counts;
	public GoldStandardComparisonCounts () {
		counts = new int [NUM_ROWS_COUNTS][15];
	}
	
	public void update(int firstRow, int lastRow, int column) {
		for(int i=firstRow;i<=lastRow;i++) {
			counts[i][column]++;
		}	
	}
	
	public void print(PrintStream out) {
		for (int i=0;i<counts.length; i++) {
			out.print(""+(i*10));
			for(int j=0;j<counts[0].length;j++) {
				out.print("\t"+counts[i][j]);
			}
			out.println();
		}
	}
	
}
