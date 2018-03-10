package ngsep.benchmark;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.AlleleCompatibilityGenomicVariantComparator;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.SNV;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;


public class VCFGoldStandardComparator {
	private ReferenceGenome genome;
	private String outVCF = null;
	private int mode = 0;
	private short minQuality = 0;

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

	private void compareFiles(String vcfGS, String vcfTest) throws IOException {
		AlleleCompatibilityGenomicVariantComparator comparator = new AlleleCompatibilityGenomicVariantComparator(genome);
		Map<Byte, GoldStandardComparisonCounts> countsPerType = new HashMap<>();
		countsPerType.put(GenomicVariant.TYPE_BIALLELIC_SNV, new GoldStandardComparisonCounts());
		countsPerType.put(GenomicVariant.TYPE_INDEL, new GoldStandardComparisonCounts());
		countsPerType.put(GenomicVariant.TYPE_STR, new GoldStandardComparisonCounts());
		//countsPerType.put(GenomicVariant.TYPE_UNDETERMINED, new GoldStandardComparisonCounts());
		int lastRowCounts = GoldStandardComparisonCounts.NUM_ROWS_COUNTS-1;
		PrintStream out = null;
		try (VCFFileReader inGS = new VCFFileReader(vcfGS);
		     VCFFileReader inTest = new VCFFileReader(vcfTest);){ 
			QualifiedSequenceList seqNames = genome.getSequencesMetadata();
			inGS.setSequences(seqNames);
			inTest.setSequences(seqNames);
			
			inGS.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			if(outVCF!=null) {
				out = new PrintStream(outVCF);
				writer.printHeader(inGS.getHeader(), out);
			}
			Iterator<VCFRecord> itGS = inGS.iterator();
			Iterator<VCFRecord> itTest = inTest.iterator();
			
			VCFRecord rGS = null;
			VCFRecord rTest = null;
			if(itGS.hasNext()) rGS = itGS.next();
			if(itTest.hasNext()) rTest = itTest.next();
			while (rGS!=null && rTest!=null) {
				CalledGenomicVariant callGS = rGS.getCalls().get(0);
				CalledGenomicVariant callTest = rTest.getCalls().get(0);
				if(callGS.isUndecided() ) {
					rGS = safeNext(itGS);
					continue;
				}
				if(callTest.isUndecided() ) {
					rTest = safeNext(itTest);
					continue;
				}
				int cmp = comparator.compare(rGS.getVariant(), rTest.getVariant());
				byte typeGS = loadType(rGS);
				byte typeTest = loadType(rTest);
				
				short qualTest = callTest.getGenotypeQuality();
				int n1 = getGenotypeNumber(callGS);
				int n2 = getGenotypeNumber(callTest);
				if(cmp<0) {
					//Gold standard variant not found
					countsPerType.get(typeGS).update(0,lastRowCounts,9+n1);
					//if(typeGS==GenomicVariant.TYPE_INDEL) System.out.println("Printing variant");
					if(mode == 1) writer.printVCFRecord(rGS, out); 
					rGS = safeNext(itGS);
				} else if (cmp>0) {
					//False positive variant
					countsPerType.get(typeTest).update(0,Math.min(qualTest/10, lastRowCounts),12+n2);
					if(mode == 2 && qualTest>=minQuality) writer.printVCFRecord(rTest, out);
					rTest = safeNext(itTest);
				} else {
					int k = 3*n1+n2;
					int row = Math.min(qualTest/10, lastRowCounts);
					GoldStandardComparisonCounts counts = countsPerType.get(typeGS); 
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
					rGS = safeNext(itGS);
					rTest = safeNext(itTest);
				}
			}
			while(rGS!=null) {
				//Gold standard variant not found
				CalledGenomicVariant callGS = rGS.getCalls().get(0);
				if(!callGS.isUndecided()) {
					byte typeGS = loadType(rGS);
					int n1 = getGenotypeNumber(callGS);
					countsPerType.get(typeGS).update(0,lastRowCounts,9+n1);
					if(mode == 1) writer.printVCFRecord(rGS, out);
				}
				rGS = safeNext(itGS);
			}
			while(rTest!=null) {
				//False positive variant
				CalledGenomicVariant callTest = rTest.getCalls().get(0);
				if(!callTest.isUndecided()) {
					byte typeTest = loadType(rTest);
					short qualTest = callTest.getGenotypeQuality();
					int n2 = getGenotypeNumber(callTest);
					countsPerType.get(typeTest).update(0,Math.min(qualTest/10, lastRowCounts),12+n2);
					if(mode == 2 && qualTest>=minQuality) writer.printVCFRecord(rTest, out);
				}
				rTest = safeNext(itTest);
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

	
	private VCFRecord safeNext(Iterator<VCFRecord> itGS) {
		VCFRecord rGS;
		if(itGS.hasNext()) rGS = itGS.next();
		else rGS = null;
		return rGS;
	}
	
	private byte loadType(VCFRecord record) {
		GenomicVariant v = record.getVariant();
		byte type = v.getType();
		if(type == GenomicVariant.TYPE_BIALLELIC_SNV) return type;
		if(type == GenomicVariant.TYPE_INDEL) return type;
		if(type == GenomicVariant.TYPE_STR) return type;
		if(v instanceof SNV) return GenomicVariant.TYPE_BIALLELIC_SNV;
		String [] alleles = v.getAlleles();
		if(alleles.length==2 && alleles[0].length()!=alleles[1].length()) return GenomicVariant.TYPE_INDEL;
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
