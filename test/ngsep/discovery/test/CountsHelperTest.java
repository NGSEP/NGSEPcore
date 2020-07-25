package ngsep.discovery.test;

import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;
import ngsep.discovery.CountsHelper;
import ngsep.discovery.PileupAlleleCall;
import ngsep.sequences.DNASequence;

public class CountsHelperTest extends TestCase {
	public void testConditionals() {
		System.out.println("Helper 0.01");
		String [] alleles = new String[2];
		alleles[0] = DNASequence.BASES_ARRAY[0];
		alleles[1] = DNASequence.BASES_ARRAY[1];
		CountsHelper helper = new CountsHelper(alleles);
		int minorCount = 10;
		int totalCount = 1000;
		byte maxQS = 30;
		helper.setMaxBaseQS(maxQS);
		helper.setHeterozygousProportion((double)minorCount/totalCount);
		List<PileupAlleleCall> calls = createSNVAlleleCalls(totalCount, minorCount);
		for(PileupAlleleCall call:calls ) {
			byte q = (byte)(Math.min(CountsHelper.DEF_MAX_BASE_QS, call.getQualityScores().charAt(0)-33));
			helper.updateCounts(call.getAlleleString().substring(0,1), q, call.isNegativeStrand());
			System.out.println(helper.getLogConditionalProbs()[0][0]+" "+helper.getLogConditionalProbs()[0][1]);
		}
		System.out.println("Helper 0.25");
		CountsHelper helper2 = new CountsHelper(alleles);
		helper2.setMaxBaseQS(maxQS);
		helper2.setHeterozygousProportion(0.25);
		for(PileupAlleleCall call:calls ) {
			byte q = (byte)(Math.min(CountsHelper.DEF_MAX_BASE_QS, call.getQualityScores().charAt(0)-33));
			helper2.updateCounts(call.getAlleleString().substring(0,1), q, call.isNegativeStrand());
			System.out.println(helper2.getLogConditionalProbs()[0][0]+" "+helper2.getLogConditionalProbs()[0][1]);
		}
		System.out.println("Helper 0.5");
		CountsHelper helper3 = new CountsHelper(alleles);
		helper3.setMaxBaseQS(maxQS);
		helper3.setHeterozygousProportion(0.5);
		for(PileupAlleleCall call:calls ) {
			byte q = (byte)(Math.min(CountsHelper.DEF_MAX_BASE_QS, call.getQualityScores().charAt(0)-33));
			helper3.updateCounts(call.getAlleleString().substring(0,1), q, call.isNegativeStrand());
			System.out.println(helper3.getLogConditionalProbs()[0][0]+" "+helper3.getLogConditionalProbs()[0][1]);
		}
		System.out.println("Summary");
		System.out.println(helper.getLogConditionalProbs()[0][0]+" "+helper.getLogConditionalProbs()[0][1]+" "+helper.getLogConditionalProbs()[1][0]);
		System.out.println(helper2.getLogConditionalProbs()[0][0]+" "+helper2.getLogConditionalProbs()[0][1]+" "+helper2.getLogConditionalProbs()[1][0]);
		System.out.println(helper3.getLogConditionalProbs()[0][0]+" "+helper3.getLogConditionalProbs()[0][1]+" "+helper3.getLogConditionalProbs()[1][0]);
	}
	
	public void testConditionalsStaticMethods() {
		System.out.println("Helper 0.01");
		String [] alleles = new String[2];
		alleles[0] = DNASequence.BASES_ARRAY[0];
		alleles[1] = DNASequence.BASES_ARRAY[1];
		int minorCount = 10;
		int totalCount = 1000;
		List<PileupAlleleCall> calls = createSNVAlleleCalls(totalCount, minorCount);
		CountsHelper helper = CountsHelper.calculateCountsGTSNV(alleles, calls, (byte)30, (double)minorCount/totalCount);
		System.out.println("Helper 0.25");
		CountsHelper helper2 = CountsHelper.calculateCountsGTSNV(alleles, calls, (byte)30, 0.25);
		System.out.println("Helper 0.5");
		CountsHelper helper3 = CountsHelper.calculateCountsGTSNV(alleles, calls, (byte)30, 0.5);
		System.out.println("Summary");
		System.out.println(helper.getLogConditionalProbs()[0][0]+" "+helper.getLogConditionalProbs()[0][1]+" "+helper.getLogConditionalProbs()[1][0]);
		System.out.println(helper2.getLogConditionalProbs()[0][0]+" "+helper2.getLogConditionalProbs()[0][1]+" "+helper2.getLogConditionalProbs()[1][0]);
		System.out.println(helper3.getLogConditionalProbs()[0][0]+" "+helper3.getLogConditionalProbs()[0][1]+" "+helper3.getLogConditionalProbs()[1][0]);
	}

	private List<PileupAlleleCall> createSNVAlleleCalls(int totalCount, int minorCount) {
		List<PileupAlleleCall> calls = new ArrayList<PileupAlleleCall>();
		int i=0;
		while(i<totalCount-minorCount) {
			PileupAlleleCall call = new PileupAlleleCall(DNASequence.BASES_ARRAY[0], "A");
			calls.add(call);
			i++;
		}
		while(i<totalCount) {
			PileupAlleleCall call = new PileupAlleleCall(DNASequence.BASES_ARRAY[1], "A");
			calls.add(call);
			i++;
		}
		return calls;
	}
}
