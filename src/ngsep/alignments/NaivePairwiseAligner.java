package ngsep.alignments;

import java.util.Arrays;

import ngsep.sequences.LimitedSequence;
import ngsep.sequences.PairwiseAligner;

class NaivePairwiseAligner implements PairwiseAligner {
	private boolean gapsLeft = false;
	
	public NaivePairwiseAligner(boolean gapsLeft) {
		super();
		this.gapsLeft = gapsLeft;
	}


	@Override
	public String[] calculateAlignment(CharSequence sequence1, CharSequence sequence2) {
		int diff = sequence1.length()- sequence2.length();
		
		StringBuilder aln1 = new StringBuilder();
		StringBuilder aln2 = new StringBuilder();
		if(gapsLeft) appendGaps(diff, aln1, aln2);
		aln1.append(sequence1);
		aln2.append(sequence2);
		if(!gapsLeft) appendGaps(diff, aln1, aln2);
		
		String [] answer = {aln1.toString(),aln2.toString()};
		return answer;
	}


	private void appendGaps(int diff, StringBuilder aln1, StringBuilder aln2) {
		char [] gaps = new char [Math.abs(diff)];
		Arrays.fill(gaps, LimitedSequence.GAP_CHARACTER);
		if(diff>0) {
			aln2.append(gaps);
		} else if(diff < 0){
			aln1.append(gaps);
		}
	}
	
}