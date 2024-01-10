package ngsep.alignments;

import ngsep.sequences.LimitedSequence;

public class PairwiseAlignerSubstring implements PairwiseAligner {

	@Override
	public String[] calculateAlignment(CharSequence sequence1, CharSequence sequence2) {
		int n1 = sequence1.length();
		int n2 = sequence2.length();
		if(n1<n2) {
			String[] alnRev = calculateAlignment(sequence2, sequence1);
			if(alnRev == null) return null;
			String [] answer = {alnRev[1].toString(),alnRev[0].toString()};
			return answer;
		}
		String strSeq1 = sequence1.toString();
		String strSeq2 = sequence2.toString();
		int pos = strSeq1.indexOf(strSeq2);
		if(pos==-1) return null;
		StringBuilder aln2 = new StringBuilder();
		int i;
		for(i=0;i<pos;i++) {
			aln2.append(LimitedSequence.GAP_CHARACTER);	
		}
		aln2.append(strSeq2);
		for(i=aln2.length();i<n1;i++) {
			aln2.append(LimitedSequence.GAP_CHARACTER);	
		}
		String [] answer = {strSeq1,aln2.toString()};
		return answer;
	}

}
