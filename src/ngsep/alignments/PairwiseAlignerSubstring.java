/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.alignments;

import ngsep.sequences.LimitedSequence;

public class PairwiseAlignerSubstring implements PairwiseAligner {

	@Override
	public PairwiseAlignment calculateAlignment(CharSequence sequence1, CharSequence sequence2) {
		int n1 = sequence1.length();
		int n2 = sequence2.length();
		if(n1<n2) {
			PairwiseAlignment alnRev = calculateAlignment(sequence2, sequence1);
			if(alnRev == null) return null;
			PairwiseAlignment aln = new PairwiseAlignment(sequence1, sequence2);
			aln.setScore(alnRev.getScore());
			aln.setStartLimits(alnRev.getStart2(), alnRev.getStart1());
			aln.setEndLimits(alnRev.getEnd2(), alnRev.getEnd1());
			aln.setAlignedSequences(alnRev.getAlignedSequence2(), alnRev.getAlignedSequence1());
			return aln;
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
		PairwiseAlignment aln = new PairwiseAlignment(sequence1, sequence2);
		aln.setAlignedSequences(strSeq1, aln2.toString());
		return aln;
	}
}
