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

import java.util.Arrays;

import ngsep.sequences.LimitedSequence;

class PairwiseAlignerNaive implements PairwiseAligner {
	private boolean gapsLeft = false;
	
	public PairwiseAlignerNaive(boolean gapsLeft) {
		super();
		this.gapsLeft = gapsLeft;
	}


	@Override
	public PairwiseAlignment calculateAlignment(CharSequence sequence1, CharSequence sequence2) {
		int diff = sequence1.length()- sequence2.length();
		
		StringBuilder aln1 = new StringBuilder();
		StringBuilder aln2 = new StringBuilder();
		if(gapsLeft) appendGaps(diff, aln1, aln2);
		aln1.append(sequence1);
		aln2.append(sequence2);
		if(!gapsLeft) appendGaps(diff, aln1, aln2);
		PairwiseAlignment aln = new PairwiseAlignment(sequence1, sequence2);
		aln.setAlignedSequences(aln1.toString(), aln2.toString());
		return aln;
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