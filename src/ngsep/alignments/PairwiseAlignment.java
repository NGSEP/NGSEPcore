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

public class PairwiseAlignment {
	private CharSequence sequence1;
	private CharSequence sequence2;
	private String alignedSequence1;
	private String alignedSequence2;
	private int start1 = 0;
	private int start2 = 0;
	private int end1 = 0;
	private int end2 = 0;
	private int score = 0;
	private int mismatches = 0;

	
	public PairwiseAlignment(CharSequence sequence1, CharSequence sequence2) {
		super();
		this.sequence1 = sequence1;
		this.sequence2 = sequence2;
		end1 = sequence1.length();
		end2 = sequence2.length();
	}
	
	public CharSequence getSequence1() {
		return sequence1;
	}

	public CharSequence getSequence2() {
		return sequence2;
	}

	public int getScore() {
		return score;
	}
	public void setScore(int score) {
		this.score = score;
	}
	public String getAlignedSequence1() {
		return alignedSequence1;
	}
	public String getAlignedSequence2() {
		return alignedSequence2;
	}
	public int getStart1() {
		return start1;
	}
	public int getStart2() {
		return start2;
	}
	public int getEnd1() {
		return end1;
	}
	public int getEnd2() {
		return end2;
	}
	public int getMismatches() {
		return mismatches;
	}

	public void setEndLimits(int end1, int end2) {
		this.end1 = end1;
		this.end2 = end2;
	}
	public void setStartLimits(int start1, int start2) {
		this.start1 = start1;
		this.start2 = start2;
	}
	public void setAlignedSequences (String seq1, String seq2) {
		this.alignedSequence1 = seq1;
		this.alignedSequence2 = seq2;
		mismatches = countMismatches();
	}
	private int countMismatches() {
		int answer = 0;
		int nextGapMismatches = 0;
		boolean started = false;
		for(int i=0;i<alignedSequence1.length();i++) {
			char c1 = alignedSequence1.charAt(i);
			char c2 = alignedSequence2.charAt(i);
			if(c1==LimitedSequence.GAP_CHARACTER || c2 == LimitedSequence.GAP_CHARACTER) {
				nextGapMismatches++;
				if(nextGapMismatches==1) nextGapMismatches++;
			} else {
				if(c1!=c2) answer++;
				if(started) answer+=nextGapMismatches;
				else started = true;
				nextGapMismatches=0;	
			}
		}
		return answer;
	}
}
