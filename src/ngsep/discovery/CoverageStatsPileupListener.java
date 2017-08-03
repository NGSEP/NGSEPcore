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
package ngsep.discovery;

import java.io.PrintStream;
import java.util.Arrays;

import ngsep.sequences.QualifiedSequence;


public class CoverageStatsPileupListener implements PileupListener {
	private int maxCoverage = 300;
	private int [] coverageCounts;
	private int highCoverageCount = 0;
	private int [] coverageCountUniqueAlignments;
	private int highCoverageCountUniqueAlignments = 0;
	
	public CoverageStatsPileupListener () {
		reset();
	}
	public void reset() {
		this.coverageCounts = new int[maxCoverage];
		Arrays.fill(coverageCounts, 0);
		highCoverageCount = 0;
		this.coverageCountUniqueAlignments = new int[maxCoverage];
		Arrays.fill(coverageCountUniqueAlignments, 0);
		highCoverageCountUniqueAlignments = 0;
	}
	public int[] getCoverageCounts() {
		return coverageCounts;
	}
	public int getMaxCoverage() {
		return maxCoverage;
	}
	public void setMaxCoverage(int maxCoverage) {
		this.maxCoverage = maxCoverage;
		reset();
	}


	public int getHighCoverageCount() {
		return highCoverageCount;
	}

	public void addPileup (PileupRecord record) {
		int calls = record.getNumAlignments();
		if(calls < coverageCounts.length) {
			coverageCounts[calls]++;
		} else {
			highCoverageCount++;
		}
		int callsUniqueAlns = record.getNumUniqueAlns();
		if(callsUniqueAlns < coverageCountUniqueAlignments.length) {
			coverageCountUniqueAlignments[callsUniqueAlns]++;
		} else {
			highCoverageCountUniqueAlignments++;
		}
	}
	
	public int getCoverageMaxCount() {
		int maxIndex = 1;
		for(int i=1;i<coverageCounts.length;i++) {
			if(coverageCounts[maxIndex]<coverageCounts[i]) {
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	public void printCoverageStats(PrintStream out) {
		for(int i=1;i<coverageCounts.length;i++) {
			out.println(""+i+"\t"+coverageCounts[i]+"\t"+coverageCountUniqueAlignments[i]);
		}
		out.println("More\t"+highCoverageCount+"\t"+highCoverageCountUniqueAlignments);
		out.flush();
	}

	@Override
	public void onPileup(PileupRecord pileup) {
		addPileup(pileup);
		
	}

	
	@Override
	public void onSequenceStart(QualifiedSequence sequenceName) {
		
	}
	

	@Override
	public void onSequenceEnd(QualifiedSequence sequenceName) {
		
	}

}
