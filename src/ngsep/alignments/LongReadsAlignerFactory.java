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

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;

public class LongReadsAlignerFactory {
	
	private Logger log = Logger.getAnonymousLogger();
	private int maxAlnsPerRead = ReadsAligner.DEF_MAX_ALNS_PER_READ;
	private int kmerLength = ReadsAligner.DEF_KMER_LENGTH;
	private int windowLength = ReadsAligner.DEF_WINDOW_LENGTH;
	private int numThreads = 1;
	private ReferenceGenome genome;

	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	public int getMaxAlnsPerRead() {
		return maxAlnsPerRead;
	}
	public void setMaxAlnsPerRead(int maxAlnsPerRead) {
		this.maxAlnsPerRead = maxAlnsPerRead;
	}
	public int getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
	}
	public int getWindowLength() {
		return windowLength;
	}
	public void setWindowLength(int windowLength) {
		this.windowLength = windowLength;
	}
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	private List<MinimizersTableReadAlignmentAlgorithm> longReadsAligners = new ArrayList<MinimizersTableReadAlignmentAlgorithm>();
	private int lastReadsAlignerIndex = 0;
	public synchronized MinimizersTableReadAlignmentAlgorithm requestLongReadsAligner()  {
		return requestLongReadsAligner(MinimizersTableReadAlignmentAlgorithm.ALIGNMENT_ALGORITHM_AFFINE_GAP);
	}
	public synchronized MinimizersTableReadAlignmentAlgorithm requestLongReadsAligner(int alignmentAlgorithm)  {
		if(longReadsAligners.size()==0) {
			createFirstLongReadAligner(alignmentAlgorithm);
			return longReadsAligners.get(0);
		} else if (alignmentAlgorithm!=MinimizersTableReadAlignmentAlgorithm.ALIGNMENT_ALGORITHM_AFFINE_GAP || longReadsAligners.size()<2*numThreads) {
			MinimizersTableReadAlignmentAlgorithm aligner = new MinimizersTableReadAlignmentAlgorithm(alignmentAlgorithm);
			MinimizersTableReadAlignmentAlgorithm first = longReadsAligners.get(0);
			aligner.setLog(log);
			aligner.setMaxAlnsPerRead(maxAlnsPerRead);
			if(genome!=null) aligner.setMinimizersTable(genome, first.getMinimizersTable());
			longReadsAligners.add(aligner);
			lastReadsAlignerIndex=longReadsAligners.size()-1;
			//Runtime runtime = Runtime.getRuntime();
			//long usedMemory = runtime.totalMemory()-runtime.freeMemory();
			//usedMemory/=1000000000;
			//System.out.println("Created long reads aligner "+(lastReadsAlignerIndex+1)+". Memory: "+usedMemory);
			return aligner;
		}
		lastReadsAlignerIndex++;
		if(lastReadsAlignerIndex==longReadsAligners.size()) lastReadsAlignerIndex=0;
		return longReadsAligners.get(lastReadsAlignerIndex);
	}
	private void createFirstLongReadAligner(int alignmentAlgorithm) {
		Runtime runtime = Runtime.getRuntime();
		long startTime = System.currentTimeMillis();
		MinimizersTableReadAlignmentAlgorithm longReadsAligner = new MinimizersTableReadAlignmentAlgorithm(alignmentAlgorithm);
		longReadsAligner.setLog(log);
		longReadsAligner.setMaxAlnsPerRead(maxAlnsPerRead);
		if(genome!=null) {
			longReadsAligner.loadGenome (genome, kmerLength, windowLength, numThreads,false);
		}
		long usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time2 = System.currentTimeMillis();
		long diff = (time2-startTime)/1000;
		log.info("Created first long reads aligner. Time (s): "+diff+". Memory: "+usedMemory);
		longReadsAligners.add(longReadsAligner);
	}
}
