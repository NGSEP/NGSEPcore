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

import ngsep.alignments.ReadAlignment.Platform;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.ReferenceGenomeFMIndex;

public class ReadAlignmentObjectsFactory {
	
	private Logger log = Logger.getAnonymousLogger();
	private int kmerLength = ReadsAligner.DEF_KMER_LENGTH;
	private int windowLength = ReadsAligner.DEF_WINDOW_LENGTH;
	private int numThreads = 1;
	private Platform platform;
	private int alignmentAlgorithm=UngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS;
	private ReferenceGenome genome;
	private ReferenceGenomeFMIndex fmIndex;

	
	public ReadAlignmentObjectsFactory(ReferenceGenome genome) {
		super();
		this.genome = genome;
	}
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
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
	
	public Platform getPlatform() {
		return platform;
	}
	public void setPlatform(Platform platform) {
		this.platform = platform;
	}
	public int getAlignmentAlgorithm() {
		return alignmentAlgorithm;
	}
	public void setAlignmentAlgorithm(int alignmentAlgorithm) {
		this.alignmentAlgorithm = alignmentAlgorithm;
	}
	
	public ReferenceGenomeFMIndex getFmIndex() {
		return fmIndex;
	}
	public void setFmIndex(ReferenceGenomeFMIndex fmIndex) {
		this.fmIndex = fmIndex;
	}


	UngappedSearchHitsClustersFinder first=null;
	public synchronized UngappedSearchHitsClustersFinder requestClustersFinder()  {
		if(first == null) {
			createFirstFinder();
			return first;
		} else if (first instanceof MinimizersUngappedSearchHitsClustersFinder){
			MinimizersUngappedSearchHitsClustersFinder firstI = (MinimizersUngappedSearchHitsClustersFinder) first;
			MinimizersUngappedSearchHitsClustersFinder next = new MinimizersUngappedSearchHitsClustersFinder();
			next.setLog(log);
			next.setMinRawHits(firstI.getMinRawHits());
			next.setMinProportionReadLength(firstI.getMinProportionReadLength());
			next.setKmerCodesTable(genome, firstI.getKmerCodesTable());
			return next;
		} else {
			FMIndexUngappedSearchHitsClustersFinder next = new FMIndexUngappedSearchHitsClustersFinder(fmIndex, kmerLength);
			return next;
		}
	}
	private void createFirstFinder() {
		Runtime runtime = Runtime.getRuntime();
		long startTime = System.currentTimeMillis();
		if(platform.isLongReads()) {
			MinimizersUngappedSearchHitsClustersFinder finder = new MinimizersUngappedSearchHitsClustersFinder();
			finder.setLog(log);
			//if(!platform.isLongReads()) first.setMinRawHits(1);
			finder.loadGenome (genome, kmerLength, windowLength, numThreads,false);
			first = finder; 
		} else {
			if (fmIndex!=null) {
				log.info("Aligning reads using built index with "+fmIndex.getSequencesMetadata().size()+" sequences");
			} else {
				log.info("Calculating FM-index from genome file: "+genome.getFilename());
				fmIndex = new ReferenceGenomeFMIndex(genome, log);
			}
			FMIndexUngappedSearchHitsClustersFinder finder = new FMIndexUngappedSearchHitsClustersFinder(fmIndex, kmerLength);
			first = finder;
		}
		
		long usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time2 = System.currentTimeMillis();
		long diff = (time2-startTime)/1000;
		if(diff>5) log.info("Created first clusters finder. Time (s): "+diff+". Memory: "+usedMemory);
	}
	private List<LongReadsUngappedSearchHitsClusterAligner> aligners = new ArrayList<LongReadsUngappedSearchHitsClusterAligner>();
	private int lastAlignerIndex = 0;
	public synchronized UngappedSearchHitsClusterAligner requestAligner()  {
		if(UngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_SHORT_READS==alignmentAlgorithm) return new ShortReadsUngappedSearchHitsClusterAligner();
		if(aligners.size()==0) {
			LongReadsUngappedSearchHitsClusterAligner aligner = new LongReadsUngappedSearchHitsClusterAligner(alignmentAlgorithm);
			aligners.add(aligner);
			return aligner;
		} else if (alignmentAlgorithm!=UngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_AFFINE_GAP || aligners.size()<2*numThreads) {
			LongReadsUngappedSearchHitsClusterAligner aligner = new LongReadsUngappedSearchHitsClusterAligner(alignmentAlgorithm);
			aligners.add(aligner);
			lastAlignerIndex=aligners.size()-1;
			return aligner;
		}
		lastAlignerIndex++;
		if(lastAlignerIndex==aligners.size()) lastAlignerIndex=0;
		return aligners.get(lastAlignerIndex);
	}
}
