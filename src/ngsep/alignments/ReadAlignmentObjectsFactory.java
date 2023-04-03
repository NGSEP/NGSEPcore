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

public class ReadAlignmentObjectsFactory {
	
	private Logger log = Logger.getAnonymousLogger();
	private int kmerLength = ReadsAligner.DEF_KMER_LENGTH;
	private int windowLength = ReadsAligner.DEF_WINDOW_LENGTH;
	private int numThreads = 1;
	private Platform platform;
	private int alignmentAlgorithm=UngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_AFFINE_GAP;
	private ReferenceGenome genome;

	
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
	MinimizersUngappedSearchHitsClustersFinder first=null;
	public synchronized UngappedSearchHitsClustersFinder requestClustersFinder()  {
		if(first == null) {
			createFirstFinder();
			return first;
		} else {
			MinimizersUngappedSearchHitsClustersFinder next = new MinimizersUngappedSearchHitsClustersFinder();
			next.setLog(log);
			next.setMinRawHits(first.getMinRawHits());
			next.setMinProportionReadLength(first.getMinProportionReadLength());
			next.setKmerCodesTable(genome, first.getKmerCodesTable());
			return next;
		}
	}
	private void createFirstFinder() {
		Runtime runtime = Runtime.getRuntime();
		long startTime = System.currentTimeMillis();
		first = new MinimizersUngappedSearchHitsClustersFinder();
		first.setLog(log);
		if(!platform.isLongReads()) first.setMinRawHits(2);
		first.loadGenome (genome, kmerLength, windowLength, numThreads,false);
		long usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time2 = System.currentTimeMillis();
		long diff = (time2-startTime)/1000;
		log.info("Created first clusters finder. Time (s): "+diff+". Memory: "+usedMemory);
		/*if (fMIndex!=null) {
						log.info("Aligning reads using built index with "+fMIndex.getSequencesMetadata().size()+" sequences");
					} else if (fmIndexFile!=null) {
						log.info("Loading reference index from file: "+fmIndexFile);
						fMIndex = ReferenceGenomeFMIndex.load(genome, fmIndexFile);
					} else {
						log.info("Calculating FM-index from genome file: "+genome.getFilename());
						fMIndex = new ReferenceGenomeFMIndex(genome, log);
					}
					createFMIndexReadsAligner();
				}*/
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
