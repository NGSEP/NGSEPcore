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
package ngsep.genome;

import java.io.PrintStream;
import java.util.List;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.SimpleSequenceListLoader;

/**
 * Script that compares sets of genomic regions in the same reference genome
 * @author Jorge Duitama
 *
 */

public class GenomicRegionSetsComparator {
	private double maxPCT = 1000;
	/**
	 * compares sets of genomic regions in the same reference genome.
	 * Files to compare should have at least three columns: chromosome, first position and last position
	 * Coordinates should be 1-based
	 * @param args Usage parameters.
	 * args[0] text file with chromosome names. The fai file can be used
	 * args[1] first text file with genomic regions to compare
	 * args[2] second text file with genomic regions to compare
	 * args[3] output prefix. One file will be generated for each dataset
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		GenomicRegionSetsComparator instance = new GenomicRegionSetsComparator();
		SimpleSequenceListLoader seqNamesHandler = new SimpleSequenceListLoader();
		QualifiedSequenceList seqNames = seqNamesHandler.loadSequences(args[0]);
		String filename1 = args[1];
		String filename2 = args[2];
		SimpleGenomicRegionFileHandler fh = new SimpleGenomicRegionFileHandler();
		GenomicRegionSortedCollection<GenomicRegion> regionSet1 = new GenomicRegionSortedCollection<GenomicRegion>(seqNames);
		regionSet1.addAll(fh.loadRegions(filename1));
		GenomicRegionSortedCollection<GenomicRegion> regionSet2 = new GenomicRegionSortedCollection<GenomicRegion>(seqNames);
		regionSet2.addAll(fh.loadRegions(filename2));
		String outPrefix = args[3];
		if(args.length>4) instance.maxPCT = Double.parseDouble(args[4]);
		int [] stats1;
		try (PrintStream out1 = new PrintStream(outPrefix+"_overlapSet1.txt")) {
			stats1 = instance.calculateOverlap(regionSet1.asList(),regionSet2,out1);
		}
		int [] stats2;
		try (PrintStream out2 = new PrintStream(outPrefix+"_overlapSet2.txt")) {
			stats2 = instance.calculateOverlap(regionSet2.asList(),regionSet1,out2);
		}
		System.out.println("Total set 1: "+stats1[0]);
		System.out.println("Covered set 1: "+stats1[1]);
		System.out.println("PCT covered set 1: "+(100.0*stats1[1]/stats1[0]));
		System.out.println("Total set 2: "+stats2[0]);
		System.out.println("Covered set 2: "+stats2[1]);
		System.out.println("PCT covered set 2: "+(100.0*stats2[1]/stats2[0]));
		
	}

	private int[] calculateOverlap(List<GenomicRegion> statsRegions, GenomicRegionSortedCollection<GenomicRegion> queryRegions, PrintStream out) {
		int totalStats = 0;
		int coveredStats = 0;
		for(GenomicRegion r:statsRegions) {
			totalStats+=r.length();
			GenomicRegionSortedCollection<GenomicRegion> overlapSet = queryRegions.findSpanningRegions(r);
			GenomicRegionSpanComparator spanCMP = GenomicRegionSpanComparator.getInstance();
			
			int maxOverlap = 0;
			double sumOverlap = 0;
			int lengthRegMaxOverlap = 0;
			for(GenomicRegion r2:overlapSet) {
				int overlap = spanCMP.getSpanLength(r.getFirst(), r.getLast(), r2.getFirst(), r2.getLast());
				sumOverlap+=overlap;
				if(overlap > maxOverlap) {
					maxOverlap = overlap;
					lengthRegMaxOverlap = r2.length();
				}
			}
			coveredStats += maxOverlap;
			double pct1 = 100.0*maxOverlap/r.length();
			double pct2 = 0;
			if(lengthRegMaxOverlap>0) pct2 = 100.0*maxOverlap/lengthRegMaxOverlap;
			double pct3 = 100.0*sumOverlap/r.length();
			if(pct1<=maxPCT && pct3<=maxPCT) {
				out.println(r.getSequenceName()+"\t"+r.getFirst()+"\t"+r.getLast()+"\t"+r.length()+"\t"+lengthRegMaxOverlap+"\t"+maxOverlap+"\t"+pct1+"\t"+pct2+"\t"+pct3);
			}
		}
		int [] answer = {totalStats,coveredStats};
		return answer;
	}
}
