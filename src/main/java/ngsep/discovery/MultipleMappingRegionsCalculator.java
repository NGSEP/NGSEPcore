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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.math.PhredScoreHelper;
import ngsep.variants.CalledCNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;


public class MultipleMappingRegionsCalculator {
	
	public static final String SOURCE_MULTIPLE_ALNS = "MultiAlns";
	private int minMQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	/**
	 * @return the minMQ
	 */
	public int getMinMQ() {
		return minMQ;
	}
	/**
	 * @param minMQ the minMQ to set
	 */
	public void setMinMQ(int minMQ) {
		this.minMQ = minMQ;
	}
	
	public List<CalledCNV> calculateMultipleMappingRegions(String alnsFile) throws IOException {
		List<CalledCNV> multipleMappingRegions = new ArrayList<CalledCNV>();
		GenomicRegionImpl lastRegion = null;
		int nonUniqueLastRegion = 0;
		int minReadLength=-1;
		LinkedList<Integer> uniqueStarts = new LinkedList<Integer>();
		
		try (ReadAlignmentFileReader reader = new ReadAlignmentFileReader(alnsFile);) {
			reader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_ALIGNMENT);
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			reader.setFilterFlags(filterFlags);
			reader.setMinMQ(minMQ);
			String currentSeqName = null;
			Iterator<ReadAlignment> it = reader.iterator();
			while(it.hasNext()) {
				ReadAlignment aln = it.next();
				if(aln.isPartialAlignment(10)) continue;
				//if(aln.getReadLength()<100) System.out.println("Small read alignment: "+aln.getReadName()+" length: "+aln.getReadLength()+" CIGAR: "+aln.getCigarString());
				if(minReadLength==-1 || minReadLength>aln.getReadLength()) minReadLength = aln.getReadLength();
				boolean sequenceChange = !aln.getSequenceName().equals(currentSeqName);
				if(lastRegion!=null && (sequenceChange || lastRegion.getLast() < aln.getFirst()-5)) {
					CalledCNV cnv = makeCNVCall(lastRegion, nonUniqueLastRegion, uniqueStarts, minReadLength);
					if(cnv!=null) multipleMappingRegions.add(cnv);
					lastRegion = null;
				}
				if(sequenceChange) {
					uniqueStarts.clear();
					currentSeqName = aln.getSequenceName();
				}
				else if (lastRegion==null && uniqueStarts.size()>100000) purgeList(uniqueStarts, aln.getFirst());
				boolean isUnique = aln.isUnique();
				if(!isUnique) {
					if(lastRegion == null) {
						lastRegion = new GenomicRegionImpl(aln.getSequenceName(), aln.getFirst(), aln.getLast());
						nonUniqueLastRegion=1;
					} else {
						nonUniqueLastRegion++;
						if (lastRegion.getLast()<aln.getLast()) lastRegion.setLast(aln.getLast());
					}
					
				} else {
					uniqueStarts.add(aln.getFirst());
				}
				
			}
		}
		
		if(lastRegion!=null) {
			CalledCNV cnv = makeCNVCall(lastRegion, nonUniqueLastRegion, uniqueStarts, minReadLength);
			if(cnv!=null) multipleMappingRegions.add(cnv);
		}
		return multipleMappingRegions;
	}
	public CalledCNV makeCNVCall(GenomicRegion region, int nonUniqueAlns, LinkedList<Integer> uniqueStarts, int minReadLength) {
		CalledCNV cnv = new CalledCNV(new GenomicVariantImpl(region.getSequenceName(), region.getFirst(), region.getLast(), GenomicVariant.TYPE_REPEAT));
		cnv.setSource(SOURCE_MULTIPLE_ALNS);
		cnv.setNonUniqueAlns(nonUniqueAlns);
		int uniqueAlns = 0;
		while(uniqueStarts.size()>0) {
			int nextFirst = uniqueStarts.peekFirst();
			int nextLast = nextFirst+minReadLength-1;
			//if(region.getFirst() == 2195287) System.out.println("Non unique alns: "+nonUniqueAlns+" Next first unique: "+nextFirst+" nextLast unique: "+nextLast+" readlength: "+minReadLength);
			if(nextFirst<=region.getLast()) {
				if(nextFirst>=region.getFirst() && nextLast<=region.getLast()) uniqueAlns++;
				uniqueStarts.removeFirst();
			} else {
				break;
			}
		}
		cnv.setUniqueAlns(uniqueAlns);
		//TODO: Make real p-value
		if(nonUniqueAlns<5) return null;
		double pValue = 1.0/(1.0+nonUniqueAlns);
		if(uniqueAlns>0) pValue = (double)uniqueAlns/(nonUniqueAlns+uniqueAlns);
		cnv.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(pValue));
		return cnv;
	}
	private void purgeList(LinkedList<Integer> uniqueMidPoints, int end) {
		while(uniqueMidPoints.size()>0) {
			int next = uniqueMidPoints.peekFirst();
			if(next<end) {
				uniqueMidPoints.removeFirst();
			} else {
				break;
			}
		}
	}
}
