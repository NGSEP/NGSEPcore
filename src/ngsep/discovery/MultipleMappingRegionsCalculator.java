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
	private short minGenotypeQuality = 7;
	private int minMQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	
	
	
	public short getMinGenotypeQuality() {
		return minGenotypeQuality;
	}
	public void setMinGenotypeQuality(short minGenotypeQuality) {
		this.minGenotypeQuality = minGenotypeQuality;
	}
	
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
		LinkedList<Integer> uniqueMidPoints = new LinkedList<Integer>();
		ReadAlignmentFileReader reader = null;
		try {
			reader = new ReadAlignmentFileReader(alnsFile);
			reader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_MINIMAL);
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			reader.setFilterFlags(filterFlags);
			reader.setMinMQ(minMQ);
			String currentSeqName = null;
			Iterator<ReadAlignment> it = reader.iterator();
			while(it.hasNext()) {
				ReadAlignment aln = it.next();
				boolean sequenceChange = !aln.getSequenceName().equals(currentSeqName);
				if(lastRegion!=null && (sequenceChange || lastRegion.getLast() < aln.getFirst()-5)) {
					CalledCNV cnv = makeCNVCall(lastRegion, nonUniqueLastRegion, uniqueMidPoints);
					if (cnv.getGenotypeQuality()>=minGenotypeQuality) multipleMappingRegions.add(cnv);
					lastRegion = null;
				}
				if(sequenceChange) {
					uniqueMidPoints.clear();
					currentSeqName = aln.getSequenceName();
				}
				else if (lastRegion==null && uniqueMidPoints.size()>100000) purgeList(uniqueMidPoints, aln.getFirst());
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
					int midPoint = aln.getFirst()+aln.getReadLength()/2;
					uniqueMidPoints.add(midPoint);
				}
				
			}
		} finally {
			if (reader!=null) reader.close();
		}
		
		if(lastRegion!=null) {
			CalledCNV cnv = makeCNVCall(lastRegion, nonUniqueLastRegion, uniqueMidPoints);
			if (cnv.getGenotypeQuality()>=minGenotypeQuality) multipleMappingRegions.add(cnv);
		}
		return multipleMappingRegions;
	}
	public CalledCNV makeCNVCall(GenomicRegion region, int nonUniqueAlns, LinkedList<Integer> uniqueMidPoints) {
		CalledCNV cnv = new CalledCNV(new GenomicVariantImpl(region.getSequenceName(), region.getFirst(), region.getLast(), GenomicVariant.TYPE_REPEAT));
		cnv.setSource(SOURCE_MULTIPLE_ALNS);
		cnv.setNonUniqueAlns(nonUniqueAlns);
		int uniqueAlns = 0;
		while(uniqueMidPoints.size()>0) {
			int next = uniqueMidPoints.peekFirst();
			if(next<=region.getLast()) {
				if(next>=region.getFirst())uniqueAlns++;
				uniqueMidPoints.removeFirst();
			} else {
				break;
			}
		}
		cnv.setUniqueAlns(uniqueAlns);
		//TODO: Improve p-value
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
