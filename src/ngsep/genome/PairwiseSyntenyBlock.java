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

import java.util.List;

/**
 * 
 * @author Laura Gonzalez
 *
 */
public class PairwiseSyntenyBlock {

	private int genomeId1;
	private GenomicRegion regionGenome1;
	private int genomeId2;
	private GenomicRegion regionGenome2;
	private List<SyntenyVertex> homologies;
	
	
	
	public PairwiseSyntenyBlock(List<SyntenyVertex> homologies) {
		String seqName1 = null;
		int first1 = Integer.MAX_VALUE;
		int last1 = 0;
		String seqName2 = null;
		int first2 = Integer.MAX_VALUE;
		int last2 = 0;
		int lastStart = -1;
		int votesPositive = 0;
		int votesNegative = 0;
		
		for (SyntenyVertex sv : homologies) {
			LocalHomologyCluster g1 = sv.getLocalRegion1();
			LocalHomologyCluster g2 = sv.getLocalRegion2();
			if(seqName1==null) {
				genomeId1 = g1.getGenomeId();
				seqName1 = g1.getSequenceName();
			}
			first1 = Math.min(first1, g1.getFirst());
			first2 = Math.min(first2, g2.getFirst());
			if(seqName2==null) {
				genomeId2 = g2.getGenomeId();
				seqName2 = g2.getSequenceName();
			}
			last1 = Math.max(last1, g1.getLast());
			last2 = Math.max(last2, g2.getLast());
			if(lastStart>0) {
				if( g2.getFirst()>=lastStart) votesPositive++;
				else votesNegative++;
			}
			//if("chrVII".equals(seqName1)) System.out.println("Next vertex "+g1.getSequenceName()+":"+g1.getFirst()+"-"+g1.getLast()+" to "+g2.getSequenceName()+":"+g2.getFirst()+"-"+g2.getLast()+" lastStart: "+lastStart+" votes: "+votesPositive+" "+votesNegative);
			lastStart = g2.getFirst();
		}
		regionGenome1 = new GenomicRegionImpl(seqName1, first1, last1);
		regionGenome2 = new GenomicRegionImpl(seqName2, first2, last2);
		((GenomicRegionImpl)regionGenome2).setNegativeStrand(votesNegative>votesPositive);
		this.homologies = homologies;
	}

	public List<SyntenyVertex> getHomologies() {
		return homologies;
	}

	
	public int getGenomeId1() {
		return genomeId1;
	}

	public int getGenomeId2() {
		return genomeId2;
	}

	public GenomicRegion getRegionGenome1() {
		return regionGenome1;
	}

	public GenomicRegion getRegionGenome2() {
		return regionGenome2;
	}
}
