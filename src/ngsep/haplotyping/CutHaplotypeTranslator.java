/*******************************************************************************
 * SingleIndividualHaplotyper - Efficient heuristic algorithms for the SIH problem
 * Copyright 2011 Jorge Duitama
 *
 * This file is part of SingleIndividualHaplotyper.
 *
 *     SingleIndividualHaplotyper is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     SingleIndividualHaplotyper is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with SingleIndividualHaplotyper.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.haplotyping;

import ngsep.variants.CalledGenomicVariant;

public class CutHaplotypeTranslator {
	public static final int CONSENSUS_GROUP_1 = 1;
	public static final int CONSENSUS_GROUP_2 = 2;
	public static final int CONSENSUS_COMBINED = 3;
	public static final int CONSENSUS_ALL = 4;
	public static byte [] getHaplotype(HaplotypeBlock block, boolean[] cut, int consensusType) {
		int numFragments = block.getNumFragments();
		
		int numVariants = block.getNumVariants();
		byte [] haplotype = new byte [numVariants];
		
		int i=0;
		for(int j=0;j <numVariants;j++) {
			while(i<numFragments && block.getLastColumn(i)<j) {
				i++;
			}
			int votesAlleleRef = 0;
			int votes=0;
			for(int k=i; k<numFragments;k++) {
				
				if(block.getFirstColumn(k)>j) {
					break;
				}
				byte call = block.getAllele(k, j);
				if(call != CalledGenomicVariant.ALLELE_UNDECIDED) {
					if(consensusType == CONSENSUS_COMBINED) {
						votes++;
						if((call == CalledGenomicVariant.ALLELE_REFERENCE && !cut[k]) || (call != CalledGenomicVariant.ALLELE_REFERENCE && cut[k]) ) {
							votesAlleleRef++;
						}
					} else if (consensusType == CONSENSUS_GROUP_1 && !cut[k]) {
						votes++;
						if(call == CalledGenomicVariant.ALLELE_REFERENCE) {
							votesAlleleRef++;
						}
					} else if (consensusType == CONSENSUS_GROUP_2 && cut[k]) {
						votes++;
						if(call == CalledGenomicVariant.ALLELE_REFERENCE) {
							votesAlleleRef++;
						}
					} else if (consensusType == CONSENSUS_ALL) {
						votes++;
						if(call == CalledGenomicVariant.ALLELE_REFERENCE) {
							votesAlleleRef++;
						}
					}
				}
			}
			if (2*votesAlleleRef < votes) {
				haplotype[j] = CalledGenomicVariant.ALLELE_ALTERNATIVE;
			} else if (2*votesAlleleRef > votes){
				haplotype[j] = CalledGenomicVariant.ALLELE_REFERENCE;
			} else {
				haplotype[j] = CalledGenomicVariant.ALLELE_UNDECIDED;
			}
		}
		return haplotype;
	}
}
