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

import ngsep.variants.HaplotypeBlock;

public class CutHaplotypeTranslator {
	public static final int CONSENSUS_GROUP_1 = 1;
	public static final int CONSENSUS_GROUP_2 = 2;
	public static final int CONSENSUS_COMBINED = 3;
	public static final int CONSENSUS_ALL = 4;
	public static byte [] getHaplotype(HaplotypeBlock block, boolean[] cut, int consensusType) {
		int numFragments = block.getNumFragments();
		int numVariants = block.getNumVariants();
		byte [] haplotype = new byte [numVariants];
		/*int pos=b.getFirstPos();
		int lastPos = b.getLastPos();
		int i=0;
		for(;pos <=lastPos;pos++) {
			while(i<numFragments && fragments.get(i).getLastPos()<pos) {
				i++;
			}
			int votesAllele1 = 0;
			int votes=0;
			for(int j=i;j<numFragments;j++) {
				
				if(f2.getFirstPos()>pos) {
					break;
				}
				byte call = block.getAllele(j, pos);
				if(call != CalledGenomicVariant.ALLELE_UNDECIDED) {
					if(consensusType == CONSENSUS_COMBINED) {
						votes++;
						if((call == CalledGenomicVariant.ALLELE_REFERENCE && !cut[j]) || (call != CalledGenomicVariant.ALLELE_REFERENCE && cut[j]) ) {
							votesAllele1++;
						}
					} else if (consensusType == CONSENSUS_GROUP_1 && !cut[j]) {
						votes++;
						if(call == CalledGenomicVariant.ALLELE_REFERENCE) {
							votesAllele1++;
						}
					} else if (consensusType == CONSENSUS_GROUP_2 && cut[j]) {
						votes++;
						if(call == CalledGenomicVariant.ALLELE_REFERENCE) {
							votesAllele1++;
						}
					} else if (consensusType == CONSENSUS_ALL) {
						votes++;
						if(call == CalledGenomicVariant.ALLELE_REFERENCE) {
							votesAllele1++;
						}
					}
				}
			}
			//TODO: Check index
			if (2*votesAllele1 < votes) {
				haplotype[i] = CalledGenomicVariant.ALLELE_ALTERNATIVE;
			} else if (2*votesAllele1 > votes){
				haplotype[i] = CalledGenomicVariant.ALLELE_REFERENCE;
			} else {
				haplotype[i] = CalledGenomicVariant.ALLELE_UNDECIDED;
			}
		}*/
		return haplotype;
	}
	/*public static String getHaplotype(Block b, boolean[] cut, double[] fragmentScores) {
		
		int pos=b.getFirstPos();
		int lastPos = b.getLastPos();
		StringBuilder haplotype = new StringBuilder(lastPos-pos+1);
		List <Fragment> fragments = b.getFragments();
		int i=0;
		for(;pos <=lastPos;pos++) {
			while(i<fragments.size() && fragments.get(i).getLastPos()<pos) {
				i++;
			}
			double scoreAllele1=0;
			double scoreAllele2=0;
			int votesAllele1 = 0;
			int votes=0;
			for(int j=i;j<fragments.size();j++) {
				Fragment f2 = fragments.get(j);
				if(f2.getFirstPos()>pos) {
					break;
				}
				char call = f2.getCall(pos);
				if(call != Fragment.NODATACHAR) {
					votes++;
					if((call == Fragment.ALLELE1CHAR && !cut[j])|| (call!=Fragment.ALLELE1CHAR&& cut[j]) ) {
						scoreAllele1+=fragmentScores[j];
						votesAllele1++;
					} else {
						scoreAllele2+=fragmentScores[j];
					}
				}
			}
			if (scoreAllele1 > scoreAllele2) {
				haplotype.append(Fragment.ALLELE1CHAR);
			} else if (scoreAllele2 > scoreAllele1 ) {
				haplotype.append(Fragment.ALLELE2CHAR);
			} else if (2*votesAllele1 < votes) {
				haplotype.append(Fragment.ALLELE2CHAR);
			} else if (2*votesAllele1 > votes){
				haplotype.append(Fragment.ALLELE1CHAR);
			} else {
				haplotype.append(Fragment.NODATACHAR);
			}
		}
		return haplotype.toString();
	}*/
}
