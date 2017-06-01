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

import java.util.Arrays;

import ngsep.variants.CalledGenomicVariant;


public class DGSSIHAlgorithm implements SIHAlgorithm {
	private boolean [] cut;
	
	@Override
	public void buildHaplotype(HaplotypeBlock block) {
		byte [] haplotype = new byte[0]; 
		cut = new boolean [block.getNumFragments()];
		System.err.println("Running DGS for: "+cut.length+" fragments");
		initCut(block);
		for (int i=0;i<1000;i++) {
			System.err.println("DGS iteration: "+i);
			byte [] newhap=CutHaplotypeTranslator.getHaplotype(block, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);
			//TODO: check if this works
			if(newhap.equals(haplotype)) {
				break;
			}
			haplotype = newhap;
			updateCut( block, haplotype);
		}
		block.setHaplotype(haplotype);

	}
	private void initCut(HaplotypeBlock b) {
		byte [] hap = new byte [b.getNumVariants()];
		 
		for(int i=0;i<hap.length;i++) {
			hap[i] = CalledGenomicVariant.ALLELE_UNDECIDED;
		}
		int n = b.getNumFragments();
		boolean [] assigned = new boolean [n];
		Arrays.fill(assigned, false);
		//Get the fragment with larger number of calls
		int maxC=0;
		int maxI=0;
		for(int i=0;i<n;i++) {
			int calls = b.getFragmentCalls(i);
			if(calls>maxC) {
				maxC=calls;
				maxI=i;
			}
		}
		assigned[maxI] = true;
		cut[maxI] = false;
		updateHaplotype(hap, b, maxI, false);
		//Assign the other fragments
		for(int h=0;h<n-1;h++) {
			int maxAbsScore =0;
			int maxScore =0;
			//Find the best next fragment
			for(int i=0;i<n;i++) {
				if(!assigned[i]) {
					int score = b.getHamming2(hap,i);
					int absScore = Math.abs(score);
					if(absScore > maxAbsScore) {
						maxAbsScore = absScore;
						maxI = i;
						maxScore = score;
					}
				}
			}
			assigned[maxI] = true;
			cut[maxI] = (maxScore==maxAbsScore);
			updateHaplotype(hap, b, maxI, cut[maxI]);
		}
	}
	private void updateHaplotype(byte [] hap, HaplotypeBlock b, int row, boolean reverse) {
		int firstJ = b.getFirstColumn(row);
		int lastJ = b.getLastColumn(row);
		for(int j=firstJ;j<=lastJ;j++) {
			byte allele = b.getAllele(row, j);
			if(hap[j]==CalledGenomicVariant.ALLELE_UNDECIDED && allele!=CalledGenomicVariant.ALLELE_UNDECIDED) {
				//TODO: Fix for multiallelic
				if(!reverse) {
					hap[j] = allele;
				} else if (allele == CalledGenomicVariant.ALLELE_REFERENCE) {
					hap[j] = CalledGenomicVariant.ALLELE_ALTERNATIVE;
				} else {
					hap[j] = CalledGenomicVariant.ALLELE_REFERENCE;
				}
				
			}
		}
		
	}
	private void updateCut(HaplotypeBlock b, byte [] haplotype) {
		int n = b.getNumFragments();
		for(int i=0;i<n;i++) {
			int score = b.getHamming2(haplotype,i);
			if(score != 0) {
				cut[i] = score > 0;
			}
			
		}
	}

}
