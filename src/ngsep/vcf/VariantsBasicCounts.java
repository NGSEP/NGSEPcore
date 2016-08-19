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
package ngsep.vcf;


public class VariantsBasicCounts {
	
	public static final int GENOTYPE_STATUS_UNDECIDED = -1;
	public static final int GENOTYPE_STATUS_HOMOREF = 0;
	public static final int GENOTYPE_STATUS_HETEROZYGOUS = 1;
	public static final int GENOTYPE_STATUS_HOMOALT = 2;
	
	public static final int FUNCTIONAL_STATUS_GENERALCODING = 1;
	public static final int FUNCTIONAL_STATUS_SYNONYMOUS = 2;
	public static final int FUNCTIONAL_STATUS_MISSENSE = 3;
	public static final int FUNCTIONAL_STATUS_NONSENSE = 4;
	public static final int FUNCTIONAL_STATUS_FRAMESHIFT = 5;
	
	public static final int POPULATION_STATUS_GENPOP = 1;
	public static final int POPULATION_STATUS_RARE = 2;
	public static final int POPULATION_STATUS_UNIQUE = 3;
	
	
	private int genotyped = 0;
	private int nonReference = 0;
	private int homozygousAlternative = 0;
	private int heterozygous = 0;
	
	private int coding = 0;
	private int synonymous = 0;
	private int missense = 0;
	private int nonsense = 0;
	private int frameshift = 0;
	private int codingHeterozygous = 0;
	private int synonymousHeterozygous = 0;
	private int missenseHeterozygous = 0;
	private int nonsenseHeterozygous = 0;
	private int frameshiftHeterozygous = 0;
	
	private int genotypedPopCounts = 0;
	private int rareAllele = 0;
	private int uniqueAllele = 0;
	
	private int transitions = 0;
	private int homozygousAlternativeTransitions = 0;
	private int heterozygousTransitions = 0;
	private int synonymousTransitions = 0;
	private int nonSynonymousTransitions = 0;
	
	
	public void processGenotypeCall (int genotypingStatus, boolean isTransition, int functionalStatus, int populationStatus) {
		if(genotypingStatus != GENOTYPE_STATUS_UNDECIDED) {
			genotyped++;
			if(genotypingStatus!=GENOTYPE_STATUS_HOMOREF) {
				nonReference++;
				if(isTransition) transitions++;
				if(functionalStatus>0) {
					coding++;
				}
				if(functionalStatus == FUNCTIONAL_STATUS_SYNONYMOUS) {
					synonymous++;
					if(isTransition) synonymousTransitions++;
				}
				if(functionalStatus == FUNCTIONAL_STATUS_MISSENSE) {
					missense++;
					if(isTransition) nonSynonymousTransitions++;
				}
				if(functionalStatus == FUNCTIONAL_STATUS_NONSENSE) {
					nonsense++;
					if(isTransition) nonSynonymousTransitions++;
				}
				if(functionalStatus == FUNCTIONAL_STATUS_FRAMESHIFT) frameshift++;
				if(genotypingStatus==GENOTYPE_STATUS_HETEROZYGOUS) {
					heterozygous++;
					if(isTransition) heterozygousTransitions++;
					if(functionalStatus>0) codingHeterozygous++;
					if(functionalStatus == FUNCTIONAL_STATUS_SYNONYMOUS) synonymousHeterozygous++;
					if(functionalStatus == FUNCTIONAL_STATUS_MISSENSE) missenseHeterozygous++;
					if(functionalStatus == FUNCTIONAL_STATUS_NONSENSE) nonsenseHeterozygous++;
					if(functionalStatus == FUNCTIONAL_STATUS_FRAMESHIFT) frameshiftHeterozygous++;
				} else {
					homozygousAlternative++;
					if(isTransition) homozygousAlternativeTransitions++;
				}
			}
			if(populationStatus > 0) genotypedPopCounts++;
			if(populationStatus == POPULATION_STATUS_RARE || populationStatus == POPULATION_STATUS_UNIQUE) rareAllele++;
			if(populationStatus == POPULATION_STATUS_UNIQUE) uniqueAllele++;
		}
	}


	public int getGenotyped() {
		return genotyped;
	}


	public int getNonReference() {
		return nonReference;
	}


	public int getHomozygousAlternative() {
		return homozygousAlternative;
	}


	public int getHeterozygous() {
		return heterozygous;
	}

	public int getCoding() {
		return coding;
	}

	public int getSynonymous() {
		return synonymous;
	}


	public int getMissense() {
		return missense;
	}


	public int getNonsense() {
		return nonsense;
	}


	public int getFrameshift() {
		return frameshift;
	}

	public int getCodingHeterozygous() {
		return codingHeterozygous;
	}


	public int getSynonymousHeterozygous() {
		return synonymousHeterozygous;
	}


	public int getMissenseHeterozygous() {
		return missenseHeterozygous;
	}


	public int getNonsenseHeterozygous() {
		return nonsenseHeterozygous;
	}

	public int getFrameshiftHeterozygous() {
		return frameshiftHeterozygous;
	}

	public int getGenotypedPopCounts() {
		return genotypedPopCounts;
	}


	public int getRareAllele() {
		return rareAllele;
	}


	public int getUniqueAllele() {
		return uniqueAllele;
	}
	
	public int getTransitions() {
		return transitions;
	}

	public int getHomozygousAlternativeTransitions() {
		return homozygousAlternativeTransitions;
	}

	public int getHeterozygousTransitions() {
		return heterozygousTransitions;
	}
	
	public int getSynonymousTransitions() {
		return synonymousTransitions;
	}

	public int getNonSynonymousTransitions() {
		return nonSynonymousTransitions;
	}

	public double getNonSynonymousToSynonymousRate () {
		return safeDoubleRatio(missense+nonsense,synonymous);
	}
	
	public double getNonSynonymousToSynonymousRateHeterozygous () {
		return safeDoubleRatio(missenseHeterozygous+nonsenseHeterozygous,synonymousHeterozygous);
	}
	
	public double getPCTFrameshift() {
		return 100.0*safeDoubleRatio(frameshift, coding);
	}
	
	public double getPCTFrameshiftHeterozygous() {
		return 100.0*safeDoubleRatio(frameshiftHeterozygous, codingHeterozygous);
	}
	
	public double getTrTvRatio() {
		int transversions = nonReference-transitions;
		return safeDoubleRatio(transitions, transversions); 
	}
	
	public double getTrTvRatioHeterozygous() {
		int transversions = heterozygous-heterozygousTransitions;
		return safeDoubleRatio(heterozygousTransitions, transversions); 
	}
	
	public double getTrTvRatioHomozygousAlternative() {
		int transversions = homozygousAlternative-homozygousAlternativeTransitions;
		return safeDoubleRatio(homozygousAlternativeTransitions, transversions); 
	}
	
	public double getTrTvRatioSynonymous() {
		int transversions = synonymous-synonymousTransitions;
		return safeDoubleRatio(synonymousTransitions, transversions);
	}
	
	public double getTrTvRatioNonSynonymous() {
		int transversions = missense+nonsense-nonSynonymousTransitions;
		return safeDoubleRatio(nonSynonymousTransitions, transversions);
	}
	
	private double safeDoubleRatio(int numerator, int denominator) {
		if(denominator==0) return 0;
		return ((double)numerator)/denominator;
	}
	
}
