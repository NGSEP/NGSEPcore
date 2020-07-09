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

import java.util.HashMap;
import java.util.Map;

import ngsep.transcriptome.VariantFunctionalAnnotation;
import ngsep.transcriptome.VariantFunctionalAnnotationType;

public class VariantsBasicCounts {
	
	public static final int GENOTYPE_STATUS_UNDECIDED = -1;
	public static final int GENOTYPE_STATUS_HOMOREF = 0;
	public static final int GENOTYPE_STATUS_HETEROZYGOUS = 1;
	public static final int GENOTYPE_STATUS_HOMOALT = 2;
	
	public static final int POPULATION_STATUS_GENPOP = 1;
	public static final int POPULATION_STATUS_RARE = 2;
	public static final int POPULATION_STATUS_UNIQUE = 3;
	
	
	private int genotyped = 0;
	private int nonReference = 0;
	private int homozygousAlternative = 0;
	private int heterozygous = 0;
	private int transitions = 0;
	private int homozygousAlternativeTransitions = 0;
	private int heterozygousTransitions = 0;
	
	private Map<String, Integer> totalCountsPerAnnotation = new HashMap<>();
	private Map<String, Integer> hetCountsPerAnnotation = new HashMap<>();
	private Map<String, Integer> transitionCountsPerAnnotation = new HashMap<>();
	
	private int genotypedPopCounts = 0;
	private int rareAllele = 0;
	private int uniqueAllele = 0;
	
	
	public void processGenotypeCall (int genotypingStatus, boolean isTransition, VariantFunctionalAnnotation annotation, int populationStatus) {
		if(genotypingStatus != GENOTYPE_STATUS_UNDECIDED) {
			genotyped++;
			if(genotypingStatus!=GENOTYPE_STATUS_HOMOREF) {
				nonReference++;
				String key = null;
				if(annotation!=null) key = annotation.getTypeName();
				add1(totalCountsPerAnnotation,key);
				if(isTransition) {
					transitions++;
					add1(transitionCountsPerAnnotation,key);
				}
				
				if(genotypingStatus==GENOTYPE_STATUS_HETEROZYGOUS) {
					heterozygous++;
					if(isTransition) heterozygousTransitions++;
					add1(hetCountsPerAnnotation,key);
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

	private static void add1(Map<String, Integer> countsMap, String key) {
		if(key==null) return;
		Integer count = countsMap.get(key);
		if(count == null) count = 0;
		count++;
		countsMap.put(key, count);
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

	public int getTotalCount(String annotation) {
		return getCount(totalCountsPerAnnotation, annotation);
	}
	
	public int getHeterozygousCount(String annotation) {
		return getCount(hetCountsPerAnnotation, annotation);
	}
	
	public int getTransitionCount(String annotation) {
		return getCount(transitionCountsPerAnnotation, annotation);
	}
	
	private static int getCount (Map<String, Integer> countsMap, String key) {
		if(key==null) return 0;
		Integer count = countsMap.get(key);
		if(count==null) return 0;
		return count;
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
	public int getCodingTotalCount () {
		return getCodingCount(totalCountsPerAnnotation);
	}
	public int getCodingHeterozygousCount () {
		return getCodingCount(hetCountsPerAnnotation);
	}
	public int getCodingTransitionCount () {
		return getCodingCount(transitionCountsPerAnnotation);
	}
	public int getNonSynonymousTotalCount () {
		return getNonSynonymousCount(totalCountsPerAnnotation);
	}
	public int getNonSynonymousHeterozygousCount () {
		return getNonSynonymousCount(hetCountsPerAnnotation);
	}
	public int getNonSynonymousTransitionCount () {
		return getNonSynonymousCount(transitionCountsPerAnnotation);
	}
	
	public int getTotalSpliceRegions() {
		int count = getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_SPLICE_DONOR);
		count += getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_SPLICE_ACCEPTOR);
		count += getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_SPLICE_REGION);
		count += getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_EXONIC_SPLICE_REGION);
		return count;
	}
	
	public int getHeterozygousSpliceRegions() {
		int count = getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_SPLICE_DONOR);
		count += getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_SPLICE_ACCEPTOR);
		count += getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_SPLICE_REGION);
		count += getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_EXONIC_SPLICE_REGION);
		return count;
	}
	private static int getCodingCount(Map<String, Integer> countsMap) {
		int answer = 0;
		for(Map.Entry<String, Integer> entry:countsMap.entrySet()) {
			if(VariantFunctionalAnnotationType.isTypeCoding(entry.getKey())) {
				answer+= entry.getValue();
			}
		}
		return answer;
	}
	public static int getNonSynonymousCount(Map<String, Integer> countsMap) {
		int answer = 0;
		for(Map.Entry<String, Integer> entry:countsMap.entrySet()) {
			if(VariantFunctionalAnnotationType.isTypeNonSynonymous(entry.getKey())) {
				answer+= entry.getValue();
			}
		}
		return answer;
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
		int synTransitions = getTransitionCount(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS);
		int synTransversions = getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS)-synTransitions;
		return safeDoubleRatio(synTransitions, synTransversions);
	}
	
	public double getTrTvRatioNonSynonymous() {
		int nonSynonymousTransitions = getNonSynonymousTransitionCount();
		int transversions = getNonSynonymousTotalCount()-nonSynonymousTransitions;
		return safeDoubleRatio(nonSynonymousTransitions, transversions);
	}
	
	private double safeDoubleRatio(int numerator, int denominator) {
		if(denominator==0) return 0;
		return ((double)numerator)/denominator;
	}
	
}
