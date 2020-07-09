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
package ngsep.variants;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

public class CalledGenomicVariantImpl implements CalledGenomicVariant {
	private GenomicVariant variant;
	private String sampleId;
	//Positions in the list of called alleles
	private byte [] indexesCalledAlleles = new byte[0];
	private int totalReadDepth=0;
	private short genotypeQuality=0;
	private short totalCopyNumber=DEFAULT_PLOIDY;
	private short [] allelesCopyNumber;
	private VariantCallReport callReport;
	private int [] allCounts;
	private byte [] indexesPhasedAlleles = new byte [0];
	private byte strandBiasScore = INVALID_STRAND_BIAS_SCORE;
	
	public CalledGenomicVariantImpl (String sequenceName, int position, List<String> alleles, byte[] indexesCalledAlleles) {
		variant = new GenomicVariantImpl(sequenceName, position, alleles);
		allelesCopyNumber = new short [alleles.size()];
		setIndexesCalledAlleles(indexesCalledAlleles);
	}
	
	public CalledGenomicVariantImpl (GenomicVariant variant, byte[] indexesCalledAlleles) {
		this.variant = variant;
		String [] alleles = variant.getAlleles();
		allelesCopyNumber = new short [alleles.length];
		setIndexesCalledAlleles(indexesCalledAlleles);
	}
	
	/**
	 * Creates a new call for a biallelic genomic variant with a default encoding of genotypes for a biallelic variant
	 * representing the number of alternative alleles
	 * @param variant to call. It must be either a reference only site or a biallelic variant
	 * @param genotype -1 for undecided, 0 for homozygous reference, 1 for heterozygous, 2 for homozygous alternative
	 */
	public CalledGenomicVariantImpl (GenomicVariant variant, int genotype) {
		this.variant = variant;
		String [] alleles = variant.getAlleles();
		allelesCopyNumber = new short [alleles.length];
		byte [] indexesCalledAlleles = new byte [0];
		if(genotype==0 || genotype == 2) {
			indexesCalledAlleles = new byte[1];
			indexesCalledAlleles[0] = (byte) (genotype/2);
		} else if (genotype == 1) {
			indexesCalledAlleles = new byte [2];
			indexesCalledAlleles[0] = 0;
			indexesCalledAlleles[1] = 1;
		}
		setIndexesCalledAlleles(indexesCalledAlleles);
	}
	/**
	 * Creates a new called genomic variant for the given variant and the genotype given by the calledAlleles
	 * @param variant Genomic variant that was called
	 * @param calledAlleles Alleles called in this variant. PRE: all alleles are different and appear in the list of alleles for the given variant
	 */
	public CalledGenomicVariantImpl (GenomicVariant variant, String [] calledAlleles) {
		this(variant,Arrays.asList(calledAlleles));
	}
	
	/**
	 * Creates a new called genomic variant for the given variant and the genotype given by the calledAlleles
	 * @param variant Genomic variant that was called
	 * @param calledAlleles Alleles called in this variant. PRE: all alleles are different and appear in the list of alleles for the given variant
	 */
	public CalledGenomicVariantImpl (GenomicVariant variant, Collection<String> calledAlleles) {
		this.variant = variant;
		String [] alleles = variant.getAlleles();
		byte [] idxCalledAlleles = new byte[calledAlleles.size()];
		Set<String> callsSet = new HashSet<>();
		Arrays.fill(idxCalledAlleles, (byte)-1);
		Iterator<String> it = calledAlleles.iterator();
		for(int i=0;it.hasNext();i++) {
			String calledAllele = it.next();
			if(callsSet.contains(calledAllele)) throw new RuntimeException("Error creating called genomic variant at "+variant.getSequenceName()+":"+variant.getFirst()+". Duplicated called allele: "+calledAllele);
			callsSet.add(calledAllele);
			for(int j=0;j<alleles.length;j++) {
				if(calledAllele.equals(alleles[j])) {
					idxCalledAlleles[i] = (byte)j;
					break;
				}
			}
			if(idxCalledAlleles[i]==-1) throw new RuntimeException("Error creating called genomic variant at "+variant.getSequenceName()+":"+variant.getFirst()+". Called allele: "+calledAllele+" does not appear in the list of alleles of this variant");
		}
		allelesCopyNumber = new short [alleles.length];
		setIndexesCalledAlleles(idxCalledAlleles);
	}
	
	@Override
	public String getSampleId() {
		return sampleId;
	}

	@Override
	public void setSampleId(String sampleId) {
		this.sampleId = sampleId;
	}

	@Override
	public String[] getAlleles() {
		return variant.getAlleles();
	}
	@Override
	public String getReference() {
		return variant.getReference();
	}

	@Override
	public String getId() {
		return variant.getId();
	}
	
	@Override
	public void setId(String id) {
		variant.setId(id);
	}

	@Override
	public short getVariantQS() {
		return variant.getVariantQS();
	}
	
	@Override
	public void setVariantQS(short variantQS) {
		variant.setVariantQS(variantQS);
	}
	@Override
	public String getSequenceName() {
		return variant.getSequenceName();
	}

	@Override
	public int getFirst() {
		return variant.getFirst();
	}

	@Override
	public int getLast() {
		return variant.getLast();
	}

	@Override
	public int length() {
		return variant.length();
	}

	@Override
	public boolean isPositiveStrand() {
		return variant.isPositiveStrand();
	}

	@Override
	public boolean isNegativeStrand() {
		return variant.isNegativeStrand();
	}
	

	@Override
	public String[] getCalledAlleles() {
		String [] alleles = getAlleles();
		String [] answer = new String [indexesCalledAlleles.length];
		for (int i=0;i<indexesCalledAlleles.length;i++) {
			answer[i] = alleles[indexesCalledAlleles[i]];
		}
		return answer;
	}
	
	@Override
	public short [] getAllelesCopyNumber() {
		return Arrays.copyOf(allelesCopyNumber,allelesCopyNumber.length);
	}
	
	
	/**
	 * Changes the alleles copy number
	 * @param allelesCN Array with the copy number of each allele 
	 */
	@Override
	public void setAllelesCopyNumber(short [] allelesCN) {
		String [] alleles = this.getAlleles();
		if(alleles.length!=allelesCN.length)  throw new IllegalArgumentException("The size of the array must be compatible with the number of alleles");
		boolean [] allelesCalled = new boolean[alleles.length];
		Arrays.fill(allelesCalled, false);
		for(int i=0;i<indexesCalledAlleles.length;i++) allelesCalled[indexesCalledAlleles[i]] = true;
		int totalCN = 0;
		for(int i=0;i<allelesCN.length;i++) {
			if(allelesCN[i]<0)  throw new IllegalArgumentException("Inconsistent copy number detected. Copy number of the allele "+alleles[i]+" can not be the negative value "+allelesCN[i]);
			if(allelesCN[i]>0 && !allelesCalled[i]) throw new IllegalArgumentException("Inconsistent copy number detected. Copy number of the uncalled allele "+alleles[i]+" should not be "+allelesCN[i]);
			totalCN+=allelesCN[i];
		}
		if(totalCN>CalledGenomicVariant.MAX_PLOIDY_SAMPLE)  throw new IllegalArgumentException("Inconsistent copy number detected. Total copy number "+totalCN+" can not be larger than maximum allowed "+CalledGenomicVariant.MAX_PLOIDY_SAMPLE);
		this.totalCopyNumber = (short)totalCN;
		if(isUndecided()) {
			Arrays.fill(allelesCopyNumber, (short)0);
		} else {
			allelesCopyNumber = Arrays.copyOf(allelesCN, allelesCN.length);
		}
		//Remove phasing
		indexesPhasedAlleles = new byte [0];
	}
	
	@Override
	public void updateAllelesCopyNumberFromCounts(short totalCopyNumber) {
		this.totalCopyNumber = totalCopyNumber;
		Arrays.fill(allelesCopyNumber, (short)0);
		if(this.isUndecided()) return;
		if(this.isHomozygousReference()) {
			allelesCopyNumber[0] = totalCopyNumber;
			return;
		}
		int nCalledAlleles = indexesCalledAlleles.length;
		if(totalCopyNumber<=nCalledAlleles) {
			for(int i=0;i<nCalledAlleles;i++) allelesCopyNumber[indexesCalledAlleles[i]]=(short)1;
			return;
		}
		
		
		if(callReport==null || !callReport.countsPresent()) {
			//Distribute equally
			short defAllelePloidy = (short) (totalCopyNumber/nCalledAlleles);
			for(int i=0;i<nCalledAlleles;i++) allelesCopyNumber[indexesCalledAlleles[i]] = defAllelePloidy;
			int diff = totalCopyNumber - defAllelePloidy*indexesCalledAlleles.length;
			allelesCopyNumber[indexesCalledAlleles[0]]+= diff;
			return;
		}
		String [] alleles = getAlleles();
		
		int totalReadCount=0;
		int [] calledAlleleReadCounts = new int [nCalledAlleles];
		for(int i=0;i<nCalledAlleles;i++) {
			String nextAllele = alleles[indexesCalledAlleles[i]];
			calledAlleleReadCounts[i] = callReport.getCount(nextAllele);
			if(calledAlleleReadCounts[i]==0) calledAlleleReadCounts[i]=1; 
			totalReadCount+=calledAlleleReadCounts[i];
		}
		//Normalize by total copy number
		int totalCount=0;
		for(int i=0;i<nCalledAlleles;i++) {
			int j = indexesCalledAlleles[i];
			allelesCopyNumber[j] = (short)Math.max(1, Math.round((double)totalCopyNumber*calledAlleleReadCounts[i]/totalReadCount));
			totalCount+=allelesCopyNumber[j];
		}
		//Assign leftovers to the first allele or remove excess from the last alleles 
		if(totalCount<totalCopyNumber) {
			allelesCopyNumber[indexesCalledAlleles[0]]+=(totalCopyNumber-totalCount);
		} else {
			int excess = totalCount-totalCopyNumber;
			for(int i=nCalledAlleles-1;excess>0 && i>=0;i--) {
				int j = indexesCalledAlleles[i];
				int toRemove = Math.min(excess, allelesCopyNumber[j]-1);
				allelesCopyNumber[j]-=toRemove;
				excess -=toRemove;
			}
		}
		//Remove phasing
		indexesPhasedAlleles = new byte [0];
	}

	@Override
	public byte[] getIndexesCalledAlleles() {
		return Arrays.copyOf(indexesCalledAlleles, indexesCalledAlleles.length);
	}
	
	public void setIndexesCalledAlleles(byte[] indexesCalledAlleles) {
		this.indexesCalledAlleles= Arrays.copyOf(indexesCalledAlleles, indexesCalledAlleles.length);
		short totalCopyNumber = (short)Math.max(indexesCalledAlleles.length, getCopyNumber());
		updateAllelesCopyNumberFromCounts(totalCopyNumber);
	}

	@Override
	public short getCopyNumber () {
		return totalCopyNumber;
	}

	@Override
	public int getTotalReadDepth() {
		return totalReadDepth;
	}
	
	@Override
	public void setTotalReadDepth(int totalReadDepth) {
		this.totalReadDepth = totalReadDepth;
	}

	@Override
	public short getGenotypeQuality() {
		return genotypeQuality;
	}
	@Override
	public void setGenotypeQuality(short genotypeQuality) {
		this.genotypeQuality = genotypeQuality;
	}
	
	@Override
	public void makeUndecided() {
		indexesCalledAlleles = new byte[0];
		genotypeQuality = 0;
		short totalCopyNumber = getCopyNumber();
		updateAllelesCopyNumberFromCounts(totalCopyNumber);
	}

	@Override
	public boolean isUndecided () {
		return indexesCalledAlleles.length==0;
	}
	@Override
	public boolean isHeterozygous() {
		return indexesCalledAlleles.length>1;
	}

	@Override
	public boolean isHomozygous() {
		return indexesCalledAlleles.length == 1;
	}

	@Override
	public boolean isHomozygousReference() {
		return isHomozygous() && indexesCalledAlleles[0]==0;
	}

	@Override
	public VariantCallReport getCallReport() {
		return callReport;
	}
	
	public void setCallReport(VariantCallReport callReport) {
		this.callReport = callReport;
	}

	@Override
	public int[] getAllCounts() {
		return allCounts;
	}

	public void setAllCounts(int[] allCounts) {
		this.allCounts = Arrays.copyOf(allCounts, allCounts.length);
	}

	@Override
	public boolean isPhased() {
		return indexesPhasedAlleles.length>0;
	}

	@Override
	public String[] getPhasedAlleles() {
		String [] alleles = getAlleles();
		String [] answer = new String [indexesPhasedAlleles.length];
		for(int i=0;i<answer.length;i++) {
			answer[i] = alleles[indexesPhasedAlleles[i]];
		}
		return answer;
	}
	@Override
	public byte[] getIndexesPhasedAlleles() {
		return indexesPhasedAlleles;
	}
	/**
	 * PRE: length of phased alleles coincides with total copy number
	 * All input alleles are valid for the variant
	 * @param phasedAlleles
	 */
	public void setIndexesPhasedAlleles(byte [] phasedAlleles) {
		if(phasedAlleles.length!=totalCopyNumber) return;
		this.indexesPhasedAlleles = Arrays.copyOf(phasedAlleles, phasedAlleles.length);
	}

	@Override
	public boolean isCompatible(GenomicVariant variant) {
		if(this == variant) return true;
		return this.variant.isCompatible(variant);
	}

	@Override
	public boolean isBiallelic() {
		return this.variant.isBiallelic();
	}

	@Override
	public boolean isSNV() {
		return this.variant.isSNV();
	}

	@Override
	public byte getType() {
		return variant.getType();
	}

	@Override
	public void setType(byte type) {
		variant.setType(type);
	}

	@Override
	public byte getStrandBiasScore() {
		return strandBiasScore;
	}

	/**
	 * @param strandBiasScore the strandBiasScore to set
	 */
	public void setStrandBiasScore(byte strandBiasScore) {
		this.strandBiasScore = strandBiasScore;
	}
	

	
}
