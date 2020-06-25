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
import java.util.LinkedList;


/**
 * Information for a SNV called for an individual
 * @author Jorge Duitama
 *
 */
public class CalledSNV implements CalledGenomicVariant {
	
	private SNV snv;
	private String sampleId;
	private byte genotype = GENOTYPE_UNDECIDED; //-1 for undecided, 0 for homozygous reference, 1 for heterozygous, 2 for homozygous variant
	private short genotypeQuality = 0;
	private int totalReadDepth = 0; //Total read depth of this variant in the sample. It can be more than the sum of the allele counts
	private int countA = 0;
	private int countC = 0;
	private int countG = 0;
	private int countT = 0;
	private float homoRefLogCond=0;
	private float homoAltLogCond=0;
	private float refAltLogCond=0;
	private float altRefLogCond=0;
 
	//The copy number of the called variant. 
	private short totalCopyNumber = 0;
	
	//The copy number of the reference allele. 
	private short refCopyNumber = 0;
	
	//Phasing as a bit array
	private byte phasing = -1;
	
	//Strand bias score
	private byte strandBiasScore =INVALID_STRAND_BIAS_SCORE;
	
	
	/**
	 * Creates a new CalledSNV with the given information
	 * @param sequenceName Name of the sequence where the SNV is located
	 * @param absolutePosition Position in the sequence where the SNp is located
	 * @param refBase Base present in the reference sequence
	 * @param altBase Base representing the minor allele 
	 * @param genotype 0 for homozygous reference, 1 for heterozygous, 2 for homozygous alternative, -1 for undecided
	 */
	public CalledSNV(String sequenceName, int position, char refBase, char altBase, byte genotype) {
		snv = new SNV(sequenceName, position, refBase, altBase);
		refCopyNumber = DEFAULT_PLOIDY;
		this.setGenotype(genotype);
	}
	/**
	 * Calls the given genotype for the given snv
	 * @param snv Variant to call
	 * @param genotype 0 for homozygous reference, 1 for heterozygous, 2 for homozygous alternative, -1 for undecided
	 */
	public CalledSNV(SNV snv, byte genotype) {
		this.snv = snv;
		refCopyNumber = DEFAULT_PLOIDY;
		this.setGenotype(genotype);
	}

	@Override
	public String getSampleId() {
		return sampleId;
	}

	@Override
	public void setSampleId(String sampleId) {
		this.sampleId = sampleId;
	}

	public byte getGenotype() {
		return genotype;
	}
	
	/**
	 * Changes the called genotype. If the genotype is changed from homozygous to heterozygous and
	 * the total copy number is larger than 2, the two alleles are equally distributed.
	 * @param genotype New genotype. 0 for homozygous reference, 1 for heterozygous, 2 for homozygous alternative
	 */
	public void setGenotype(byte genotype) {
		if(this.genotype == genotype) return;
		if(genotype <GENOTYPE_UNDECIDED || genotype>GENOTYPE_HOMOALT) throw new IllegalArgumentException("Invalid genotype for a called SNV: "+genotype);
		this.genotype = genotype;
		short totalCopyNumber = getCopyNumber();
		updateAllelesCopyNumberFromCounts(totalCopyNumber);
	}
	
	/**
	 * Changes the alleles copy number
	 * @param allelesCN Array with the copy number of each allele 
	 */
	public void setAllelesCopyNumber(short [] allelesCN) {
		if(allelesCN.length!=2) throw new IllegalArgumentException("For Called SNVs the length of the array should be 2");
		if(allelesCN[0]<0) throw new IllegalArgumentException("Inconsistent copy number detected. Copy number of the reference allele can not be the negative value "+allelesCN[0]);
		if(allelesCN[1]<0) throw new IllegalArgumentException("Inconsistent copy number detected. Copy number of the alternative allele can not be the negative value "+allelesCN[1]);
		int totalCN = allelesCN[0] + allelesCN[1];
		if(totalCN>CalledGenomicVariant.MAX_PLOIDY_SAMPLE)  throw new IllegalArgumentException("Inconsistent copy number detected. Total copy number "+totalCN+" can not be larger than maximum allowed "+CalledGenomicVariant.MAX_PLOIDY_SAMPLE);
		if(this.isHomozygousReference() && allelesCN[1]!=0) throw new IllegalArgumentException("Inconsistent copy number detected. Copy number of the alternative allele can not be the non-zero value "+allelesCN[1]+" for a homozygous reference call");
		else if(this.genotype == GENOTYPE_HOMOALT && allelesCN[0]!=0) throw new IllegalArgumentException("Inconsistent copy number detected. Copy number of the reference allele can not be the non-zero value "+allelesCN[0]+" for a homozygous non-reference call");
		this.totalCopyNumber = (short)(totalCN);
		if(this.isUndecided()) {
			this.refCopyNumber = (short)0;
		} else {
			this.refCopyNumber = allelesCN[0];
		}
		//Unphase variant because the copy number changed
		phasing = -1;
	}
	
	@Override
	public void updateAllelesCopyNumberFromCounts(short totalCopyNumber) {
		this.totalCopyNumber = totalCopyNumber;
		refCopyNumber = (short)0;
		if (this.isUndecided() ) return;
		if (this.isHomozygousReference()) {
			refCopyNumber = totalCopyNumber;
			return;
		} else if (this.isHomozygous()) {
			return;
		} else if (totalCopyNumber<=2) {
			this.totalCopyNumber = (short)2;
			refCopyNumber = (short)1;
			return;
		}
		double countRef = getCountReference();
		double sumCounts = getCountReference() + getCountAlternative();
		double refProp = 0.5;
		if(sumCounts > 0) refProp = countRef/sumCounts;
		if(refProp>1) refProp = 1;
		this.refCopyNumber = (short) Math.round(refProp*totalCopyNumber);
		if(this.refCopyNumber==0) this.refCopyNumber = 1;
		else if (this.refCopyNumber>=totalCopyNumber) this.refCopyNumber = (short) (totalCopyNumber-1);
		//Unphase variant because the copy number changed
		phasing = -1;
	}
	
	@Override
	public short [] getAllelesCopyNumber() {
		short [] answer = new short [2];
		Arrays.fill(answer, (short)0);
		if(this.isUndecided()) return answer;
		answer[0] = refCopyNumber;
		answer[1] = (short) (totalCopyNumber - refCopyNumber);
		return answer;
	}
	
	@Override
	public short getCopyNumber() {
		return totalCopyNumber;
	}
	
	/**
	 * @return int Number of times that the reference allele was observed
	 */
	public int getCountReference() {
		char base = snv.getReferenceBase();
		switch(base) {
			case 'A': return countA;
			case 'C': return countC;
			case 'G': return countG;
			case 'T': return countT;
		}
		return 0;
	}

	/**
	 * Sets the number of times that the reference allele was observed
	 * @param countReference New number of times that the reference allele was observed
	 */
	public void setCountReference(int countReference) {
		char base = snv.getReferenceBase();
		switch(base) {
			case 'A': countA = countReference;break;
			case 'C': countC = countReference;break;
			case 'G': countG = countReference;break;
			case 'T': countT = countReference;break;
		}
	}

	/**
	 * @return int Number of times that the alternative allele was observed
	 */
	public int getCountAlternative() {
		char base = snv.getAlternativeBase();
		switch(base) {
			case 'A': return countA;
			case 'C': return countC;
			case 'G': return countG;
			case 'T': return countT;
		}
		return 0;
	}

	/**
	 * Sets the number of times that the alternative allele was observed
	 * @param countAlternative New number of times that the alternative allele was observed
	 */
	public void setCountAlternative(int countAlternative) {
		char base = snv.getAlternativeBase();
		switch(base) {
			case 'A': countA = countAlternative;break;
			case 'C': countC = countAlternative;break;
			case 'G': countG = countAlternative;break;
			case 'T': countT = countAlternative;break;
		}
	}

	@Override
	public short getGenotypeQuality() {
		return genotypeQuality;
	}
	@Override
	public void setGenotypeQuality(short genotypeQuality) {
		this.genotypeQuality = genotypeQuality;
	}


	/**
	 * Changes the number of As, Cs, Gs and Ts observed in the sample
	 * @param allCounts Array of length 4 with the number of times that each allele appears. 
	 * The order must be the count for A, the count for C, the count for G, and the count for T 
	 */
	public void setAllBaseCounts(int [] allCounts) {
		this.countA = allCounts[0];
		this.countC = allCounts[1];
		this.countG = allCounts[2];
		this.countT = allCounts[3];
	}
	/**
	 * Changes the log conditional probabilities based on the full matrix of probabilities for the 16 possible sorted genotypes
	 * Indexes are 0 for A, 1 for C, 2 for G, and 3 for T
	 * @param allProbabilities Matrix of conditional probabilities
	 */
	public void setAllGenotypeLogConditionals(double [][] allLogConditionals) {
		int refIdx = snv.getRefBaseDNAIndex();
		int altIdx = snv.getAltBaseDNAIndex();
		this.homoRefLogCond = (float)allLogConditionals[refIdx][refIdx];
		this.homoAltLogCond = (float)allLogConditionals[altIdx][altIdx];
		this.refAltLogCond = (float)allLogConditionals[refIdx][altIdx];
		this.altRefLogCond = (float)allLogConditionals[altIdx][refIdx];
	}
	/**
	 * Changes the log conditional probabilities based on a 2x2 matrix with the probabilities of the four sorted genotypes
	 * product of this reference and alternative allele
	 * @param logConditionals 2x2 matrix with the posterior probabilities
	 */
	public void setRefAltGenotypeLogConditionals(double [][] logConditionals) {
		this.homoRefLogCond = (float)logConditionals[0][0];
		this.homoAltLogCond = (float)logConditionals[1][1];
		this.refAltLogCond = (float)logConditionals[0][1];
		this.altRefLogCond = (float)logConditionals[1][0];
	}
	@Override
	public void makeUndecided() {
		genotype = GENOTYPE_UNDECIDED;
		genotypeQuality = 0;
		refCopyNumber = (short)0;
		phasing = -1;
	}
	@Override
	public boolean isUndecided () {
		return genotype==GENOTYPE_UNDECIDED;
	}
	/**
	 * Tells if the genotype is homozygous
	 * @return boolean true if only one allele is called, false otherwise
	 */
	public boolean isHomozygous() {
		return genotype==GENOTYPE_HOMOREF || genotype==GENOTYPE_HOMOALT;
	}
	/**
	 * Tells if the genotype is homozygous for the reference allele
	 * @return boolean true if only the reference allele is called, false otherwise
	 */
	public boolean isHomozygousReference() {
		return genotype==GENOTYPE_HOMOREF;
	}
	/**
	 * Tells if the genotype is heterozygous
	 * @return boolean true if both alleles are called, false otherwise
	 */
	public boolean isHeterozygous() {
		return genotype==GENOTYPE_HETERO;
	}
	@Override
	public String[] getAlleles() {
		return snv.getAlleles();
	}
	@Override
	public String getReference() {
		return snv.getReference();
	}
	@Override
	public String getId() {
		return snv.getId();
	}
	
	@Override
	public void setId(String id) {
		snv.setId(id);
		
	}
	@Override
	public short getVariantQS() {
		return snv.getVariantQS();
	}
	@Override
	public void setVariantQS(short variantQS) {
		snv.setVariantQS(variantQS);
	}
	@Override
	public String getSequenceName() {
		return snv.getSequenceName();
	}

	@Override
	public int getFirst() {
		return snv.getFirst();
	}
	@Override
	public int getLast() {
		return snv.getLast();
	}
	@Override
	public int length() {
		return snv.length();
	}
	@Override
	public boolean isPositiveStrand() {
		return snv.isPositiveStrand();
	}
	@Override
	public boolean isNegativeStrand() {
		return snv.isNegativeStrand();
	}
	
	public int getPosition() {
		return snv.getPosition();
	}
	public char getReferenceBase() {
		return snv.getReferenceBase();
	}
	public char getAlternativeBase() {
		return snv.getAlternativeBase();
	}
	
	@Override
	public String[] getCalledAlleles() {
		if(this.isUndecided()) return new String[0];
		String [] answer;
		String [] alleles = snv.getAlleles();
		if(this.isHomozygous()) { 
			answer = new String [1];
			answer[0] = genotype==0?alleles[0]:alleles[1];
		} else {
			answer = new String [2];
			answer[0] = alleles[0];
			answer[1] = alleles[1];
		}
		return answer;
	}
	@Override
	public byte[] getIndexesCalledAlleles() {
		if(this.isUndecided()) return new byte[0];
		byte [] answer;
		if(this.isHomozygous()) { 
			answer = new byte [1];
			answer[0] = genotype==0?(byte)0:(byte)1;
		} else {
			answer = new byte [2];
			answer[0] = (byte)0;
			answer[1] = (byte)1;
		}
		return answer;
	}
	
	@Override
	public int getTotalReadDepth() {
		return totalReadDepth;
	}
	
	@Override
	public void setTotalReadDepth(int depth) {
		this.totalReadDepth = depth;
	}
	
	@Override
	public VariantCallReport getCallReport() {
		String [] alleles = snv.getAlleles();
		int [] counts = null;
		if(countA+countC+countG+countT > 0) {
			counts = new int[2];
			counts[0] = getCountReference();
			counts[1] = getCountAlternative();
		}
		double [][] logConditionals = null;
		if(homoRefLogCond+refAltLogCond+altRefLogCond+homoAltLogCond!=0) {
			logConditionals = new double[2][2];
			logConditionals[0][0] = homoRefLogCond;
			logConditionals[0][1] = refAltLogCond;
			logConditionals[1][0] = altRefLogCond;
			logConditionals[1][1] = homoAltLogCond;
		}
		return new VariantCallReport(alleles, counts, logConditionals);
	}
	@Override
	public boolean isPhased() {
		return phasing >=0;
	}
	@Override
	public String[] getPhasedAlleles() {
		if(!isPhased()) return null;
		
		String [] alleles = snv.getAlleles();
		LinkedList<String> answer = new LinkedList<String>();
		byte p = this.phasing;
		for(int i=0;i<totalCopyNumber;i++) {
			answer.add(0,alleles[p%2]);
			p/=2;
		}
		return answer.toArray(new String [0]);
	}
	@Override
	public byte[] getIndexesPhasedAlleles() {
		if(!isPhased()) return null;
		byte[] answer = new byte [totalCopyNumber];
		byte p = this.phasing;
		for(int i=0;i<totalCopyNumber;i++) {
			answer[totalCopyNumber - i - 1] = (byte)(p%2);
			p/=2;
		}
		return answer;
	}
	public void setPhasingCN2 (boolean alternativeFirst) {
		if(totalCopyNumber!=2) {
			throw new RuntimeException("This method should not be called with a copy number different than 2");
		}
		if(genotype == GENOTYPE_UNDECIDED) return;
		else if (genotype == GENOTYPE_HOMOREF) phasing = 0;
		else if (genotype == GENOTYPE_HOMOALT) phasing = 3;
		else if (alternativeFirst) phasing = 2;
		else phasing = 1;
	}
	@Override
	public int [] getAllCounts () {
		int [] counts = new int[4];
		counts[0] = countA;
		counts[1] = countC;
		counts[2] = countG;
		counts[3] = countT;
		return counts;
	}
	@Override
	public boolean isCompatible(GenomicVariant variant) {
		if(this == variant) return true;
		return this.snv.isCompatible(variant);
	}
	@Override
	public boolean isBiallelic() {
		return this.snv.isBiallelic();
	}
	
	@Override
	public boolean isSNV() {
		return this.snv.isSNV();
	}
	
	public byte getRefBaseDNAIndex() {
		return snv.getRefBaseDNAIndex();
	}
	public byte getAltBaseDNAIndex() {
		return snv.getAltBaseDNAIndex();
	}
	@Override
	public byte getType() {
		return snv.getType();
	}
	@Override
	public void setType(byte type) {
		snv.setType(type);
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
