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
package ngsep.benchmark;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import ngsep.math.Distribution;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;


/**
 * Distributes variants in an F1 population in 4 types corresponding to 4 possible genotype configuration in the parents and count genotype errors based on the category assignment. 
 * Variants are classified according to estimations of minor allele frequency (MAF) and observed heterozygosity (OH):
 * 1. Conserved site (OH<0.2): This sites should not be present. Genotype calls different than homozygous for the major allele are counted as errors
 * 2. Homozygous/Heterozygous (0.2<=OH<=0.8, MAF<0.37): One parent is heterozygous and the other is homozygous. Homozygous genotype calls to the minor allele are counted as errors
 * 3. Double heterozygous (0.2<=OH<=0.8, MAF>=0.37): Both parents are heterozygous. No errors can be detected in the population for genotype calls in this category
 * 4. Homozygous contrasting: (OH>0.8): Both parents are homozygous but for different alleles. Genotypes different than heterozygous are counted as errors
 * @author Jorge Duitama
 *
 */
public class QualityStatisticsOutbredF1Families {
	private Distribution ohDist = new Distribution(0, 1, 0.01);
	private Distribution mafOH0 = new Distribution(0, 0.51, 0.01);
	private Distribution mafOH50 = new Distribution(0, 0.51, 0.01);
	private Distribution mafOH100 = new Distribution(0, 0.51, 0.01);
	private Distribution mafOH0NonMendelian = new Distribution(0, 0.51, 0.01);
	private int parentErrors [] = new int [4];
	private int progenyErrors [] = new int [4];
	private int totalVars = 0;
	private int totalDatapoints = 0;
	private int varsPerCategory [] = new int [4];
	private int datapointsPerCategory [] = new int [4];
	private boolean conservative = false;
	
	public static void main(String[] args) throws Exception {
		QualityStatisticsOutbredF1Families instance = new QualityStatisticsOutbredF1Families();
		String filename = args[0];
		String p1 = null;
		if(args.length>1) p1 = args[1];
		String p2 = null;
		if(args.length>2) p2 = args[2];
		instance.conservative = args.length>3;
		instance.processFile(filename,p1,p2);
		instance.printStatistics();
	}

	private void processFile(String filename, String parent1, String parent2) throws IOException {
		Arrays.fill(parentErrors, 0);
		Arrays.fill(progenyErrors, 0);
		 
		try (VCFFileReader reader = new VCFFileReader(filename)){ 
			VCFFileHeader header = reader.getHeader();
			int idx1 = -1;
			if(parent1!=null) {
				idx1 = header.getIndexSampleId(parent1);
				if(idx1 == -1) throw new IOException("Parent "+parent1+" not found in vcf file");
			}
			int idx2 = -1;
			if(parent2!=null) {
				idx2 = header.getIndexSampleId(parent2);
				if(idx2 == -1) throw new IOException("Parent "+parent2+" not found in vcf file");
			}
			Iterator<VCFRecord> it = reader.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				GenomicVariant variant = record.getVariant();
				if(!variant.isBiallelic()) continue;
				List<CalledGenomicVariant> genotypeCalls = record.getCalls();
				double n00 = 0;
				double n01 = 0;
				double n11 = 0;
				double n0 = 0;
				double n1 = 0;
				double n = 0;
				double nC = 0;
				double nC0 = 0;
				double nC1 = 0;
				for(int i=0;i<genotypeCalls.size();i++) {
					boolean isParent = (i==idx1 || i==idx2);
					CalledGenomicVariant genotypeCall = genotypeCalls.get(i);
					if(genotypeCall.isUndecided()) continue;
					nC++;
					if(!isParent)n++;
					if(genotypeCall.isHomozygousReference()) {
						nC0+=2;
						if(!isParent) {
							n00++;
							n0+=2;
						}
					} else if (genotypeCall.isHeterozygous()) {
						nC0++;
						nC1++;
						if(!isParent) {
							n01++;
							n0++;
							n1++;
						}
					} else {
						nC1+=2;
						if(!isParent) {
							n11++;
							n1+=2;
						}
						
					}
				}
				if(n<2) continue;
				if(nC0==0) continue;
				if(nC1==0) continue;
				
				
				int gP1 = CalledGenomicVariant.GENOTYPE_UNDECIDED;
				if(idx1>=0) gP1 = getGenotype(genotypeCalls.get(idx1));
				int gP2 = CalledGenomicVariant.GENOTYPE_UNDECIDED;
				if(idx2>=0) gP2 = getGenotype(genotypeCalls.get(idx2));
				totalVars++;
				totalDatapoints+=nC;
				double oh = n01/n;
				double minAF = Math.min(n0, n1);
				byte maxHomoGenotype = CalledGenomicVariant.GENOTYPE_HOMOREF;
				if(n00<n11) maxHomoGenotype = CalledGenomicVariant.GENOTYPE_HOMOALT;
				double maf = minAF/(n0+n1);
				//double propHomo = Math.min(n00, n11)/(n00+n11);
				ohDist.processDatapoint(oh);
				if(oh<0.2) {
					mafOH0.processDatapoint(maf);
					if(gP1!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP2!=CalledGenomicVariant.GENOTYPE_UNDECIDED && (gP1!=gP2 || gP1==CalledGenomicVariant.GENOTYPE_HETERO)) {
						mafOH0NonMendelian.processDatapoint(maf);
					}
					if(!conservative || (gP1==gP2 && gP1!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP1==maxHomoGenotype)) {
						//Monomorphic site	
						varsPerCategory[0]++;
						datapointsPerCategory[0]+=nC;
						
						//Genotype calls different than homozygous for major allele are errors
						if(gP1!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP1!=maxHomoGenotype) parentErrors[0]++;
						if(gP2!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP2!=maxHomoGenotype) parentErrors[0]++;
						
						progenyErrors[0]+=n01;
						progenyErrors[0]+=Math.min(n00, n11);
						//System.out.println("Monomorphic at "+variant.getSequenceName()+":"+variant.getFirst()+" Parent 1: "+gP1+" Parent 2: "+gP2+" MAF: "+maf+" n0: "+n0+" n1: "+n1+" nC0: "+nC0+" nC1: "+nC1);
					}
					
				} else if (oh<=0.8) {
					mafOH50.processDatapoint(maf);
					if(maf<0.37) {
						if(!conservative || (gP1!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP2!=CalledGenomicVariant.GENOTYPE_UNDECIDED && Math.abs(gP1-gP2)==1)) {
							//One parent should be homozygous and the other heterozygous
							varsPerCategory[1]++;
							datapointsPerCategory[1]+=nC;
							if(gP1!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP2!=CalledGenomicVariant.GENOTYPE_UNDECIDED && Math.abs(gP1-gP2)!=1) parentErrors[1]++;
							progenyErrors[1]+=Math.min(n00, n11);
						}
						
					} else {
						if(!conservative || (gP1==CalledGenomicVariant.GENOTYPE_HETERO && gP2==CalledGenomicVariant.GENOTYPE_HETERO)) {
							//Both parents should be heterozygous. Only parent errors can be detected because the three genotypes are possible
							varsPerCategory[2]++;
							datapointsPerCategory[2]+=nC;
							if(gP1!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP1!=CalledGenomicVariant.GENOTYPE_HETERO) parentErrors[2]++;
							if(gP2!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP2!=CalledGenomicVariant.GENOTYPE_HETERO) parentErrors[2]++;
						}
					}
				} else {
					mafOH100.processDatapoint(maf);
					if(!conservative || (gP1!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP2!=CalledGenomicVariant.GENOTYPE_UNDECIDED && Math.abs(gP1-gP2)==2)) {
						//Parents should be homozygous with different alleles
						varsPerCategory[3]++;
						datapointsPerCategory[3]+=nC;
						if(gP1==CalledGenomicVariant.GENOTYPE_HETERO) parentErrors[3]++;
						if(gP2==CalledGenomicVariant.GENOTYPE_HETERO) parentErrors[3]++;
						if(gP1!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP2!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP1!=CalledGenomicVariant.GENOTYPE_HETERO && gP2!=CalledGenomicVariant.GENOTYPE_HETERO && gP1 ==gP2) parentErrors[3]++;
						progenyErrors[3]+=n00;
						progenyErrors[3]+=n11;
					}
				}
			}
		}
	}

	private void printStatistics() {
		System.out.print("SUMMARY COUNTS\t"+totalVars+"\t"+totalDatapoints);
		for(int i=0;i<4;i++) {
			System.out.print("\t"+varsPerCategory[i]+"\t"+datapointsPerCategory[i]);
			System.out.print("\t"+parentErrors[i]+"\t"+progenyErrors[i]);
		}
		System.out.println();
		System.out.println("Observed heterozygosity distribution");
		ohDist.printDistribution(System.out);
		System.out.println();
		System.out.println("MAF distributions");
		double [] mafOH0D = mafOH0.getDistribution();
		double [] mafOH50D = mafOH50.getDistribution();
		double [] mafOH100D = mafOH100.getDistribution();
		double [] mafOH0NM = mafOH0NonMendelian.getDistribution();
		for(int i=0;i<mafOH0D.length;i++) {
			double minMAF = 0.01*i;
			System.out.println(""+minMAF+"\t"+mafOH0D[i]+"\t"+mafOH50D[i]+"\t"+mafOH100D[i]+"\t"+mafOH0NM[i]);
		}
	}

	private int getGenotype(CalledGenomicVariant call) {
		if(call.isUndecided()) return CalledGenomicVariant.GENOTYPE_UNDECIDED;
		if(call.isHomozygousReference()) return CalledGenomicVariant.GENOTYPE_HOMOREF;
		if(call.isHeterozygous()) return CalledGenomicVariant.GENOTYPE_HETERO;
		return CalledGenomicVariant.GENOTYPE_HOMOALT;
	}
}
