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
 * Variants are classified according to estimations of minor allele frequency (MAF):
 * 1. Conserved site (MAF<0.1): This sites should not be present. Genotype calls different than homozygous for the major allele are counted as errors
 * 2. Homozygous contrasting: (MAF>0.1, OH<0.1): Both parents are homozygous but for different alleles. Heterozygous genotypes are counted as errors
 * 3. Parent heterozygous:  (MAF> 0.1, OH>0.1): At least one parent is heterozygous. No errors are counted in this case
 * @author Jorge Duitama
 *
 */
public class QualityStatisticsInbredBiparentalFamilies {
	private Distribution ohDist = new Distribution(0, 1, 0.01);
	private Distribution mafDist = new Distribution(0, 0.51, 0.01);
	private int parentErrors [] = new int [3];
	private int progenyErrors [] = new int [3];
	private int totalVars = 0;
	private int totalDatapoints = 0;
	private int varsPerCategory [] = new int [3];
	private int datapointsPerCategory [] = new int [3];
	private boolean conservative = false;
	
	public static void main(String[] args) throws Exception {
		QualityStatisticsInbredBiparentalFamilies instance = new QualityStatisticsInbredBiparentalFamilies();
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
				//Basic information from the variant
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
				mafDist.processDatapoint(maf);
				//double propHomo = Math.min(n00, n11)/(n00+n11);
				ohDist.processDatapoint(oh);
				if(maf<0.1) {
					if(!conservative || (gP1==gP2 && gP1==maxHomoGenotype)) {
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
				} else if(oh<0.1) {
					//System.out.println("MAF: "+maf+" GP1: "+gP1+" GP2: "+gP2);
					if(!conservative || (gP1!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP2!=CalledGenomicVariant.GENOTYPE_UNDECIDED && Math.abs(gP1-gP2)==2)) {
						varsPerCategory[1]++;
						datapointsPerCategory[1]+=nC;
						
						//Heterozygous genotype calls are probably errors
						if(gP1==CalledGenomicVariant.GENOTYPE_HETERO) parentErrors[1]++;
						if(gP2==CalledGenomicVariant.GENOTYPE_HETERO) parentErrors[1]++;
						progenyErrors[1]+=n01;
						
						//Parents should be homozygous with different alleles
						if(gP1!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP2!=CalledGenomicVariant.GENOTYPE_UNDECIDED && gP1!=CalledGenomicVariant.GENOTYPE_HETERO && gP2!=CalledGenomicVariant.GENOTYPE_HETERO && gP1==gP2) parentErrors[1]++;
					}
				} else {
					varsPerCategory[2]++;
					datapointsPerCategory[2]+=nC;
				}
			}
		}
	}

	private void printStatistics() {
		System.out.print("SUMMARY COUNTS\t"+totalVars+"\t"+totalDatapoints);
		for(int i=0;i<3;i++) {
			System.out.print("\t"+varsPerCategory[i]+"\t"+datapointsPerCategory[i]);
			System.out.print("\t"+parentErrors[i]+"\t"+progenyErrors[i]);
		}
		System.out.println();
		System.out.println("MAF distribution");
		mafDist.printDistribution(System.out);
		System.out.println();
		System.out.println("Observed heterozygosity distribution");
		ohDist.printDistribution(System.out);
	}

	private int getGenotype(CalledGenomicVariant call) {
		if(call.isUndecided()) return CalledGenomicVariant.GENOTYPE_UNDECIDED;
		if(call.isHomozygousReference()) return CalledGenomicVariant.GENOTYPE_HOMOREF;
		if(call.isHeterozygous()) return CalledGenomicVariant.GENOTYPE_HETERO;
		return CalledGenomicVariant.GENOTYPE_HOMOALT;
	}
}
