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
import java.util.Iterator;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.DiversityStatistics;
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
public class QualityStatisticsTetraploidF1Families {
	public static final int NUM_ROWS_COUNTS = 51;
	//Rows are MAF values. Column is number of minor allele count
	private int [][] countsMinorAlleleDosage;
	
	public static void main(String[] args) throws Exception {
		QualityStatisticsTetraploidF1Families instance = new QualityStatisticsTetraploidF1Families();
		String filename = args[0];
		String p1 = null;
		if(args.length>1) p1 = args[1];
		String p2 = null;
		if(args.length>2) p2 = args[2];
		
		instance.processFile(filename,p1,p2);
		instance.printStatistics();
	}

	private void processFile(String filename, String parent1, String parent2) throws IOException {
		countsMinorAlleleDosage = new int [NUM_ROWS_COUNTS][5];
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
				DiversityStatistics divStats = DiversityStatistics.calculateDiversityStatistics(genotypeCalls, false);
				if(!divStats.isPolymorphic()) continue;
				int minorAlleleIndex = divStats.getMafIndex();
				if(minorAlleleIndex<0) continue;
				double maf = divStats.getMaf();
				int row = (int) Math.round(100*maf);
				for(CalledGenomicVariant call:genotypeCalls) {
					short [] acn = call.getAllelesCopyNumber();
					short minorAlleleCount = acn[minorAlleleIndex];
					(countsMinorAlleleDosage[row][minorAlleleCount])++;
				}
			}
		}
	}

	private void printStatistics() {
		System.out.println("MAF\tHomozygous major\tMinor1\tMinor2\tMinor3\tHomozygous minor");
		for(int i=0;i<countsMinorAlleleDosage.length;i++) {
			double maf = ((double)i)/100.0;
			System.out.print(""+maf);
			for(int j=0;j<countsMinorAlleleDosage[i].length;j++) {
				System.out.print("\t"+countsMinorAlleleDosage[i][j]);
			}
			System.out.println();
		}
	}
}
