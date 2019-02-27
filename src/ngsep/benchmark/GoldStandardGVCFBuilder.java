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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

/**
 * Builds a gold standard genomic VCF from a reference genome, a set of confidence regions,
 * a VCF file with gold standard variants and a file with STRs
 * @author Jorge Duitama
 *
 */
public class GoldStandardGVCFBuilder {

	public static void main(String[] args) throws Exception {
		//In all variants mode, variants outside confidence regions are printed
		boolean allVariants = "-a".equals(args[0]);
		int iA=0;
		if(allVariants) iA++;
		String referenceGenomeFile = args[iA++];
		String referenceRegionsFile = args[iA++];
		String variantsFile = args[iA++];
		String strsFile = args[iA++];
		ReferenceGenome genome = new ReferenceGenome(referenceGenomeFile);
		GenomicRegionComparator comparator = new GenomicRegionComparator(genome.getSequencesMetadata());
		SimpleGenomicRegionFileHandler fh = new SimpleGenomicRegionFileHandler();
		List<GenomicRegion> referenceRegions = fh.loadRegions(referenceRegionsFile);
		List<GenomicRegion> strs = fh.loadRegions(strsFile);
		try (VCFFileReader vcfReader = new VCFFileReader(variantsFile)) {
			VCFFileHeader header = vcfReader.getHeader();
			VCFFileWriter writer = new VCFFileWriter();
			writer.printHeader(header, System.out);
			GenomicRegion overlappingRegion = null;
			int i=0;
			int k=0;
			Iterator<VCFRecord> it = vcfReader.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				boolean recordFullyCovered = false;
				//Process first last overlapping region
				if(overlappingRegion!=null) {
					int cmp = comparator.compare(overlappingRegion, record);
					if(cmp<-1) {
						//Current record is after the last overlap region
						VCFRecord refRecord = makeReferenceRecord(genome,header,overlappingRegion);
						writer.printVCFRecord(refRecord, System.out);
						overlappingRegion = null;
					} else if (cmp<2) {
						int spanVar = Math.max(record.length(), record.getLast()-record.getFirst()+1);
						if(overlappingRegion.getFirst() <= record.getFirst()-spanVar) {
							//Print confidence region before variant
							GenomicRegion regBefore = new GenomicRegionImpl(overlappingRegion.getSequenceName(), overlappingRegion.getFirst(), record.getFirst()-spanVar);
							VCFRecord refRecord = makeReferenceRecord(genome,header,regBefore);
							writer.printVCFRecord(refRecord, System.out);
						}
						recordFullyCovered = overlappingRegion.getFirst()<=record.getFirst() && record.getLast()<=overlappingRegion.getLast();
						if(record.getLast()+spanVar<=overlappingRegion.getLast()) {
							//Start the next confidence region a number of bp right of the variant that corresponds with its potential span
							overlappingRegion = new GenomicRegionImpl(overlappingRegion.getSequenceName(), record.getLast()+spanVar, overlappingRegion.getLast());
						} else {
							overlappingRegion = null;
						}
					}
				}
				if(overlappingRegion==null) {
					for(;i<referenceRegions.size();i++) {
						GenomicRegion refRegion = referenceRegions.get(i);
						int cmp = comparator.compare(refRegion, record);
						if(cmp<-1) {
							//Current record is after this confidence region
							VCFRecord refRecord = makeReferenceRecord(genome,header,refRegion);
							writer.printVCFRecord(refRecord, System.out);
						} else if (cmp<2) {
							int spanVar = Math.max(record.length(), record.getLast()-record.getFirst()+1);
							if(refRegion.getFirst() <= record.getFirst()-spanVar) {
								//Print confidence region before variant
								GenomicRegion regBefore = new GenomicRegionImpl(refRegion.getSequenceName(), refRegion.getFirst(), record.getFirst()-spanVar);
								VCFRecord refRecord = makeReferenceRecord(genome,header,regBefore);
								writer.printVCFRecord(refRecord, System.out);
							}
							recordFullyCovered = recordFullyCovered || refRegion.getFirst()<=record.getFirst() && record.getLast()<=refRegion.getLast();
							if(record.getLast()+spanVar<=refRegion.getLast()) {
								//Start the next confidence region a number of bp right of the variant that corresponds with its potential span
								overlappingRegion = new GenomicRegionImpl(refRegion.getSequenceName(), record.getLast()+spanVar, refRegion.getLast());
							}
						} else {
							break;
						}
					}
				}
				//Check intersection with STRs
				GenomicRegion overlappingSTR = null;
				for(;k<strs.size();k++) {
					GenomicRegion str = strs.get(k);
					int cmp = comparator.compare(str, record);
					if(cmp>=-1) {
						if(cmp<2) overlappingSTR = str;
						break;
					}
				}
				if(overlappingSTR!=null) {
					GenomicVariant variant = record.getVariant();
					if(variant.isSNV()) variant.setType(GenomicVariant.TYPE_EMBEDDED_SNV);
					else variant.setType(GenomicVariant.TYPE_STR);
					record.addAnnotation(new GenomicVariantAnnotation(variant, "STRLEN", overlappingSTR.length()));
				}
				if(allVariants || recordFullyCovered) writer.printVCFRecord(record, System.out);
			}
			for(;i<referenceRegions.size();i++) {
				GenomicRegion refRegion = referenceRegions.get(i);
				VCFRecord refRecord = makeReferenceRecord(genome,header,refRegion);
				writer.printVCFRecord(refRecord, System.out);
			}
		}
		

	}

	private static VCFRecord makeReferenceRecord(ReferenceGenome genome, VCFFileHeader header, GenomicRegion refRegion) {
		CharSequence reference = genome.getReference(refRegion);
		List<String> alleles = new ArrayList<>();
		alleles.add(reference.toString());
		GenomicVariantImpl variant = new GenomicVariantImpl(refRegion.getSequenceName(), refRegion.getFirst(), alleles);
		CalledGenomicVariantImpl refCall = new CalledGenomicVariantImpl(variant, 0);
		return new VCFRecord(variant, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL, refCall, header);
	}

}
