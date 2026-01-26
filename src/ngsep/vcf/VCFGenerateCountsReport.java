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

import java.io.PrintStream;
import java.util.Iterator;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.VariantCallReport;

/**
 * Extracts allele counts for biallelic variants from a VCF file.
 * The format of the counts is REF/ALT for each sample
 * @author Jorge Duitama
 */
public class VCFGenerateCountsReport {

	public static void main(String[] args) throws Exception {
		PrintStream out = System.out;
		try (VCFFileReader reader = new VCFFileReader(args[0])) {
			VCFFileHeader header = reader.getHeader();
			List<String> sampleIds = header.getSampleIds();
			out.print("Chr\tPos\tRef\tAlt");
			for(String sampleId:sampleIds) out.print("\t"+sampleId);
			out.println();
			Iterator<VCFRecord> it = reader.iterator();
			while (it.hasNext()) {
				VCFRecord record = it.next();
				GenomicVariant var = record.getVariant();
				if(!var.isBiallelic()) continue;
				String ref = var.getReference();
				String alt = var.getAlleles()[1];
				out.print(var.getSequenceName()+"\t"+var.getFirst()+"\t"+var.getReference()+"\t"+var.getAlleles()[1]);
				List<CalledGenomicVariant> calls = record.getCalls();
				for(int i=0;i<calls.size();i++) {
					CalledGenomicVariant call = calls.get(i);
					VariantCallReport report = call.getCallReport();
					out.print("\t"+report.getCount(ref)+"/"+report.getCount(alt));
				}
				out.println();
			}
		}

	}

}
