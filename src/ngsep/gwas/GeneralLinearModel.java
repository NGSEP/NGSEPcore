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
package ngsep.gwas;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import JSci.maths.statistics.FDistribution;
import ngsep.math.FactorialAnova;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;

/**
 * 
 * @author Andrea Parra
 *
 */
public class GeneralLinearModel {

	public static void main(String[] args) throws Exception {
		GeneralLinearModel instance = new GeneralLinearModel();
		instance.process (args[0], args[1], System.out);

	}

	public void process(String vcfFile, String phenotypesFile, PrintStream out) throws IOException {
		List<Double> phenotypes = readPhenotypes(phenotypesFile);
		
			
		try (VCFFileReader reader = new VCFFileReader(vcfFile)){
				
			//Ids of the samples
			List<String> sampleIds = reader.getSampleIds();
			//If you want the full header
			VCFFileHeader header = reader.getHeader();
			
			Iterator<VCFRecord> it = reader.iterator();
			while(it.hasNext()) {
				
				VCFRecord record = it.next();
				//Basic information from the variant
				GenomicVariant variant = record.getVariant();
				
				List<CalledGenomicVariant> genotypeCalls = record.getCalls();
				
				List<String> encodedGenotypes = encodeGenotypes(genotypeCalls);
				
				// Create a new genotype object to run ANOVA.
				FactorialAnova anova = new FactorialAnova(encodedGenotypes, phenotypes);
				printSummary(variant, anova,out);
			}	
		}	
	}

	private List<Double> readPhenotypes(String phenotypesFile) throws IOException {
		
        final String fileToParse = phenotypesFile;
         
        //Delimiter used in CSV file
        final String DELIMITER = ",";
        
        //Initiailze headers and file info.
        ArrayList<String> headers = new ArrayList<String>();
        List<Double> phenotypes = new ArrayList<Double>();
        
        try (FileReader fr = new FileReader(fileToParse);
        	 BufferedReader fileReader = new BufferedReader(fr);)
        {
            String line = "";
            //Create the file reader
            
             
            //Read the file line by line
            int counter = 0;
           
            while ((line = fileReader.readLine()) != null)
            {	
            	//Get all tokens available in line
                String[] tokens = line.split(DELIMITER);
                
                //Save the headers into an ArrayList
            	if(counter == 0) {
            		for(String token : tokens) {
            			headers.add(token);
            		}
            	} else {
            		for(String token : tokens)
                    {	
            			double phenotype = Double.valueOf(token);
                        phenotypes.add(phenotype);
                    }	
            	}
            	counter++;
            }
        }
        
		return phenotypes;
	}

	private List<String> encodeGenotypes(List<CalledGenomicVariant> genotypeCalls) {
		List<String> allelesList = new ArrayList<String>();
		for(int i=0;i<genotypeCalls.size();i++) {
			CalledGenomicVariant genotypeCall = genotypeCalls.get(i);
			String allele = String.join("", genotypeCall.getCalledAlleles());
			allelesList.add(allele);
		}

		return allelesList;

	}

	private void printSummary(GenomicVariant variant, FactorialAnova anova, PrintStream out) {
		double pValue = anova.getpValue();
		FDistribution fdist = anova.getpFdist();
		out.println(variant.getSequenceName() + "\t" + variant.getFirst() + "\t" + variant.getLast() + "\t" + pValue);
		
	}
}
