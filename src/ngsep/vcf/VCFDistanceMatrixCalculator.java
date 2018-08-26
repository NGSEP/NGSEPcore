package ngsep.vcf;

import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.clustering.DistanceMatrix;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.VariantCallReport;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;


public class VCFDistanceMatrixCalculator {

	public static final int DISTANCE_SOURCE_GENOTYPES_SIMPLE=0;
	public static final int DISTANCE_SOURCE_GENOTYPES_COPY_NUMBER=1;
	public static final int DISTANCE_SOURCE_COPY_NUMBER=2;
	public static final int DISTANCE_SOURCE_ALLELE_DEPTH=3;
	
	private Logger log = Logger.getLogger(VCFDistanceMatrixCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private int ploidy = 2;
	private int matrixType = 0;
	private int distanceSource = DISTANCE_SOURCE_GENOTYPES_SIMPLE;

	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}
	
	
	 public static void main (String [ ] args) throws Exception {

		VCFDistanceMatrixCalculator dmCalculator = new VCFDistanceMatrixCalculator();
		//Parameters
		int k=CommandsDescriptor.getInstance().loadOptions(dmCalculator, args);
		
		String vcfFile = args[k++];
		DistanceMatrix dm = dmCalculator.generateMatrix(vcfFile);
		dm.printMatrix(System.out);
		
	 }
	
	 /**
	  * Process a vcf file to generate a SNVs distance matrix.
	  * @param vcfFile VCF filename.
	  * @throws IOException
	  */
	 public DistanceMatrix generateMatrix (String vcfFile) throws IOException {
		 List<String> samples;
		 double distanceMatrix[][];
		 
		 int numSamples;
		 try (VCFFileReader vcfFileReader = new VCFFileReader(vcfFile)) {
			Iterator<VCFRecord> iteratorRecords = vcfFileReader.iterator();
			samples = vcfFileReader.getHeader().getSampleIds();
			numSamples = samples.size();
			distanceMatrix = new double[numSamples][numSamples];
			int genotypePerSamplesComparison[][] = new int[numSamples][numSamples];
			
			int n = ploidy;
			if(n<2){
				n = 2;
			}
			
			
			float ploidyLevels[] = new float[n+1];
			//generate ploidy range for individual from dosage data
			for(int y=0; y <= n;y++){
				ploidyLevels[y] = (1.0f/n) * y;
			}
			
			//Iterate over every variant in VCF file
			while(iteratorRecords.hasNext()){
				VCFRecord vcfRecord = iteratorRecords.next();
				GenomicVariant var = vcfRecord.getVariant();
				String [] alleles = var.getAlleles();
					
				List<CalledGenomicVariant> genotypeCalls = vcfRecord.getCalls();
				float numericGenotypes[] = new float[genotypeCalls.size()];
				Arrays.fill(numericGenotypes, CalledSNV.GENOTYPE_UNDECIDED);
		    	//Calculate dosage for each sample
		    	for (int i=0;i<genotypeCalls.size();i++) {
		    		CalledGenomicVariant call = genotypeCalls.get(i);
		    		if(call.isUndecided()) continue;
		    		if(distanceSource == DISTANCE_SOURCE_GENOTYPES_SIMPLE) {
		    			byte [] idxCalledAlleles = call.getIndexesCalledAlleles();
		    			//TODO: Improve for heterozygous in multiallelic
		    			if (idxCalledAlleles.length==1) numericGenotypes[i] = idxCalledAlleles[0];
		    			else numericGenotypes[i] = (idxCalledAlleles[0]+idxCalledAlleles[1])/alleles.length;
		    		} else if(distanceSource == DISTANCE_SOURCE_GENOTYPES_COPY_NUMBER) {
		    			byte [] acn = call.getAllelesCopyNumber();
		    			numericGenotypes[i] = 0;
		    			for(int j=0;j<acn.length;j++) {
		    				numericGenotypes[i]+=j*acn[j];
		    			}
		    			numericGenotypes[i]/=2.0;
		    		} else if(distanceSource == DISTANCE_SOURCE_COPY_NUMBER) {
		    			numericGenotypes[i] = call.getCopyNumber();
		    		} else if(distanceSource == DISTANCE_SOURCE_ALLELE_DEPTH) {
		    			if(!var.isBiallelic()) continue;
		    			VariantCallReport report = call.getCallReport();
		    			if(report == null) continue;
		    			float countRef = report.getCount(alleles[0]);
			    		float countAlt = report.getCount(alleles[1]);
			    		//Depends of ploidy assign a value to dosage
			    		if((countRef + countAlt) > 0){
			    			float dosage = countRef / (countRef + countAlt);
					    	numericGenotypes[i] = roundToArray(dosage, ploidyLevels);
		    			}
		    		}
		    	}
		    	
		    	//calculate distance samples x samples for all SNVs
		    	for(int j=0;j<numSamples;j++){
		    		for(int k=0;k<numSamples;k++){
		    			if(numericGenotypes[j]==CalledSNV.GENOTYPE_UNDECIDED || numericGenotypes[k]==CalledSNV.GENOTYPE_UNDECIDED ) continue;
		    			//distance between pair of genotypes for a single variant
	    			    distanceMatrix[j][k] += Math.abs(numericGenotypes[j]-numericGenotypes[k]);
	    				//matrix needed to save value by how divide
	    				genotypePerSamplesComparison[j][k]++;
			    	}
		    	}		    	
			}
			//Normalize genetic distance value depending number of samples x samples per Variant found genotyped (Omit missing values)
			for(int j=0;j<numSamples;j++){
	    		for(int k=0;k<numSamples;k++){
					if(genotypePerSamplesComparison[j][k] > 0){
						distanceMatrix[j][k] = distanceMatrix[j][k]/genotypePerSamplesComparison[j][k];
					}
		    	}
	    	}
		}
		
		
		
		
		DistanceMatrix dMatrix = new DistanceMatrix(samples, distanceMatrix);
		dMatrix.setMatrixOutputType(matrixType);
		
		return dMatrix;
	
	}
	 

	/**
	  * Assign real dosage value depending of ploidy.
	  * @param value Dosage ratio calculate from raw base pair count.
	  * @param array Ploidy levels array to fit value.
	  * @throws IOException
	*/
	public float roundToArray(float value, float [] array){
		 
		 float best_abs = 1.0f;
		 float rounded = value;
		 
		 for(int i=0; i < array.length; i++){
			 if(best_abs > Math.abs(value-array[i]))
			 {
				 best_abs = value-array[i];
				 rounded = array[i];
			 }		 
		 }
		 return rounded;
	 }
	

	
	public int getPloidy() {
		return ploidy;
	}
	public void setPloidy(int ploidy) {
		this.ploidy = ploidy;
	}
	public void setPloidy(Integer ploidy) {
		this.setPloidy(ploidy.intValue());
	}	
	
	public int getMatrixType() {
		return matrixType;
	}
	public void setMatrixType(int matrixType) {
		this.matrixType = matrixType;
	}
	public void setMatrixType(Integer matrixType) {
		this.setMatrixType(matrixType.intValue());
	}
	public int getDistanceSource() {
		return distanceSource;
	}
	public void setDistanceSource(int distanceSource) {
		this.distanceSource = distanceSource;
	}
	public void setDistanceSource(Integer distanceSource) {
		this.setDistanceSource(distanceSource.intValue());
	}
}
