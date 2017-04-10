package ngsep.vcf;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.DistanceMatrix;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;


public class VCFDistanceMatrixCalculator {

	private Logger log = Logger.getLogger(VCFDistanceMatrixCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private int ploidy = 2;
	private int matrixType = 0;

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
	
	
	 public static void main (String [ ] args) throws IOException {

		VCFDistanceMatrixCalculator dmCalculator = new VCFDistanceMatrixCalculator();
		
		if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")){
			CommandsDescriptor.getInstance().printHelp(VCFDistanceMatrixCalculator.class);
			return;
		}
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
	 public DistanceMatrix generateMatrix (String vcfFile) throws IOException{

		VCFFileReader vcfFileReader = new VCFFileReader(vcfFile);
		Iterator<VCFRecord> iteratorRecords = vcfFileReader.iterator();
		List<Sample> samples = vcfFileReader.getHeader().getSamples();
		float distanceMatrix[][] = new float[samples.size()][samples.size()];
		int genotypePerSamplesComparison[][] = new int[samples.size()][samples.size()];
		
		if(ploidy<2){
			ploidy = 2;
		}
		
		
		float ploidyLevels[] = new float[ploidy+1];
		//generate ploidy range for individual from dosage data
		for(int y=0; y <= ploidy;y++){
			ploidyLevels[y] = (1.0f/ploidy) * y;
		}
		
		//Iterate over every variant in VCF file
		while(iteratorRecords.hasNext()){
			VCFRecord vcfRecord = iteratorRecords.next();
			//The variant is a SNV
			if(!(vcfRecord.getVariant() instanceof SNV)){
				//Process only biallelic SNVs
				continue;
			}
				
			List<CalledGenomicVariant> genotypeCalls = vcfRecord.getCalls();
			float genotypes[] = new float[genotypeCalls.size()];
	    	//Calculate dosage for each sample
	    	for (int i=0;i<genotypeCalls.size();i++) {
	    		CalledGenomicVariant call = genotypeCalls.get(i);
	    		if(!(call instanceof CalledSNV)) {
	    			continue;
	    		}
	    		CalledSNV snv = (CalledSNV) call;
	    		genotypes[i] = snv.getGenotype();
	    		if(genotypes[i] != CalledSNV.GENOTYPE_UNDECIDED){
	    			float countRef = snv.getCountReference();
		    		float countAlt = snv.getCountAlternative();
		    		//Depends of ploidy assign a value to dosage
		    		if((countRef + countAlt) == 0){
	    				genotypes[i]=CalledSNV.GENOTYPE_UNDECIDED;
	    			} else {
	    				
	    				if(ploidy != 2){
	    					float dosage = countRef / (countRef + countAlt);
				    	    genotypes[i] = roundToArray(dosage, ploidyLevels);
	    				} else {
	    					if(genotypes[i]==CalledSNV.GENOTYPE_HOMOREF){
		    					genotypes[i] = 1.0f;
		    				}else if(genotypes[i]==CalledSNV.GENOTYPE_HOMOALT){
		    					genotypes[i] = 0.0f;
		    				} else if(genotypes[i]==CalledSNV.GENOTYPE_HETERO){
		    					genotypes[i] = 0.5f;
		    				}
	    				}

	    			}
		    	    
	    		}
	    	}
	    	
	    	//calculate distance samples x samples for all SNVs
	    	for(int j=0;j<samples.size();j++){
	    		for(int k=0;k<samples.size();k++){
	    			if(!(genotypes[j]==CalledSNV.GENOTYPE_UNDECIDED || 
	    				genotypes[k]==CalledSNV.GENOTYPE_UNDECIDED )){
	    				//distance between pair of genotypes for a single SNV
	    				
	    			    if(ploidy != 2){
	    			    	distanceMatrix[j][k] += Math.abs(genotypes[j]-genotypes[k]);
	    			    } else {
	    			    	if(j==k){
	    						distanceMatrix[j][k] += 1;
	    					} else if(genotypes[j]==genotypes[k] && genotypes[k]!=0.5){
		    					distanceMatrix[j][k] += 1;
		    				} else if(genotypes[j]==0.0 && genotypes[k] == 1.0 ||
		    						genotypes[j]==1.0 && genotypes[k] == 0.0 ){
		    					distanceMatrix[j][k] += 0;
		    				} else if(genotypes[j] == 0.5 || genotypes[k]== 0.5){
		    					distanceMatrix[j][k] += 0.5;
		    				}
	    			    }

	    				//matrix need to save value by how divide
	    				genotypePerSamplesComparison[j][k]++;
	    			}
		    	}
	    	}		    	
		}
		
		//Normalize genetic distance value depending number of samples x samples per Variant found genotyped (Omit missing values)
		for(int j=0;j<samples.size();j++){
    		for(int k=0;k<samples.size();k++){
    			
    				if(genotypePerSamplesComparison[j][k] > 0){
    					if(ploidy!=2){
    						distanceMatrix[j][k] = distanceMatrix[j][k]/genotypePerSamplesComparison[j][k];
    					} else {
    						distanceMatrix[j][k] = 1 - (distanceMatrix[j][k]/genotypePerSamplesComparison[j][k]);
    					}

    				}
		
	    	}
    	}
		
		DistanceMatrix dMatrix = new DistanceMatrix(samples, distanceMatrix);
		dMatrix.setMatrixType(matrixType);
		vcfFileReader.close();
		
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
	
	
}
