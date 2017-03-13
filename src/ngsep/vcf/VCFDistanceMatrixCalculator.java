package ngsep.vcf;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.DistanceMatrix;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;


public class VCFDistanceMatrixCalculator {

	private Logger log = Logger.getLogger(VCFDistanceMatrixCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private int ploidy = 2;

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
		dmCalculator.printMatrix(dm, System.out);
		
	 }
	
	 /**
		 * Process a vcf file to generate a SNVs distance matrix.
		 * @param vcfFile VCF filename.
		 * @throws IOException
	 */
	 public DistanceMatrix generateMatrix (String vcfFile) throws IOException{

		VCFFileReader vcfFileReader = new VCFFileReader(vcfFile);
		Iterator<VCFRecord> iteratorRecords = vcfFileReader.iterator();
		List<String> samples = vcfFileReader.getSampleIds();
		float distanceMatrix[][] = new float[samples.size()][samples.size()];
		int genotypePerSamplesComparison[][] = new int[samples.size()][samples.size()];
		
		
		float ploidyLevels[] = new float[ploidy];
		
		//generate ploidy range for individual from dosage data
		for(int y=0; y < ploidy;y++){
			ploidyLevels[y] = (1/ploidy) * y;
		}
		
		//Iterate over every variant in VCF file
		while(iteratorRecords.hasNext()){
			VCFRecord vcfRecord = iteratorRecords.next();
			//The variant is a SNV
			if(vcfRecord.getVariant().isSNV()){
				
				List<CalledGenomicVariant> genotyped = vcfRecord.getCalls();
				float genotypes[] = new float[genotyped.size()];
		    	//Calculate dosage for each sample
		    	for (int i=0;i<genotyped.size();i++) {
		    		CalledSNV snv = (CalledSNV) genotyped.get(i);
		    		genotypes[i] = snv.getGenotype();
		    		if(genotypes[i] == 1){
		    			int countRef = snv.getCountReference();
			    		int countAlt = snv.getCountAlternative();
			    		
			    		//Depends of ploidy assign a value to dosage
			    	    float dosage = countRef / (countRef + countAlt);
			    	    genotypes[i] = roundToArray(dosage, ploidyLevels);
			    		
		    		}
		    	}
		    	
		    	//calculate distance samples x samples for all SNVs
		    	for(int j=0;j<samples.size();j++){
		    		for(int k=0;k<samples.size();k++){
		    			if(!(genotypes[j]==-1||genotypes[k]==-1)){
		    				//distance between pair of genotypes for a single SNV
		    				distanceMatrix[j][k] += (float)Math.abs(genotypes[j]-genotypes[k]);
		    				//matrix need to save value by how divide
		    				genotypePerSamplesComparison[j][k]++;
		    			}
			    	}
		    	}		    	
				
			}
		}
		
		//Normalize genetic distance value depending number of samples x samples per Variant found genotyped (Omit missing values)
		for(int j=0;j<samples.size();j++){
    		for(int k=0;k<samples.size();k++){
    			distanceMatrix[j][k] = distanceMatrix[j][k]/genotypePerSamplesComparison[j][k];
	    	}
    	}
		
		DistanceMatrix dMatrix = new DistanceMatrix(samples, distanceMatrix);
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
		 
		 for(int i=0; i < array.length; i++){
			 if(best_abs > Math.abs(value-array[i]))
			 {
				 best_abs = array[i];
			 }		 
		 }
		 
		 return best_abs;
	 }
	
	/**
	  * Print distance matrix.
	  * @param out matrix in generic format.
	  * @throws IOException
	*/
	public void printMatrix (DistanceMatrix dm, PrintStream out) {
		//print number of samples of the matrix
	    out.println(dm.getnSamples());
	    // print samples x samples distance matrix
	    for(int j=0;j<dm.getnSamples();j++){
	    	String row = "";
    		for(int k=0;k<dm.getnSamples();k++){
    			row += dm.getDistanceMatrix()[j][k];
    			row += " ";
	    	}
    		
	    	out.println(dm.getSamples().get(j)+" "+row);
    	}

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
	
	
	
}
