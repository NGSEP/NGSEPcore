package ngsep.vcf;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.KmersCounter;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;


public class DistanceMatrixCalculator {

	private Logger log = Logger.getLogger(KmersCounter.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private ArrayList<String> matrix = new ArrayList<>();
	private List<String> samples;
	
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

		DistanceMatrixCalculator dmCalculator = new DistanceMatrixCalculator();
		
		if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")){
			CommandsDescriptor.getInstance().printHelp(DistanceMatrixCalculator.class);
			return;
		}
		//Parameters
		int k=CommandsDescriptor.getInstance().loadOptions(dmCalculator, args);
		
		String vcfFile = args[k++];
		dmCalculator.generateMatrix(vcfFile);
		dmCalculator.printMatrix(System.out);
		
	 }
	
	 /**
		 * Process a vcf file to generate a SNVs distance matrix.
		 * @param vcfFile VCF filename.
		 * @throws IOException
	 */
	 public float[][] generateMatrix (String vcfFile) throws IOException{

		VCFFileReader vcfFileReader = new VCFFileReader(vcfFile);
		Iterator<VCFRecord> iteratorRecords = vcfFileReader.iterator();
		samples = vcfFileReader.getSampleIds();
		float distanceMatrix[][] = new float[samples.size()][samples.size()];
		int genotypePerSamplesComparison[][] = new int[samples.size()][samples.size()];
		
		
		float ploidyLevels[] = new float[ploidy];
		
		//generate ploidy range for individual from dosage data
		for(int y=0; y < ploidy;y++){
			ploidyLevels[y] = (ploidy/1) * y;
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
			String row = "";
    		for(int k=0;k<samples.size();k++){
    			distanceMatrix[j][k] = distanceMatrix[j][k]/genotypePerSamplesComparison[j][k];
    			
    			row+=distanceMatrix[j][k];
    			row+=" ";
    			
	    	}
    		matrix.add(row);
    	}

		vcfFileReader.close();
		
		return distanceMatrix;
	
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
	public void printMatrix (PrintStream out) {
		//print number of samples of the matrix
	    out.println(matrix.size());
	    // print samples x samples distance matrix
	    for(int p=0;p<matrix.size();p++){
	    	out.println(samples.get(p)+" "+matrix.get(p));
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
