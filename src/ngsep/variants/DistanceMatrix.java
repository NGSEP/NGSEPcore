package ngsep.variants;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class DistanceMatrix {

	private List<Sample> samples;
	private float distanceMatrix[][];
	private int nSamples;
	private int matrixOutputType;
	
	/**
	 * Construct a DistanceMatrix object from a file which represent a matrix in a generic format.
	 * @param filename matrix filename.
	 * @throws IOException
	 * @throws NumberFormatException
	*/
	public DistanceMatrix(String filename) throws IOException, NumberFormatException{
		BufferedReader br = new BufferedReader(new FileReader(filename));
		samples = new ArrayList<Sample>();
		int fileMatrixType = 0;
		try {
			this.setnSamples(Integer.parseInt(br.readLine()));
			this.setDistanceMatrix(new float[this.getnSamples()][this.getnSamples()]);
			
		    String matrixRow = br.readLine();
		    int row = 0;
		    while (matrixRow != null) {
		        
		    	String[] matrixCell  = matrixRow.split("\\s+");
		    	this.getSamples().add(new Sample(matrixCell[0]));
		    	

		    	
		    	if(row ==0 && matrixCell.length < (this.getnSamples()+1)){
		    		fileMatrixType = 3; //needs to determinate which type of matrix, upper or lower triangle
		    	}
		    	
		    	if(fileMatrixType == 3){
		    		if(matrixCell.length > 2){ // row three of matrix have more than 2 columns ?
		    			fileMatrixType = 2; //upper
		    		} else {
		    			fileMatrixType = 1; //lower
		    		}
		    	}
		    	
		    	
		    	//indent to convert build upper matrix to lower
	    		int indent = (this.getnSamples()+1) - matrixCell.length;
	    		
	    		for(int column = 1; column < matrixCell.length; column++){
	    			if(fileMatrixType == 0){ // full matrix
	    				this.getDistanceMatrix()[row][column-1]=Float.parseFloat(matrixCell[column]);
	    			} else if(fileMatrixType == 1){ // lower
	    				this.getDistanceMatrix()[row][column-1]=Float.parseFloat(matrixCell[column]);
		    			this.getDistanceMatrix()[column-1][row]=Float.parseFloat(matrixCell[column]);
		    		} else if(fileMatrixType == 2){ // upper
		    			this.getDistanceMatrix()[column-1+indent][row]=Float.parseFloat(matrixCell[column]);
		    			this.getDistanceMatrix()[row][column-1+indent]=Float.parseFloat(matrixCell[column]);
		    		}
		    	}

		    	
		    	matrixRow = br.readLine();
		    	row++;
		    }
		} finally {
		    br.close();
		}
	}
	
	/**
	 * Construct a DistanceMatrix object from two objects: a list with the samples and 
	 * an array of float values which represent the distances.
	 * @param samples Identifier of the sample in the matrix.
	 * @param distanceMatrix Values of distances between the samples
	*/
	public DistanceMatrix(List<Sample> samples,float distanceMatrix[][] ){
		this.setSamples(samples);
		this.setDistanceMatrix(distanceMatrix);
		this.setnSamples(samples.size());
	}
	
	/**
	  * Print distance matrix.
	  * matrixType:
	  *  0 = full matrix
	  *  1 = Lower-left matrix
	  *  2 = Upper-right matrix
	  * @param out matrix in generic format.
	  * @throws IOException
	*/
	public void printMatrix (PrintStream out) {
		//print number of samples of the matrix
	    out.println(this.getnSamples());
	    // print samples x samples distance matrix
	    for(int j=0;j<this.getnSamples();j++){
	    	String row = "";
    		for(int k=0;k<this.getnSamples();k++){
    			if(this.getMatrixType() == 0 || (this.getMatrixType() == 1 && j>k) || (this.getMatrixType() == 2 && k>j) ){
    				
    				row += this.getDistanceMatrix()[j][k];
        			row += " ";
        			
    			} else if(this.getMatrixType() == 2 && j>k){
        			row += " ";
    			}
	
	    	}
    		
	    	out.println(this.getSamples().get(j).getId()+" "+row);
    	}

	}
	

	

	public List<Sample> getSamples() {
		return samples;
	}

	public void setSamples(List<Sample> samples) {
		this.samples = samples;
	}

	public float[][] getDistanceMatrix() {
		return distanceMatrix;
	}

	public void setDistanceMatrix(float distanceMatrix[][]) {
		this.distanceMatrix = distanceMatrix;
	}

	public int getnSamples() {
		return nSamples;
	}

	public void setnSamples(int nSamples) {
		this.nSamples = nSamples;
	}

	public int getMatrixType() {
		return matrixOutputType;
	}

	public void setmatrixOutputType(int matrixOutputType) {
		this.matrixOutputType = matrixOutputType;
	}
	
}
