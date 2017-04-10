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
	private int matrixType;
	
	/**
	 * Construct a DistanceMatrix object from a file which represent a matrix in a generic format.
	 * @param filename matrix filename.
	 * @throws IOException
	 * @throws NumberFormatException
	*/
	public DistanceMatrix(String filename) throws IOException, NumberFormatException{
		BufferedReader br = new BufferedReader(new FileReader(filename));
		samples = new ArrayList<Sample>();
		try {
			this.setnSamples(Integer.parseInt(br.readLine()));
			this.setDistanceMatrix(new float[this.getnSamples()][this.getnSamples()]);
			
		    String matrixRow = br.readLine();
		    int row = 0;
		    while (matrixRow != null) {
		        
		    	String[] matrixCell  = matrixRow.split("\\s+");
		    	this.getSamples().add(new Sample(matrixCell[0]));
		    	for(int column = 1; column < (this.getnSamples()+1); column++){
		    		this.getDistanceMatrix()[row][column-1]=Float.parseFloat(matrixCell[column]);
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
    		
	    	out.println(this.getSamples().get(j)+" "+row);
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
		return matrixType;
	}

	public void setMatrixType(int matrixType) {
		this.matrixType = matrixType;
	}
	
	
}
