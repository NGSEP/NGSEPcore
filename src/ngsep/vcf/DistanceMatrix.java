package ngsep.vcf;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class DistanceMatrix {

	private List<String> samples;
	private float distanceMatrix[][];
	private int nSamples;
	
	/**
	 * Construct a DistanceMatrix object from a file which represent a matrix in a generic format.
	 * @param filename matrix filename.
	 * @throws IOException
	 * @throws NumberFormatException
	*/
	public DistanceMatrix(String filename) throws IOException, NumberFormatException{
		BufferedReader br = new BufferedReader(new FileReader(filename));
		samples = new ArrayList<String>();
		try {
			this.setnSamples(Integer.parseInt(br.readLine()));
			this.setDistanceMatrix(new float[this.getnSamples()][this.getnSamples()]);
			
		    String matrixRow = br.readLine();
		    int row = 0;
		    while (matrixRow != null) {
		        
		    	String[] matrixCell  = matrixRow.split("\\s+");
		    	this.getSamples().add(matrixCell[0]);
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
	public DistanceMatrix(List<String> samples,float distanceMatrix[][] ){
		this.setSamples(samples);
		this.setDistanceMatrix(distanceMatrix);
		this.setnSamples(samples.size());
	}

	List<String> getSamples() {
		return samples;
	}

	void setSamples(List<String> samples) {
		this.samples = samples;
	}

	float[][] getDistanceMatrix() {
		return distanceMatrix;
	}

	void setDistanceMatrix(float distanceMatrix[][]) {
		this.distanceMatrix = distanceMatrix;
	}

	int getnSamples() {
		return nSamples;
	}

	void setnSamples(int nSamples) {
		this.nSamples = nSamples;
	}
	
	
}
