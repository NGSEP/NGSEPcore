package ngsep.clustering;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class DistanceMatrix {

	private List<String> ids;
	private double distanceMatrix[][];
	private int matrixOutputType;
	
	/**
	 * Construct a DistanceMatrix object from a file which represent a matrix in a generic format.
	 * @param filename matrix filename.
	 * @throws IOException
	 * @throws NumberFormatException
	*/
	public DistanceMatrix(String filename) throws IOException {
		ids = new ArrayList<>();
		int fileMatrixType = 0;
		try (FileReader fr=new FileReader(filename);
			 BufferedReader br = new BufferedReader(fr)){
			int numSamples;
			try {
				numSamples = Integer.parseInt(br.readLine());
			} catch (Exception e) {
				throw new IOException("Number format error reading number of samples",e);
			}
			distanceMatrix = new double[numSamples][numSamples];
			
		    String matrixRow = br.readLine();
		    int row = 0;
		    while (matrixRow != null) {
		        
		    	String[] matrixCell  = matrixRow.split("\\s+");
		    	ids.add(matrixCell[0]);
		    	if(row ==0 && matrixCell.length < (numSamples+1)){
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
	    		int indent = (numSamples+1) - matrixCell.length;
	    		
	    		for(int column = 1; column < matrixCell.length; column++){
	    			double value;
					try {
						value = Double.parseDouble(matrixCell[column]);
					} catch (Exception e) {
						throw new IOException("Number format error at row "+row+" column: "+column+" value: "+matrixCell[column],e);
					}
	    			if(fileMatrixType == 0){ // full matrix
	    				distanceMatrix[row][column-1]=value;
	    			} else if(fileMatrixType == 1){ // lower
	    				distanceMatrix[row][column-1]=value;
	    				distanceMatrix[column-1][row]=value;
		    		} else if(fileMatrixType == 2){ // upper
		    			distanceMatrix[column-1+indent][row]=value;
		    			distanceMatrix[row][column-1+indent]=value;
		    		}
		    	}

		    	
		    	matrixRow = br.readLine();
		    	row++;
		    }
		}
	}
	
	/**
	 * Construct a DistanceMatrix object from two objects: a list with ids and 
	 * an array of double values which represent the distances.
	 * @param ids of the objects with the given distances
	 * @param distanceMatrix Values of distances between the objects
	*/
	public DistanceMatrix(List<String> ids, double distanceMatrix[][] ){
		this.setIds(ids);
		this.setDistanceMatrix(distanceMatrix);
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
	    out.println(this.getNumSamples());
	    // print samples x samples distance matrix
	    for(int j=0;j<distanceMatrix.length;j++){
	    	String row = "";
    		for(int k=0;k<distanceMatrix[j].length;k++){
    			if(this.getMatrixType() == 0 || (this.getMatrixType() == 1 && j>k) || (this.getMatrixType() == 2 && k>j) ) {
    				row += distanceMatrix[j][k];
        			row += " ";
        			
    			} else if(this.getMatrixType() == 2 && j>k){
        			row += " ";
    			}
	
	    	}
	    	out.println(ids.get(j)+" "+row);
    	}
	}
	
	

	/**
	 * @return the ids
	 */
	public List<String> getIds() {
		return ids;
	}

	/**
	 * @param ids the ids to set
	 */
	public void setIds(List<String> ids) {
		this.ids = ids;
	}

	public double[][] getDistanceMatrix() {
		return distanceMatrix;
	}

	public void setDistanceMatrix(double distanceMatrix[][]) {
		this.distanceMatrix = distanceMatrix;
	}

	public int getNumSamples() {
		return distanceMatrix.length;
	}

	public int getMatrixType() {
		return matrixOutputType;
	}

	public void setMatrixOutputType(int matrixOutputType) {
		this.matrixOutputType = matrixOutputType;
	}
	
}
