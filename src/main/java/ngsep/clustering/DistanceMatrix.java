package ngsep.clustering;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class DistanceMatrix {

	public static final int MATRIX_TYPE_FULL = 0;
	public static final int MATRIX_TYPE_LOWER_LEFT = 1;
	public static final int MATRIX_TYPE_UPPER_RIGHT = 2;
	private List<String> ids;
	private double distances[][];
	private int matrixOutputType;
	
	/**
	 * Construct a DistanceMatrix object from a file which represent a matrix in a generic format.
	 * @param filename matrix filename.
	 * @throws IOException
	 */
	public DistanceMatrix(String filename) throws IOException {
		try (FileInputStream fis = new FileInputStream(filename)) {
			loadFromFile(fis);
		}
	}
	/**
	 * Construct a DistanceMatrix object from a file which represent a matrix in a generic format.
	 * @param inputStream to read the matrix data.
	 * @throws IOException
	*/
	public DistanceMatrix(InputStream is) throws IOException {
		loadFromFile(is);
	}
	
	private void loadFromFile(InputStream is) throws IOException {
		ids = new ArrayList<>();
		int fileMatrixType = MATRIX_TYPE_FULL;
		try (InputStreamReader fr=new InputStreamReader(is);
			 BufferedReader br = new BufferedReader(fr)){
			int numSamples;
			try {
				numSamples = Integer.parseInt(br.readLine());
			} catch (Exception e) {
				throw new IOException("Number format error reading number of samples",e);
			}
			distances = new double[numSamples][numSamples];
			
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
		    			fileMatrixType = MATRIX_TYPE_UPPER_RIGHT; //upper
		    		} else {
		    			fileMatrixType = MATRIX_TYPE_LOWER_LEFT; //lower
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
	    			if(fileMatrixType == MATRIX_TYPE_FULL){ // full matrix
	    				distances[row][column-1]=value;
	    			} else if(fileMatrixType == MATRIX_TYPE_LOWER_LEFT){ // lower
	    				distances[row][column-1]=value;
	    				distances[column-1][row]=value;
		    		} else if(fileMatrixType == MATRIX_TYPE_UPPER_RIGHT){ // upper
		    			distances[column-1+indent][row]=value;
		    			distances[row][column-1+indent]=value;
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
		this.setDistances(distanceMatrix);
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
	    for(int j=0;j<distances.length;j++){
	    	String row = "";
    		for(int k=0;k<distances[j].length;k++){
    			if(this.getMatrixType() == MATRIX_TYPE_FULL || (this.getMatrixType() == MATRIX_TYPE_LOWER_LEFT && j>k) || (this.getMatrixType() == MATRIX_TYPE_UPPER_RIGHT && k>j) ) {
    				row += distances[j][k];
        			row += " ";
        			
    			} else if(this.getMatrixType() == MATRIX_TYPE_UPPER_RIGHT && j>k){
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

	public double[][] getDistances() {
		return distances;
	}

	public void setDistances(double distances[][]) {
		this.distances = distances;
	}

	public int getNumSamples() {
		return distances.length;
	}

	public int getMatrixType() {
		return matrixOutputType;
	}

	public void setMatrixOutputType(int matrixOutputType) {
		this.matrixOutputType = matrixOutputType;
	}
	
}
