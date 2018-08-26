package ngsep.clustering;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;


public class NeighborJoining {

	private Logger log = Logger.getLogger(NeighborJoining.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private DistanceMatrix distanceMatrix;
	
	
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
			
	 	NeighborJoining nj = new NeighborJoining();
		//Parameters
		int k=CommandsDescriptor.getInstance().loadOptions(nj, args);
		
		String matrixFile = args[k++];
	 	DistanceMatrix dm = new DistanceMatrix(matrixFile);
	 	nj.loadMatrix(dm);
		Dendrogram njTree = nj.constructNJTree();
		njTree.printTree(System.out);	
	}
	

	 /**
	  * Load Distance matrix.
	  * @param DistanceMatrix Distance matrix loading , will be used to construct tree.
	*/
	public void loadMatrix (DistanceMatrix distanceMatrix){
		this.distanceMatrix = distanceMatrix;
	}
	
	
	/**
	  * Construct NJ tree and make newick out structure.
	*/
	public Dendrogram constructNJTree(){
	
		double iterableMatrix[][] = distanceMatrix.getDistanceMatrix();
		
		int nSamples = distanceMatrix.getNumSamples();
		
		int nodesToAssign = nSamples;
		
		Map<String, Dendrogram> nodesMerged = new HashMap<>();
		
		Dendrogram njTree = new Dendrogram(null, null, "");
		
		ArrayList<String> nodesList = new ArrayList<>(distanceMatrix.getIds());
		
		
		while(nodesToAssign>2){
			
			double njMatrix[][]= new double[nodesToAssign][nodesToAssign];
			double totalDistance[]= new double[nodesToAssign];
			
			//Calculate total distance array
			for(int i=0;i < nodesToAssign; i++){
				for(int j=0;j < nodesToAssign; j++){
					totalDistance[i] += iterableMatrix[i][j];
				}  
			}
			
			//Compute Neighbor Joining matrix
			for(int i=0;i < nodesToAssign; i++){
				for(int j=0;j < i; j++){
					njMatrix[i][j] = njMatrix[j][i] = ((nodesToAssign - 2) * iterableMatrix[i][j]) - totalDistance[i] - totalDistance[j];
				}  
			}
			
			int rowMin = -1;
			int colMin = -1;
			double minValue = Float.POSITIVE_INFINITY;
			
			//Minimum value in matrix
			//need to be improve searching only in the triangle
			for(int i=0;i < nodesToAssign; i++){
				for(int j=0;j < nodesToAssign; j++){
					if(njMatrix[i][j]<minValue && i != j){
						rowMin = i;
						colMin = j;
						minValue = njMatrix[i][j];
					}
				}  
			}
			
						
			nodesToAssign--;
			
			if(nodesToAssign > 2){
				//branch length estimation
				double leftTreeDistance = (0.5f * iterableMatrix[rowMin][colMin]) + ((totalDistance[rowMin] - totalDistance[colMin] ) / (2 * (nodesToAssign - 1)));
				double rightTreeDistance  = (0.5f * iterableMatrix[rowMin][colMin]) + ((totalDistance[colMin] - totalDistance[rowMin])/ (2 * (nodesToAssign - 1)));

				//Central node join and newick formatting		
				String nodeMergeId = "("+nodesList.get(rowMin) +":"+leftTreeDistance+","+ nodesList.get(colMin)+":"+rightTreeDistance+")";
				
				//Try to get back existing trees
				Dendrogram leftTree = nodesMerged.get(nodesList.get(rowMin)); 
				Dendrogram rightTree = nodesMerged.get(nodesList.get(colMin)); 
				
				//If trees not exist create new ones
				if(leftTree == null){
					leftTree = new Dendrogram(null, null, nodesList.get(rowMin));
				}
				
				if(rightTree == null){
					rightTree = new Dendrogram(null, null, nodesList.get(colMin));
				}
				
				//Make a new tree with the two joined nodes
				njTree = new Dendrogram( leftTree, rightTree, nodeMergeId);			
				nodesMerged.put(nodeMergeId, njTree);
				
				//Next iteration nodes positions
				ArrayList<String> currentNodes = new ArrayList<>();
				currentNodes.add(nodeMergeId);
				for(int n=0;n<nodesList.size();n++){
					if(n != rowMin && n != colMin){ //erase joined nodes from nodes list
						currentNodes.add(nodesList.get(n));
					}

				}
				
				// distance matrix update
				double updateDistanceMatrix[][] = new double[nodesToAssign][nodesToAssign];
				
				//adding new row, col and bringing back old distances
				for(int i=1;i<nodesToAssign;i++){
					for(int j=0;j<i;j++){
						int index_i = nodesList.indexOf(currentNodes.get(i));
						if(j==0){ //if new node, calculate distance
							updateDistanceMatrix[j][i] = updateDistanceMatrix[i][j] = 0.5f * (iterableMatrix[rowMin][index_i]+iterableMatrix[colMin][index_i] - iterableMatrix[rowMin][colMin]);
						} else { //else bring old ones
							int index_j = nodesList.indexOf(currentNodes.get(j));
							updateDistanceMatrix[j][i] = updateDistanceMatrix[i][j] = iterableMatrix[index_i][index_j];
						}
						
					}
					
				}
				
			
				nodesList= currentNodes;
				
				iterableMatrix = updateDistanceMatrix;
				
			} else { // FINAL JOIN ---------------------------------------------------------
				
				
				double leftTreeDistance = (0.5f * iterableMatrix[1][2]) + (0.5f * (totalDistance[1] - totalDistance[2]));
				double rightTreeDistance  = (0.5f * iterableMatrix[1][2]) + (0.5f * (totalDistance[2] - totalDistance[1]));
				
				double centralTreeDistance = 0.0f;
				
				//Calculate final distance between join node and the last two
				if(rowMin!=0){
					if(rowMin<colMin || colMin == 0){
						centralTreeDistance = (0.5f * iterableMatrix[0][rowMin]) + (0.5f * (totalDistance[rowMin] - totalDistance[0]));
					} else{
						centralTreeDistance = (0.5f * iterableMatrix[0][colMin]) + (0.5f * (totalDistance[colMin] - totalDistance[0]));
					}
				} else{
					centralTreeDistance = (0.5f * iterableMatrix[0][colMin]) + (0.5f * (totalDistance[colMin] - totalDistance[0]));
				}
				
							
				String nodeMergeId = nodesList.get(0);
				
				// Newick format
				Dendrogram leftTree = new Dendrogram(null, null, "("+nodesList.get(2)+":"+rightTreeDistance+","+nodesList.get(1)+":"+leftTreeDistance+",");
				Dendrogram rightTree = new Dendrogram(null, null, ");\n");
			
				njTree = new Dendrogram(leftTree, rightTree, nodeMergeId+":"+centralTreeDistance);			

			}
			
		}

		
		return njTree;
	}

	
}

