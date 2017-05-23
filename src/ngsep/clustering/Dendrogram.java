package ngsep.clustering;

import java.io.PrintStream;



public class Dendrogram {

  private Dendrogram leftDendro;
  private Dendrogram rightDendro;
  //private double leftTreeDistance;
  //private double rightTreeDistance;
  private String label;
  
  /*
  Tree(Tree leftTree, Tree rightTree, double leftTreeDistance, double rightTreeDistance, String label) {
    this.leftTree = leftTree;
    this.rightTree = rightTree;
    this.leftTreeDistance = leftTreeDistance;
    this.rightTreeDistance = rightTreeDistance;
    this.label = label;
  }*/
  
  Dendrogram(Dendrogram leftDendro, Dendrogram rightDendro,  String label) {
    this.leftDendro = leftDendro;
    this.rightDendro = rightDendro;
    this.label = label;
  }

  public void printTree(final PrintStream ps) {
    
	if (leftDendro != null) {
		leftDendro.printTree(ps);
    }
    ps.print(label);
    if (rightDendro != null) {
    	rightDendro.printTree(ps);
    }
  }


  


}