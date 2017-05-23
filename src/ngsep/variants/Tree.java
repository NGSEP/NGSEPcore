package ngsep.variants;

import java.io.PrintStream;



public class Tree {

  private Tree leftTree;
  private Tree rightTree;
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
  
  Tree(Tree leftTree, Tree rightTree,  String label) {
    this.leftTree = leftTree;
    this.rightTree = rightTree;
    this.label = label;
  }

  public void printTree(final PrintStream ps) {
    
	if (leftTree != null) {
      leftTree.printTree(ps);
    }
    ps.print(label);
    if (rightTree != null) {
      rightTree.printTree(ps);
    }
  }


  


}