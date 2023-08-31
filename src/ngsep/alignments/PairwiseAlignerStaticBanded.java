/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.alignments;

import java.util.HashMap;
import java.util.Objects;

import ngsep.sequences.LimitedSequence; 

/**
 * 
 * @author Ana Sofia Castellanos
 *
 */
public class PairwiseAlignerStaticBanded implements PairwiseAligner {

    private int match =1; 

    private int mismatch = 1; 

    private int  indel = 2; 

    private HashMap<Tuple,Integer> nodes = new HashMap<>(); 

    private int k = 3; 
    
    
   

    public int getK() {
		return k;
	}

	public void setK(int k) {
		this.k = k;
	}

	@Override
    public String[] calculateAlignment(CharSequence sequence1, CharSequence sequence2 ) {
        // sequence1.length()> sequence2.length()-> Right definition of 

        if (checkminK(sequence1, sequence2)){
            calculateHashMap(sequence1, sequence2);
            return aligenedSequences(sequence1, sequence2); 
        }else{
            throw new RuntimeException("K value is not possible"); 
        }
    }
	
	public int getMaxScore(CharSequence sequence1, CharSequence sequence2 ) {
		
			return nodes.get(new Tuple(sequence1.length(),sequence2.length())); 

	}
	
    private void calculateHashMap(CharSequence sequence1, CharSequence sequence2){
    	 int i = sequence1.length();
         int j = sequence2.length();
         

        //Set node 0,0 
        nodes.put(new Tuple(0,0), 0); 
        //Set other nodes
        for (int row = 0; row <sequence1.length()+1; row++ ){
            for (int col = Math.max(0,row-this.k) ; col<Math.min(k+row+1, sequence2.length()+1); col++){
            	if (row == 0 && col == 0 ){
                 
                }else if (row == 0 && col <=this.k ){
                	
                    nodes.put(new Tuple(row,col), nodes.get(new Tuple(row, col-1))+ insertionCost(sequence1.charAt(row), sequence2.charAt(col-1))); 

                }else if (col == 0 && row != 0 ){
                	
                    nodes.put(new Tuple(row,col), nodes.get(new Tuple (row-1, col))+deletionCost(sequence1.charAt(row-1), sequence2.charAt(col))); 

                }else if(col-row ==this.k){
       
                	
                	int maxNum = Math.max(nodes.get(new Tuple(row-1,col-1))+matchMismatchCost(sequence1.charAt(row-1), sequence2.charAt(col-1)), 
      					  				  nodes.get(new Tuple(row, col-1))+insertionCost(sequence1.charAt(row-1), sequence2.charAt(col-1)));
                	nodes.put(new Tuple(row,col),maxNum); 
          
                	
                	
                }else if(row-col==this.k){
                	
                	int maxNum = Math.max(nodes.get(new Tuple(row-1,col-1))+matchMismatchCost(sequence1.charAt(row-1), sequence2.charAt(col-1)), 
                						  nodes.get(new Tuple (row-1, col))+deletionCost(sequence1.charAt(row-1), sequence2.charAt(col-1)));
                	nodes.put(new Tuple(row,col),maxNum); 
                	
                }else{
         
                	
                    int maxNum = Math.max(nodes.get(new Tuple(row-1,col-1))+matchMismatchCost(sequence1.charAt(row-1), sequence2.charAt(col-1)), 
                    					  nodes.get(new Tuple(row, col-1))+insertionCost(sequence1.charAt(row-1), sequence2.charAt(col-1))); 
                    maxNum = Math.max(maxNum, nodes.get(new Tuple (row-1, col))+deletionCost(sequence1.charAt(row-1), sequence2.charAt(col-1))); 
                    nodes.put(new Tuple(row,col),maxNum); 

                }

            }
        }

    }

    private String[] aligenedSequences(CharSequence sequence1, CharSequence sequence2){

        int i = sequence1.length();
        int j = sequence2.length();
        


        StringBuffer ns1 = new StringBuffer();
		StringBuffer ns2 = new StringBuffer();

        while (j>0 || i>0){
            int actual = nodes.get(new Tuple(i,j)); 
            if (i==0 && j!=0){
                ns1.append(LimitedSequence.GAP_CHARACTER); 
                ns2.append(sequence2.charAt(j-1)); 
                j--;
            }else if(j==0 && i!=0 ){
                ns1.append(sequence1.charAt(i-1)); 
                ns2.append(LimitedSequence.GAP_CHARACTER); 
                i--; 
            }else if (actual == nodes.get(new Tuple(i-1,j-1))+matchMismatchCost(sequence1.charAt(i-1), sequence2.charAt(j-1))){
                
            	
            	ns1.append(sequence1.charAt(i-1)); 
                ns2.append(sequence2.charAt(j-1)); 
                j--; 
                i--; 
            }else if ((j-i == this.k)&&(actual == nodes.get(new Tuple(i, j-1))+insertionCost(sequence1.charAt(i-1), sequence2.charAt(j-1)))){
                ns1.append(LimitedSequence.GAP_CHARACTER); 
                ns2.append(sequence2.charAt(j-1)); 
                j--;
            }else if ((i-j == this.k ) && (actual == nodes.get(new Tuple (i-1, j))+deletionCost(sequence1.charAt(i-1), sequence2.charAt(j-1)))){
                ns1.append(sequence1.charAt(i-1)); 
                ns2.append(LimitedSequence.GAP_CHARACTER); 
                i--; 
            }
            else if ((j-i != this.k) &&(i-j != this.k )) {
            	if(actual == nodes.get(new Tuple(i, j-1))+insertionCost(sequence1.charAt(i-1), sequence2.charAt(j-1))) {
            		ns1.append(LimitedSequence.GAP_CHARACTER); 
                    ns2.append(sequence2.charAt(j-1)); 
                    j--;
            	}else if (actual == nodes.get(new Tuple (i-1, j))+deletionCost(sequence1.charAt(i-1), sequence2.charAt(j-1))) {
            		ns1.append(sequence1.charAt(i-1)); 
                    ns2.append(LimitedSequence.GAP_CHARACTER); 
                    i--; 
            		
            	}else throw new RuntimeException("Unexpected score error at "+i+" "+j); 

            	
            }else throw new RuntimeException("Unexpected score error at "+i+" "+j); 
            }
        String[] seqs = new String[2]; 
        seqs[0] = ns1.reverse().toString();
        seqs[1] = ns2.reverse().toString();
        return seqs;

        



    }
    

    private int insertionCost(char l1, char l2){
        return -indel; 

    }
    private int deletionCost(char l1, char l2){
        return -indel; 

    }
    private int matchMismatchCost(char l1, char l2){
        if (l1 == l2){
            return match;

        }else{
            return -mismatch; 
        }

    }



    private boolean checkminK(CharSequence sequence1, CharSequence sequence2){
        return this.k!=0 &&((this.k>= sequence1.length()-sequence2.length())||(this.k>= sequence2.length()-sequence1.length()));
    }
    
}

class Tuple{
    private int num1; 
    private int num2; 

    public Tuple(int num1, int num2){
        this.num1 = num1; 
        this.num2 = num2; 
    }

    public int getNum1(){
        return num1; 
    }
    public int getNum2(){
        return num2; 
    }

    @Override
    public boolean equals(Object o) {
        if (this == o)
            return true;
        if (o == null || getClass() != o.getClass())
            return false;
        Tuple other = (Tuple) o ; 
        return num1 == other.num1 && num2 == other.num2; 

    }
    @Override
    public int hashCode(){
        return Objects.hash(num1,num2); 
    }

	@Override
	public String toString() {
		return "Tuple [num1=" + num1 + ", num2=" + num2 + "]";
	}


}
