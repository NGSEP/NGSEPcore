package ngsep.alignments;

import java.util.HashMap;

import ngsep.sequences.LimitedSequence; 

public class PairwiseAlignerStaticBanded implements PairwiseAligner {

    private int match =1; 

    private int mismatch = -1; 

    private int  indel = 0; 

    private HashMap<Tuple,Integer> nodes = new HashMap<>(); 

    private int k = 0; 

    @Override
    public String[] calculateAlignment(CharSequence sequence1, CharSequence sequence2 ) {
        // sequence2.length()> sequence1.length()-> Right definition of 

        if (checkminK(sequence1, sequence2, k)){
            calculateHashMap(sequence1, sequence2, k);
            return aligenedSequences(sequence1, sequence2, k); 
        }else{
            throw new RuntimeException("K value is not possible"); 
        }

       


    }

    private void calculateHashMap(CharSequence sequence1, CharSequence sequence2, int k ){
        //Set node 0,0 
        nodes.put(new Tuple(0,0), 0); 
        //Set other nodes
        for (int row = 0; row <sequence1.length(); row++ ){
            for (int col = Math.max(0,row-k) ; col<Math.min(k+row+1, sequence2.length()); col++){
                Tuple actual = new Tuple(row,col); 
                int valDiag = nodes.get(new Tuple(row-1,col-1)); 
                int valIzq  = nodes.get(new Tuple(row, col-1));
                int valUp = nodes.get(new Tuple (row-1, col)); 
                if (row == 0 && col == 0 ){
                 
                }else if (row == 0 && col <=k ){
                    nodes.put(actual, valIzq+ insertionCost(sequence1.charAt(row), sequence2.charAt(col-1))); 

                }else if (col == 0 && row != 0 ){
                    nodes.put(actual, valUp+deletionCost(sequence1.charAt(row-1), sequence2.charAt(col))); 

                }else{
                    int maxNum = Math.max(valDiag+matchMismatchCost(sequence1.charAt(row-1), sequence2.charAt(col)), 
                                            valIzq+insertionCost(sequence1.charAt(row), sequence2.charAt(col-1))); 
                    maxNum = Math.max(maxNum, valUp+deletionCost(sequence1.charAt(row-1), sequence2.charAt(col))); 

                    nodes.put(actual,maxNum); 

                }

            }
        }
    }

    private String[] aligenedSequences(CharSequence sequence1, CharSequence sequence2, int k ){

        int i = sequence1.length();
        int j = sequence2.length();

        StringBuffer ns1 = new StringBuffer();
		StringBuffer ns2 = new StringBuffer();

        while (j>0 || i>0){
            int actual = nodes.get(new Tuple(i,j)); 
            int valDiag = nodes.get(new Tuple(i-1,j-1)); 
            int valIzq  = nodes.get(new Tuple(i, j-1));
            int valUp = nodes.get(new Tuple (i-1, j)); 
            if (i==0 && j!=0){
                ns1.append(LimitedSequence.GAP_CHARACTER); 
                ns2.append(sequence2.charAt(j-1)); 
                j--;
            }else if(j==0 && i!=0 ){
                ns1.append(sequence1.charAt(i-1)); 
                ns2.append(LimitedSequence.GAP_CHARACTER); 
                i--; 
            }else if (actual == valDiag+matchMismatchCost(sequence1.charAt(i-1), sequence2.charAt(j-1))){
                ns1.append(sequence1.charAt(i-1)); 
                ns2.append(sequence2.charAt(j-1)); 
                j--; 
                i--; 
            }else if ((j-1 == k)&&(actual == valIzq+insertionCost(sequence1.charAt(i), sequence2.charAt(j-1)) )){
                ns1.append(LimitedSequence.GAP_CHARACTER); 
                ns2.append(sequence2.charAt(j-1)); 
                j--;
            }else if ((i-j == k ) && (actual == valUp+deletionCost(sequence1.charAt(i), sequence2.charAt(j-1)))){
                ns1.append(sequence1.charAt(i-1)); 
                ns2.append(LimitedSequence.GAP_CHARACTER); 
                i--; 
            }
            else throw new RuntimeException("Unexpected score error at "+i+" "+j); 
            }
        String[] seqs = new String[2]; 
        seqs[0] = ns1.reverse().toString();
        seqs[1] = ns2.reverse().toString();
        return seqs;

        



    }
    

    private int insertionCost(char l1, char l2){
        return indel; 

    }
    private int deletionCost(char l1, char l2){
        return indel; 

    }
    private int matchMismatchCost(char l1, char l2){
        if (l1 == l2){
            return match;

        }else{
            return mismatch; 
        }

    }



    private boolean checkminK(CharSequence sequence1, CharSequence sequence2, int k ){
        return (k>= sequence1.length()-sequence2.length());
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

    public int hashCode(){
        return this.hashCode(); 
    }


}
