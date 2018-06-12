package ngsep.assembly;

import ngsep.alignments.ReadAlignment;

public class ReadOverlap2 {
    private final static int ERRORRATE = 5;

    private final int first1;
    private final int first2;
    private int last1;
    private int last2;
    private int errorCount;
    private boolean added;

    public ReadOverlap2(int first1, int first2, int last2) {
	super();
	this.first1 = first1;
	this.first2 = first2;
	this.last2 = last2;
	this.last1 = first1 + OverlapGraph.SEARCH_KMER_LENGTH;
	errorCount = 0;
	added = true;
    }

    public boolean addAlignment(ReadAlignment aln, int kmerSize) {
	if (added())
	    return false;
	int dif = Difenrece(last2 - OverlapGraph.SEARCH_KMER_OVERLAP, aln.getFirst());
	if (dif < ERRORRATE) {
	    added = true;
	    last2 = aln.getLast();
	    last1 += OverlapGraph.SEARCH_KMER_DISTANCE * errorCount;
	    last1 += kmerSize - OverlapGraph.SEARCH_KMER_OVERLAP;
	    errorCount = 0;
	    return true;
	} else
	    return false;
    }

    public boolean added() {
	return added;
    }

    private final static int Difenrece(int a, int b) {
	return Math.abs(a - b);
    }

    /**
     * @return the errorCount
     */
    public int getErrorCount() {
	return errorCount;
    }

    /**
     * @param errorCount
     *            the errorCount to set
     */
    public void errorCountPlusPlus() {
	this.errorCount++;
    }

    /**
     * @param diference
     *            the diference to set
     */
    public void nonAdded() {
	this.added = false;
    }

    /** @return the first1 */
    public int getFirst1() {
	return first1;
    }

    /** @return the first2 */
    public int getFirst2() {
	return first2;
    }

    /**
     * @return the last1
     */
    public int getLast1() {
	return last1;
    }

    /**
     * @return the last2
     */
    public int getLast2() {
	return last2;
    }

    public int length() {
	return last1 - first1;
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Object#toString()
     */
    @Override
    public String toString() {
	return "ReadOverlap2 [first1=" + first1 + ", first2=" + first2 + ", last1=" + last1 + ", last2=" + last2 + "]";
    }

}
