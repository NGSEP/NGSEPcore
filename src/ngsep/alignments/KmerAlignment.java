package ngsep.alignments;

/**
 * Datastructure to asociate the alignment with the correponding kmer 
 * @author German Andrade
 *
 */
public class KmerAlignment {
	
	private int kmerNumber;
	
	private ReadAlignment alignment;

	public KmerAlignment(int i, ReadAlignment aln) {
		kmerNumber=i;
		alignment=aln;
	}

	public int getKmerNumber() {
		// TODO Auto-generated method stub
		return kmerNumber;
	}

	public ReadAlignment getReadAlignment() {
		// TODO Auto-generated method stub
		return alignment;
	}

}
