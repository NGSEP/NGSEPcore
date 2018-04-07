package ngsep.alignments;

import java.util.List;

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
