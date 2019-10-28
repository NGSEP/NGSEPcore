package ngsep.alignments;

public class PairEndsAlignments{
	private ReadAlignment aln1;
	public ReadAlignment getAln1() {
		return aln1;
	}
	public ReadAlignment getAln2() {
		return aln2;
	}
	private ReadAlignment aln2;
	public PairEndsAlignments(ReadAlignment pAln1, ReadAlignment pAln2) {
		aln1=pAln1;
		aln2=pAln2;
	}

}