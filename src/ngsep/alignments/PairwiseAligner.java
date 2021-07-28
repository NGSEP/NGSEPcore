package ngsep.alignments;

public interface PairwiseAligner {
	public String [] calculateAlignment (CharSequence sequence1, CharSequence sequence2);
}
