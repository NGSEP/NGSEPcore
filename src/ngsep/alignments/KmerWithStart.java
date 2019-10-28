package ngsep.alignments;

class KmerWithStart {
	private CharSequence kmer;
	private int start;
	public KmerWithStart(CharSequence kmer, int start) {
		super();
		this.kmer = kmer;
		this.start = start;
	}
	/**
	 * @return the kmer
	 */
	public CharSequence getKmer() {
		return kmer;
	}
	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}	
}