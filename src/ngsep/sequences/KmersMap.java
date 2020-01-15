package ngsep.sequences;

import java.io.PrintStream;

import ngsep.math.Distribution;

public interface KmersMap {
	/**
	 * @return int number of kmers in the map
	 */
	public int size ();
	/**
	 * Gives the count of ocurrances of a given k-mer
	 * @param kmer To search
	 * @return int number of times it appears
	 */
	public int getCount(CharSequence kmer);
	/**
	 * Add 1 to the ocurrances of the given k-mer
	 * @param kmer to modify count
	 */
	public void addOcurrance(CharSequence kmer);
	/**
	 * Sets the given count for the given k-mer
	 * @param kmer to set count
	 * @param count new k-mer count
	 */
	public void setCount(CharSequence kmer, int count);
	/**
	 * Filters the k-mers map leaving only k-mers with at least the given abundance
	 * @param minAbundance Minimum abundance to keep the k-mer
	 */
	public void filterKmers(int minAbundance);
	/**
	 * Calculates the distribution of abundances of each kmer
	 * @return Distribution
	 */
	public Distribution calculateAbundancesDistribution();
	/**
	 * Saves the kmers in this map in the given output
	 * @param out Stream to save the k-mers
	 */
	public void save(PrintStream out);
}
