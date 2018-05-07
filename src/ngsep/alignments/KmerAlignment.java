package ngsep.alignments;

import ngsep.genome.GenomicRegion;

/**
 * Data structure to associate the alignment with the corresponding kmer 
 * @author German Andrade
 *
 */
public class KmerAlignment implements GenomicRegion {
	
	private int kmerNumber;
	
	private ReadAlignment alignment;

	public KmerAlignment(int i, ReadAlignment aln) {
		kmerNumber=i;
		alignment=aln;
	}

	public int getKmerNumber() {
		return kmerNumber;
	}

	public ReadAlignment getReadAlignment() {
		return alignment;
	}

	@Override
	public String getSequenceName() {
		return alignment.getSequenceName();
	}

	@Override
	public int getFirst() {
		return alignment.getFirst();
	}

	@Override
	public int getLast() {
		return alignment.getLast();
	}

	@Override
	public int length() {
		return alignment.length();
	}

	@Override
	public boolean isPositiveStrand() {
		return alignment.isPositiveStrand();
	}

	@Override
	public boolean isNegativeStrand() {
		return alignment.isNegativeStrand();
	}
}
