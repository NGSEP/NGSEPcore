package ngsep.alignments;

import java.util.List;

import ngsep.sequences.RawRead;

public interface ReadAlignmentAlgorithm {
	/**
	 * Aligns the given read to a reference genome.
	 * It is assumed that the genome and/or a corresponding data structure is loaded in the implementation class
	 * @param read to be aligned
	 * @return List<ReadAlignment> Alignments found for the given read
	 */
	public List<ReadAlignment> alignRead (RawRead read);
}
