package ngsep.transposons;

import java.util.List;

import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequence;

public interface DeNovoTransposableElementsFinder {

	/**
	 * Finds transposable elements in a refrence genome
	 * @param genome to analyze
	 * @return List<QualifiedSequence> List of sequences identified as transposable elements
	 */
	public List<QualifiedSequence> findTransposons(ReferenceGenome genome);
}
