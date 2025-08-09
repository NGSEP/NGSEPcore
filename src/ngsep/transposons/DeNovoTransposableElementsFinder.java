package ngsep.transposons;

import java.util.List;

import ngsep.genome.ReferenceGenome;

public interface DeNovoTransposableElementsFinder {

	/**
	 * Finds transposable elements in a refrence genome
	 * @param genome to analyze
	 * @return List<TransposableElement> List of sequences identified as transposable elements
	 */
	public List<TransposableElementAnnotation> findTransposons(ReferenceGenome genome);
}
