
package ngsep.haplotyping;

/**
 * Copied from SingleIndividualHaplotyper - Efficient heuristic algorithms for the SIH problem
 * @author Jorge Duitama
 */
public interface SIHAlgorithm {
	/**
	 * Builds a haplotype for the given matrix representing a block of overlapping fragments
	 * @param block that represents the matrix of overlapping fragments 
	 */
	public void buildHaplotype (HaplotypeBlock block);
}
