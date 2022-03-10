package ngsep.genome;

import java.util.List;

public interface PairwiseSyntenyBlocksFinder {
	public static final int DEF_MIN_BLOCK_LENGTH = 100000;
	public static final int DEF_MIN_HOMOLOGY_UNITS_BLOCK = 6;
	public static final int DEF_MAX_DISTANCE_BETWEEN_UNITS = 100000;
	/**
	 * Calculates synteny blocks based on orthogroups defined by homology clusters
	 * @param g1 First genome to compare
	 * @param g2 Second genome to compare
	 * @param clusters Homology clusters
	 * @return List<PairwiseSyntenyBlock> List of pairwise synteny homolog clusters between the two groups 
	 */
	public List<PairwiseSyntenyBlock> findSyntenyBlocks(AnnotatedReferenceGenome g1, AnnotatedReferenceGenome g2, List<HomologyCluster> clusters);
}
