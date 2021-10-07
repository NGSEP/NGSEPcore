package ngsep.genome;

import java.util.List;

public interface PairwiseSyntenyBlocksFinder {
	public List<PairwiseSyntenyBlock> findSyntenyBlocks(AnnotatedReferenceGenome g1, AnnotatedReferenceGenome g2, List<HomologyCluster> clusters);
}
