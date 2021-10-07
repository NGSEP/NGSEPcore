package ngsep.genome;

import java.util.List;

public class PairwiseSyntenyBlock {

	private int genomeId1;
	private GenomicRegion regionGenome1;
	private int genomeId2;
	private GenomicRegion regionGenome2;
	private List<SyntenyVertex> homologies;
	
	
	
	public PairwiseSyntenyBlock(List<SyntenyVertex> homologies) {
		String seqName1 = null;
		int first1 = Integer.MAX_VALUE;
		int last1 = 0;
		String seqName2 = null;
		int first2 = Integer.MAX_VALUE;
		int last2 = 0;
		for (SyntenyVertex sv : homologies) {
			LocalHomologyCluster g1 = sv.getLocalRegion1();
			LocalHomologyCluster g2 = sv.getLocalRegion2();
			if(seqName1==null) {
				genomeId1 = g1.getGenomeId();
				seqName1 = g1.getSequenceName();
			}
			first1 = Math.min(first1, g1.getFirst());
			first2 = Math.min(first2, g2.getFirst());
			if(seqName2==null) {
				genomeId2 = g2.getGenomeId();
				seqName2 = g2.getSequenceName();
			}
			last1 = Math.max(last1, g1.getLast());
			last2 = Math.max(last2, g2.getLast());	
		}
		regionGenome1 = new GenomicRegionImpl(seqName1, first1, last1);
		regionGenome2 = new GenomicRegionImpl(seqName2, first2, last2);
		this.homologies = homologies;
	}

	public List<SyntenyVertex> getHomologies() {
		return homologies;
	}

	
	public int getGenomeId1() {
		return genomeId1;
	}

	public int getGenomeId2() {
		return genomeId2;
	}

	public GenomicRegion getRegionGenome1() {
		return regionGenome1;
	}

	public GenomicRegion getRegionGenome2() {
		return regionGenome2;
	}
}
