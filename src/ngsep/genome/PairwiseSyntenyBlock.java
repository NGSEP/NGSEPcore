package ngsep.genome;

import java.util.List;

public class PairwiseSyntenyBlock {

	private GenomicRegion regionGenome1;
	private GenomicRegion regionGenome2;
	private List<SyntenyVertex> homologies;
	
	
	
	public PairwiseSyntenyBlock(int genomeId1, int genomeId2, List<SyntenyVertex> homologies) {
		String seqName1 = null;
		int first1 = Integer.MAX_VALUE;
		int last1 = 0;
		String seqName2 = null;
		int first2 = Integer.MAX_VALUE;
		int last2 = 0;
		for (SyntenyVertex sv : homologies) {
			HomologyUnit unit1 = sv.getHomologyCluster().findHomologyUnit(genomeId1);
			HomologyUnit unit2 = sv.getHomologyCluster().findHomologyUnit(genomeId2);
			if(seqName1==null) seqName1 = unit1.getSequenceName();
			first1 = Math.min(first1, unit1.getFirst());
			first2 = Math.min(first2, unit2.getFirst());
			if(seqName2==null) seqName2 = unit2.getSequenceName();
			last1 = Math.max(last1, unit1.getLast());
			last2 = Math.max(last2, unit2.getLast());	
		}
		regionGenome1 = new GenomicRegionImpl(seqName1, first1, last1);
		regionGenome2 = new GenomicRegionImpl(seqName2, first2, last2);
		this.homologies = homologies;
	}

	public List<SyntenyVertex> getHomologies() {
		return homologies;
	}

	public GenomicRegion getRegionGenome1() {
		return regionGenome1;
	}

	public GenomicRegion getRegionGenome2() {
		return regionGenome2;
	}
}
