package ngsep.genome;

import java.util.List;

public class SyntenyBlock {

	private GenomicRegion regionGenome1;
	private GenomicRegion regionGenome2;
	private List<SyntenyEdge> homologies;
	
	
	
	public SyntenyBlock(List<SyntenyEdge> homologies) {
		SyntenyEdge firstEdge = homologies.get(0);
		int first1 = Integer.MAX_VALUE;
		int last1 = 0;
		int first2 = Integer.MAX_VALUE;
		int last2 = 0;
		for (SyntenyEdge se : homologies) {
			HomologyUnit source1 = se.getSource().getQueryUnit();
			HomologyUnit source2 = se.getSource().getSubjectUnit();
			HomologyUnit target1 = se.getTarget().getQueryUnit();
			HomologyUnit target2 = se.getTarget().getSubjectUnit();
			first1 = Math.min(first1, source1.getFirst());
			first1 = Math.min(first1, target1.getFirst());
			first2 = Math.min(first2, source2.getFirst());
			first2 = Math.min(first2, target2.getFirst());
			last1 = Math.max(last1, source1.getLast());
			last1 = Math.max(last1, target1.getLast());
			last2 = Math.max(last2, source2.getLast());
			last2 = Math.max(last2, target2.getLast());	
		}
		regionGenome1 = new GenomicRegionImpl(firstEdge.getSource().getQueryUnit().getSequenceName(), first1, last1);
		regionGenome2 = new GenomicRegionImpl(firstEdge.getSource().getSubjectUnit().getSequenceName(), first2, last2);
		this.homologies = homologies;
	}

	public List<SyntenyEdge> getHomologies() {
		return homologies;
	}

	public GenomicRegion getRegionGenome1() {
		return regionGenome1;
	}

	public GenomicRegion getRegionGenome2() {
		return regionGenome2;
	}
}
