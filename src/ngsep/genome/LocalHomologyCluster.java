package ngsep.genome;

import java.util.ArrayList;
import java.util.List;

public class LocalHomologyCluster implements GenomicRegion {
	private List<HomologyUnit> units = new ArrayList<HomologyUnit>(); 
	private int genomeId;
	private HomologyCluster parent;
	private String seqName;
	private int first;
	private int last;
	public LocalHomologyCluster (HomologyCluster parent, HomologyUnit unit) {
		this.parent = parent;
		genomeId = unit.getGenomeId();
		units.add(unit);
		seqName = unit.getSequenceName();
		first = unit.getFirst();
		last = unit.getLast();
	}
	public void addUnit(HomologyUnit unit) {
		if(genomeId!=unit.getGenomeId()) throw new RuntimeException("Error merging local homology unit. Current genome id: "+genomeId+" new genome id: "+unit.getGenomeId());
		if(seqName!=unit.getSequenceName()) throw new RuntimeException("Error merging local homology unit. Current sequence name: "+seqName+" new sequence name: "+unit.getSequenceName());
		units.add(unit);
		first = Math.min(first, unit.getFirst());
		last = Math.max(last, unit.getLast());
	}
	
	
	public List<HomologyUnit> getUnits() {
		return units;
	}
	public int getGenomeId() {
		return genomeId;
	}
	public List<HomologyUnit> getHomologyUnitsCluster() {
		return units;
	}
	
	public HomologyCluster getParent() {
		return parent;
	}
	public void setParent(HomologyCluster parent) {
		this.parent = parent;
	}
	@Override
	public String getSequenceName() {
		return seqName;
	}
	@Override
	public int getFirst() {
		return first;
	}
	@Override
	public int getLast() {
		return last;
	}
	@Override
	public int length() {
		return last-first+1;
	}
	@Override
	public boolean isPositiveStrand() {
		return true;
	}
	@Override
	public boolean isNegativeStrand() {
		return false;
	}
}
