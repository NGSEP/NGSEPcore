package ngsep.discovery;

import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.variants.GenomicVariant;

public class ActiveIndelRegionsFinder {
	private PriorityQueue<GenomicRegion> rawRegions = new PriorityQueue<GenomicRegion>((r1,r2)->r1.getFirst()-r2.getFirst());
	private GenomicRegionSortedCollection<GenomicRegion> activeRegions = new GenomicRegionSortedCollection<GenomicRegion>();
	
	private int minEventSupportActiveRegion = 2;
	private int flankBasepairsIndel = 1;
	
	
	private int lastActiveRegionEnd = 0;
	
	
	
	
	public int getMinEventSupportActiveRegion() {
		return minEventSupportActiveRegion;
	}

	public void setMinEventSupportActiveRegion(int minEventSupportActiveRegion) {
		this.minEventSupportActiveRegion = minEventSupportActiveRegion;
	}

	public int getFlankBasepairsIndel() {
		return flankBasepairsIndel;
	}

	public void setFlankBasepairsIndel(int flankBasepairsIndel) {
		this.flankBasepairsIndel = flankBasepairsIndel;
	}

	public void processAlignment(ReadAlignment aln) {
		Map<Integer,GenomicVariant> alnIndelCalls = aln.getIndelCalls();
		if(alnIndelCalls==null) return;
		for(GenomicVariant indelCall:alnIndelCalls.values()) {
			if(indelCall.length()>50) continue;
			if(activeRegions.findSpanningRegions(indelCall).size()>0) continue;
			GenomicRegion mononucleotideReg = aln.checkMononucleotide(indelCall.getFirst());
			if(mononucleotideReg!=null) rawRegions.add(mononucleotideReg);
			else if(indelCall.getLast()-indelCall.getFirst()>1) rawRegions.add(new GenomicRegionImpl(aln.getSequenceName(), indelCall.getFirst(), indelCall.getLast()));
			else rawRegions.add(new GenomicRegionImpl(aln.getSequenceName(), indelCall.getFirst()-flankBasepairsIndel, indelCall.getLast()+flankBasepairsIndel));
		}
		updateActiveRegions(aln.getFirst());
	}
	private int countSupport = 0;
	private GenomicRegionImpl nextRegionCandidate;
	private int lastPositionUpdate = 0;
	private void updateActiveRegions(int currentPos) {
		if(currentPos<=lastPositionUpdate) return;
		GenomicRegion rawRegion = rawRegions.peek();
		//if(currentPos>3000 && currentPos<4000) System.out.println("Updating active regions until: "+currentPos+" . Raw regions: "+rawRegions.size());
		while(rawRegion != null && rawRegion.getFirst()<currentPos) {
			rawRegions.remove();
			//if(currentPos>3000 && currentPos<4000) System.out.println("Updating active regions until: "+currentPos+" . Next region: "+rawRegion.getFirst()+"-"+rawRegion.getLast());
			if(nextRegionCandidate!=null && GenomicRegionSpanComparator.getInstance().span(nextRegionCandidate, rawRegion) ) {
				nextRegionCandidate.setLast(Math.max(nextRegionCandidate.getLast(), rawRegion.getLast()));
				countSupport++;
			} else if (nextRegionCandidate==null) {
				nextRegionCandidate = (GenomicRegionImpl)rawRegion;
			} else {
				//if(nextRegionCandidate.length()>20) System.out.println("Adding region "+nextRegionCandidate.getSequenceName()+":"+nextRegionCandidate.getFirst()+"-"+nextRegionCandidate.getLast()+" "+nextRegionCandidate.length()+" support: "+countSupport);
				addActiveRegion (nextRegionCandidate);		
				nextRegionCandidate = (GenomicRegionImpl)rawRegion;
				countSupport=1;
			}
			rawRegion = rawRegions.peek();
		}
		if(nextRegionCandidate!=null && nextRegionCandidate.getLast()<currentPos) {
			addActiveRegion(nextRegionCandidate);
			nextRegionCandidate = null;
		}
		lastPositionUpdate = currentPos;
	}
	public boolean isRegionInProgress() {
		return nextRegionCandidate!=null;
	}

	private void addActiveRegion(GenomicRegionImpl nextRegionCandidate) {
		if(countSupport>=minEventSupportActiveRegion) {
			activeRegions.add(nextRegionCandidate);
			lastActiveRegionEnd = Math.max(lastActiveRegionEnd, nextRegionCandidate.getLast());
		}
	}
	
	public void startPileups() {
		activeRegions.forceSort();
	}

	public void updatePileup(PileupRecord pileup) {
		int posPrint = -1;
		List<GenomicRegion> activePileup = activeRegions.findSpanningRegions(pileup.getSequenceName(),pileup.getPosition()).asList();
		if(pileup.getPosition()==posPrint) System.out.println("Spanning Active regions: "+activePileup.size()+" Total active regions: "+activeRegions.size());
		if(activePileup.size()==0) return;
		GenomicRegion r = activePileup.get(0);
		if(pileup.getPosition()==posPrint) System.out.println("Spanning region selected: "+r.getSequenceName()+":"+r.getFirst()+"-"+r.getLast());
		pileup.setActiveRegion(r);
	}

	public void setInputVariants(List<GenomicVariant> inputVariants) {
		for(GenomicRegion region:inputVariants) {
			if(region.getLast()-region.getFirst()>0) activeRegions.add(region);
		}
	}
}
