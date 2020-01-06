package ngsep.sequences;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class KmerHitsCluster {
	private CharSequence query;
	private List<FMIndexUngappedSearchHit> hits=new ArrayList<>();
	private int sequenceIdx;
	private String sequenceName;
	private int first;
	private int last;
	private Set<Integer> kmerNumbers = new HashSet<>();
	private boolean allConsistent = true;
	private boolean lastKmerPresent = false;

	public KmerHitsCluster(CharSequence query, FMIndexUngappedSearchHit kmerHit) {
		this.query = query;
		sequenceIdx = kmerHit.getSequenceIdx();
		sequenceName = kmerHit.getSequenceName();
		int kmerQueryStart = kmerHit.getQueryIdx();
		first = kmerHit.getStart() - kmerQueryStart +1;
		last = kmerHit.getStart()+(query.length()-kmerQueryStart);
		hits.add(kmerHit);
		kmerNumbers.add(kmerQueryStart);
		lastKmerPresent = kmerQueryStart+kmerHit.getQuery().length()==query.length();
	}

	public String getSequenceName() {
		return sequenceName;
	}
	

	/**
	 * @return the sequenceIdx
	 */
	public int getSequenceIdx() {
		return sequenceIdx;
	}

	public int getFirst() {
		return first;
	}

	public int getLast() {
		return last;
	}

	public int length() {
		return last-first+1;
	}

	public boolean addKmerHit(FMIndexUngappedSearchHit kmerHit, int toleranceChange) {
		int kmerQueryStart = kmerHit.getQueryIdx();
		int estFirst = kmerHit.getStart() - kmerQueryStart +1;
		int estLast = kmerHit.getStart()+(query.length()-kmerQueryStart);
		//System.out.println("Previous coords: "+first+"-"+last+" next cords: "+estFirst+"-"+estLast);
		if(first > estLast || last < estFirst) return false;
		if(toleranceChange>0 && Math.abs(first-estFirst)>toleranceChange) return false;
		if(toleranceChange>0 && Math.abs(last-estLast)>toleranceChange) return false;
		if(first != estFirst) allConsistent = false;
		if(last != estLast) allConsistent = false;
		
		if(!kmerNumbers.contains(kmerQueryStart)) kmerNumbers.add(kmerQueryStart);
		if(kmerQueryStart+kmerQueryStart+kmerHit.getQuery().length()==query.length()) lastKmerPresent=true;
		if(first>estFirst) first = estFirst;
		if(last<estLast) last = estLast;
		hits.add(kmerHit);
		return true;	
	}

	/**
	 * @return the query
	 */
	public CharSequence getQuery() {
		return query;
	}

	/**
	 * 
	 * @return the kmerNumbers
	 */
	public int getNumDifferentKmers() {
		return kmerNumbers.size();
	}

	/**
	 * @return the allConsistent
	 */
	public boolean isAllConsistent() {
		return allConsistent;
	}

	public boolean isFirstKmerPresent() {
		return kmerNumbers.contains(0);
	}
	/**
	 * @return the lastAlnPresent
	 */
	public boolean isLastKmerPresent() {
		return lastKmerPresent;
	}
	
	public List<FMIndexUngappedSearchHit> getHitsByQueryIdx () {
		List<FMIndexUngappedSearchHit> sortedHits = new ArrayList<FMIndexUngappedSearchHit>(hits.size());
		Collections.sort(sortedHits, new Comparator<FMIndexUngappedSearchHit>() {
			@Override
			public int compare(FMIndexUngappedSearchHit hit0, FMIndexUngappedSearchHit hit1) {
				return hit0.getQueryIdx()-hit1.getQueryIdx();
			}
		});
		return sortedHits;
	}
}