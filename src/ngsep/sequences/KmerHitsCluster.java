package ngsep.sequences;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.genome.GenomicRegion;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class KmerHitsCluster implements GenomicRegion, Serializable {

	private static final long serialVersionUID = -4473724263138366546L;
	private CharSequence query;
	//Map indexed by query start
	private Map<Integer,FMIndexUngappedSearchHit> hitsMap=new TreeMap<Integer, FMIndexUngappedSearchHit>();
	private int sequenceIdx;
	private String sequenceName;
	private int first;
	private int last;
	private int queryStart;
	private int queryEnd;
	private int numDifferentKmers = 0;
	private double averageHitsQuery;
	private int totalKmersQuery;
	private double weightedCount=0;
	private boolean allConsistent = true;
	private boolean firstKmerPresent = false;
	private boolean lastKmerPresent = false;
	private static int idxDebug = -1;
	
	public KmerHitsCluster(CharSequence query, List<FMIndexUngappedSearchHit> inputHits) {
		this.query = query;
		if(inputHits.size()==0) throw new RuntimeException("Invalid empty input hits for cluster"+inputHits.size());
		
		FMIndexUngappedSearchHit firstHit = inputHits.get(0);
		sequenceIdx = firstHit.getSequenceIdx();
		sequenceName = firstHit.getSequenceName();
		if(sequenceIdx==idxDebug) System.out.println("Clustering "+inputHits.size()+" hits. Target idx: "+sequenceIdx);
		Map<Integer,List<FMIndexUngappedSearchHit>> hitsMultiMap = new TreeMap<Integer, List<FMIndexUngappedSearchHit>>();
		
		for(FMIndexUngappedSearchHit hit:inputHits) {
			List<FMIndexUngappedSearchHit> list = hitsMultiMap.computeIfAbsent(hit.getQueryIdx(), l -> new ArrayList<FMIndexUngappedSearchHit>());
			list.add(hit);
		}
		if(sequenceIdx==idxDebug) System.out.println("Num different kmers: "+hitsMultiMap.size());
		List<Integer> subjectStarts = new ArrayList<Integer>();
		//Try first with local unique hits
		double sum = 0;
		double sum2 = 0;
		double n = 0;
		for(List<FMIndexUngappedSearchHit> hits:hitsMultiMap.values()) {
			if(hits.size()>1) continue;
			FMIndexUngappedSearchHit hit = hits.get(0);
			int estFirst = estimateSubjectFirst(hit);
			subjectStarts.add(estFirst);
			sum+=estFirst;
			sum2+=estFirst*estFirst;
			n++;
		}
		if(sequenceIdx==idxDebug) System.out.println("Num unique: "+n);
		if(n<5) {
			subjectStarts.clear();
			sum=sum2=n=0;
			for(List<FMIndexUngappedSearchHit> hits:hitsMultiMap.values()) {
				for(FMIndexUngappedSearchHit hit:hits) {
					int estFirst = estimateSubjectFirst(hit);
					subjectStarts.add(estFirst);
					sum+=estFirst;
					sum2+=estFirst*estFirst;
					n++;
				}
				
			}
		}
		Collections.sort(subjectStarts);
		int median = subjectStarts.get(subjectStarts.size()/2);
		double stdev = Math.sqrt((sum2-sum*sum/n)/(n-1));
		if(sequenceIdx==idxDebug) System.out.println("Num unique: "+n+" median: "+median+" stdev: "+stdev);
		int maxDistance = (int)stdev;
		if(maxDistance<5) maxDistance*=2;
		if(stdev<0.1*query.length()) maxDistance*=2;
		else if (stdev>0.2*query.length()) maxDistance/=2;
		
		first = -1;
		for(List<FMIndexUngappedSearchHit> hits:hitsMultiMap.values()) {
			FMIndexUngappedSearchHit hit = selectHit(hits,median, Math.min(query.length()/20, maxDistance));
			if(hit!=null) {
				if(first==-1) {
					first = estimateSubjectFirst(hit);
					last = estimateSubjectLast(hit);
					queryStart = hit.getQueryIdx();
					queryEnd = hit.getQueryIdx() + hit.getQuery().length();
				}
				addHit(hit);
			}
		}
		if(sequenceIdx==idxDebug) System.out.println("Final hits: "+hitsMap.size()+" first: "+first+" last: "+last);
	}

	private FMIndexUngappedSearchHit selectHit(List<FMIndexUngappedSearchHit> hits, int median, int maxDistance) {
		FMIndexUngappedSearchHit answer = null;
		int minDistance = 0;
		for(FMIndexUngappedSearchHit hit:hits) {
			int estFirst = estimateSubjectFirst(hit);
			int distance = Math.abs(estFirst-median);
			if(sequenceIdx==idxDebug) System.out.println("Next hit: "+hit.getQueryIdx()+" first: "+estFirst+" distance: "+distance+ " max: "+maxDistance);
			if(distance <= maxDistance) {
				if(answer == null || minDistance>distance) {
					answer = hit;
					minDistance = distance;
				}
			}
		}
		return answer;
	}



	public KmerHitsCluster(CharSequence query, FMIndexUngappedSearchHit kmerHit) {
		this.query = query;
		sequenceIdx = kmerHit.getSequenceIdx();
		sequenceName = kmerHit.getSequenceName();
		int kmerQueryStart = kmerHit.getQueryIdx();
		first = estimateSubjectFirst(kmerHit);
		last = estimateSubjectLast(kmerHit);
		queryStart = kmerHit.getQueryIdx();
		queryEnd = kmerHit.getQueryIdx() + kmerHit.getQuery().length();
		hitsMap.put(kmerQueryStart, kmerHit);
		numDifferentKmers = 1;
		firstKmerPresent = queryStart == 0;
		lastKmerPresent = queryEnd==query.length();
	}
	
	private void addHit(FMIndexUngappedSearchHit hit) {
		hitsMap.put(hit.getQueryIdx(), hit);
		int estFirst = estimateSubjectFirst(hit);
		int estLast = estimateSubjectLast(hit);
		if(estFirst!=first || estLast!=last) allConsistent = false;
		first = Math.min(first, estFirst);
		last = Math.max(last, estLast);
		queryStart = Math.min(queryStart, hit.getQueryIdx());
		queryEnd = Math.max(queryEnd, hit.getQueryIdx() + hit.getQuery().length());
		if(queryStart==0) firstKmerPresent = true;
		if (queryEnd==query.length()) lastKmerPresent = true;
		numDifferentKmers++;
	}
	
	public boolean addKmerHit(FMIndexUngappedSearchHit kmerHit, int toleranceChange) {
		int estFirst = estimateSubjectFirst(kmerHit);
		int estLast = estimateSubjectLast(kmerHit);
		//System.out.println("Hit with idx: "+kmerHit.getQueryIdx()+" Previous coords: "+first+"-"+last+" next cords: "+estFirst+"-"+estLast);
		if(first > estLast || last < estFirst) return false;
		if(toleranceChange>0 && Math.abs(first-estFirst)>toleranceChange) return false;
		if(toleranceChange>0 && Math.abs(last-estLast)>toleranceChange) return false;
		addHit(kmerHit);
		return true;	
	}
	
	private int estimateSubjectFirst(FMIndexUngappedSearchHit hit) {
		return hit.getStart() - hit.getQueryIdx() +1;
	}

	private int estimateSubjectLast(FMIndexUngappedSearchHit hit) {
		return hit.getStart()+(query.length()-hit.getQueryIdx());
	}
	
	public void summarize(double averageHitsQuery, int kmersQuery) {
		this.averageHitsQuery = averageHitsQuery;
		this.totalKmersQuery = kmersQuery;
		numDifferentKmers = hitsMap.size();
		weightedCount = 0;
		for(FMIndexUngappedSearchHit hit: hitsMap.values()) {
			double n = hit.getTotalHitsQuery(); 
			if(n<=averageHitsQuery) weightedCount++;
			else weightedCount += averageHitsQuery/n;
		}
		//Disposes detailed information about hits
		hitsMap.clear();
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

	/**
	 * @return the query
	 */
	public CharSequence getQuery() {
		return query;
	}
	
	public int getQueryStart() {
		return queryStart;
	}

	public int getQueryEnd() {
		return queryEnd;
	}
	
	public double getQueryCoverage () {
		double queryCoverage = queryEnd-queryStart;
		return queryCoverage / query.length();
	}
	
	public int getTotalKmersQuery() {
		return totalKmersQuery;
	}

	public double getAverageHitsQuery() {
		return averageHitsQuery;
	}
	public double getWeightedCount() {
		return weightedCount;
	}
	
	public double getProportionKmers () {
		return weightedCount / (double)totalKmersQuery;
	}

	/**
	 * 
	 * @return the kmerNumbers
	 */
	public int getNumDifferentKmers() {
		return numDifferentKmers;
	}

	/**
	 * @return the allConsistent
	 */
	public boolean isAllConsistent() {
		return allConsistent;
	}

	public boolean isFirstKmerPresent() {
		return firstKmerPresent;
	}
	/**
	 * @return the lastAlnPresent
	 */
	public boolean isLastKmerPresent() {
		return lastKmerPresent;
	}
	
	
	
	
	public List<FMIndexUngappedSearchHit> getHitsByQueryIdx () {
		List<FMIndexUngappedSearchHit> sortedHits = new ArrayList<FMIndexUngappedSearchHit>();
		sortedHits.addAll(hitsMap.values());
		return sortedHits;
	}

	public FMIndexUngappedSearchHit getKmerHit(int queryKmerIdx) {
		return hitsMap.get(queryKmerIdx);
	}

	@Override
	public boolean isPositiveStrand() {
		return true;
	}

	@Override
	public boolean isNegativeStrand() {
		// TODO Auto-generated method stub
		return false;
	}
	
	
}