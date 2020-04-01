package ngsep.sequences;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class KmerHitsCluster implements Serializable {

	private static final long serialVersionUID = -4473724263138366546L;
	private CharSequence query;
	//Map indexed by query start
	private Map<Integer,UngappedSearchHit> hitsMap=new TreeMap<Integer, UngappedSearchHit>();
	private int sequenceIdx;
	private String sequenceName;
	private int subjectPredictedStart;
	private int subjectPredictedEnd;
	private int queryStart;
	private int queryEnd;
	private int subjectEvidenceStart;
	private int subjectEvidenceEnd;
	private int numDifferentKmers = 0;
	private double averageHitsQuery;
	private double weightedCount=0;
	private boolean allConsistent = true;
	private boolean firstKmerPresent = false;
	private boolean lastKmerPresent = false;
	private static int idxDebug = -2;
	
	public KmerHitsCluster(CharSequence query, List<UngappedSearchHit> inputHits) {
		this.query = query;
		if(inputHits.size()==0) throw new RuntimeException("Invalid empty input hits for cluster"+inputHits.size());
		
		UngappedSearchHit firstHit = inputHits.get(0);
		sequenceIdx = firstHit.getSequenceIdx();
		sequenceName = firstHit.getSequenceName();
		if(sequenceIdx==idxDebug) System.out.println("KmerHitsCluster. Clustering "+inputHits.size()+" hits. Target idx: "+sequenceIdx);
		Map<Integer,List<UngappedSearchHit>> hitsMultiMap = new TreeMap<Integer, List<UngappedSearchHit>>();
		
		for(UngappedSearchHit hit:inputHits) {
			List<UngappedSearchHit> list = hitsMultiMap.computeIfAbsent(hit.getQueryIdx(), l -> new ArrayList<UngappedSearchHit>());
			list.add(hit);
		}
		if(sequenceIdx==idxDebug) System.out.println("KmerHitsCluster. Num different kmers: "+hitsMultiMap.size());
		List<Integer> subjectStarts = new ArrayList<Integer>();
		//Try first with local unique hits
		double sum = 0;
		double sum2 = 0;
		double n = 0;
		for(List<UngappedSearchHit> hits:hitsMultiMap.values()) {
			if(hits.size()>1) continue;
			UngappedSearchHit hit = hits.get(0);
			int estStart = estimateSubjectStart(hit);
			subjectStarts.add(estStart);
			sum+=estStart;
			sum2+=estStart*estStart;
			n++;
		}
		if(sequenceIdx==idxDebug) System.out.println("KmerHitsCluster. Num unique: "+n);
		if(n<5) {
			subjectStarts.clear();
			sum=sum2=n=0;
			for(List<UngappedSearchHit> hits:hitsMultiMap.values()) {
				for(UngappedSearchHit hit:hits) {
					int estStart = estimateSubjectStart(hit);
					subjectStarts.add(estStart);
					sum+=estStart;
					sum2+=(estStart*estStart);
					n++;
				}
				
			}
		}
		Collections.sort(subjectStarts);
		int median = subjectStarts.get(subjectStarts.size()/2);
		double variance = (sum2-sum*sum/n)/(n-1);
		double stdev = (variance>0)?Math.sqrt(variance):0;
		if(sequenceIdx==idxDebug) System.out.println("KmerHitsCluster. Num unique: "+n+" median: "+median+" variance: "+variance+" stdev: "+stdev);
		int maxDistance = (int)(stdev);
		if(maxDistance < 50) maxDistance=50;
		if(maxDistance<0.1*query.length()) maxDistance*=2;
		else if (maxDistance>0.2*query.length()) maxDistance/=2;
		
		subjectPredictedStart = -1;
		for(List<UngappedSearchHit> hits:hitsMultiMap.values()) {
			UngappedSearchHit hit = selectHit(hits,median, Math.min(query.length()/20, maxDistance));
			if(hit!=null) {
				if(subjectPredictedStart==-1) {
					subjectPredictedStart = estimateSubjectStart(hit);
					subjectPredictedEnd = estimateSubjectEnd(hit);
					queryStart = hit.getQueryIdx();
					queryEnd = hit.getQueryIdx() + hit.getQuery().length();
					subjectEvidenceStart = hit.getStart();
					subjectEvidenceEnd = subjectEvidenceStart+hit.getQuery().length();
				}
				addHit(hit);
			}
		}
		if(sequenceIdx==idxDebug) System.out.println("Final hits: "+hitsMap.size()+" start: "+subjectPredictedStart+" end: "+subjectPredictedEnd);
	}

	private UngappedSearchHit selectHit(List<UngappedSearchHit> hits, int median, int maxDistance) {
		UngappedSearchHit answer = null;
		int minDistance = 0;
		for(UngappedSearchHit hit:hits) {
			int estStart = estimateSubjectStart(hit);
			int distance = Math.abs(estStart-median);
			if(sequenceIdx==idxDebug) System.out.println("Next hit: "+hit.getQueryIdx()+" start: "+estStart+" distance: "+distance+ " max: "+maxDistance);
			if(distance <= maxDistance) {
				if(answer == null || minDistance>distance) {
					answer = hit;
					minDistance = distance;
				}
			}
		}
		return answer;
	}



	public KmerHitsCluster(CharSequence query, UngappedSearchHit kmerHit) {
		this.query = query;
		sequenceIdx = kmerHit.getSequenceIdx();
		sequenceName = kmerHit.getSequenceName();
		int kmerQueryStart = kmerHit.getQueryIdx();
		subjectPredictedStart = estimateSubjectStart(kmerHit);
		subjectPredictedEnd = estimateSubjectEnd(kmerHit);
		queryStart = kmerHit.getQueryIdx();
		queryEnd = kmerHit.getQueryIdx() + kmerHit.getQuery().length();
		subjectEvidenceStart = kmerHit.getStart();
		subjectEvidenceEnd = subjectEvidenceStart+kmerHit.getQuery().length();
		hitsMap.put(kmerQueryStart, kmerHit);
		numDifferentKmers = 1;
		firstKmerPresent = queryStart == 0;
		lastKmerPresent = queryEnd==query.length();
	}
	
	private void addHit(UngappedSearchHit hit) {
		hitsMap.put(hit.getQueryIdx(), hit);
		int estStart = estimateSubjectStart(hit);
		int estEnd = estimateSubjectEnd(hit);
		if(estStart!=subjectPredictedStart || estEnd!=subjectPredictedEnd) allConsistent = false;
		subjectPredictedStart = Math.min(subjectPredictedStart, estStart);
		subjectPredictedEnd = Math.max(subjectPredictedEnd, estEnd);
		queryStart = Math.min(queryStart, hit.getQueryIdx());
		queryEnd = Math.max(queryEnd, hit.getQueryIdx() + hit.getQuery().length());
		subjectEvidenceStart = Math.min(subjectEvidenceStart, hit.getStart());
		subjectEvidenceEnd = Math.max(subjectEvidenceEnd, hit.getStart()+hit.getQuery().length());
		if(queryStart==0) firstKmerPresent = true;
		if (queryEnd==query.length()) lastKmerPresent = true;
		numDifferentKmers++;
	}
	
	public boolean addKmerHit(UngappedSearchHit kmerHit, int toleranceChange) {
		int estStart = estimateSubjectStart(kmerHit);
		int estEnd = estimateSubjectEnd(kmerHit);
		//System.out.println("Hit with idx: "+kmerHit.getQueryIdx()+" Previous coords: "+first+"-"+last+" next cords: "+estFirst+"-"+estLast);
		if(subjectPredictedStart > estEnd || subjectPredictedEnd < estStart) return false;
		if(toleranceChange>0 && Math.abs(subjectPredictedStart-estStart)>toleranceChange) return false;
		if(toleranceChange>0 && Math.abs(subjectPredictedEnd-estEnd)>toleranceChange) return false;
		addHit(kmerHit);
		return true;
	}
	
	private int estimateSubjectStart(UngappedSearchHit hit) {
		return hit.getStart() - hit.getQueryIdx();
	}

	private int estimateSubjectEnd(UngappedSearchHit hit) {
		return hit.getStart()+(query.length()-hit.getQueryIdx());
	}
	
	public void summarize(double averageHitsQuery) {
		this.averageHitsQuery = averageHitsQuery;
		numDifferentKmers = hitsMap.size();
		weightedCount = 0;
		for(UngappedSearchHit hit: hitsMap.values()) {
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

	public int getSubjectPredictedStart() {
		return subjectPredictedStart;
	}

	public int getSubjectPredictedEnd() {
		return subjectPredictedEnd;
	}
	
	public void setSubjectPredictedLimits (int predictedStart, int predictedEnd) {
		if(predictedEnd<=predictedStart) throw new IllegalArgumentException("Predicted end "+predictedEnd+" must be larger than predicted start: "+predictedStart);
		subjectPredictedStart = predictedStart;
		subjectPredictedEnd = predictedEnd;
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

	public int getSubjectEvidenceStart() {
		return subjectEvidenceStart;
	}

	public int getSubjectEvidenceEnd() {
		return subjectEvidenceEnd;
	}

	public double getQueryCoverage () {
		double queryCoverage = queryEnd-queryStart;
		return queryCoverage / query.length();
	}

	public double getAverageHitsQuery() {
		return averageHitsQuery;
	}
	public double getWeightedCount() {
		return weightedCount;
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
	
	public List<UngappedSearchHit> getHitsByQueryIdx () {
		List<UngappedSearchHit> sortedHits = new ArrayList<UngappedSearchHit>();
		sortedHits.addAll(hitsMap.values());
		return sortedHits;
	}

	public UngappedSearchHit getKmerHit(int queryKmerIdx) {
		return hitsMap.get(queryKmerIdx);
	}
	
}