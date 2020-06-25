package ngsep.sequences;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.math.Distribution;

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
	private int sequenceLength;
	private int subjectPredictedStart;
	private int subjectPredictedEnd;
	private int queryPredictedStart;
	private int queryPredictedEnd;
	private int queryEvidenceStart;
	private int queryEvidenceEnd;
	private int subjectEvidenceStart;
	private int subjectEvidenceEnd;
	private int predictedOverlap;
	private int numDifferentKmers = 0;
	private int selfHitsCountQuery = 0;
	private double averageHitsQuery;
	private double weightedCount=0;
	private boolean allConsistent = true;
	private boolean firstKmerPresent = false;
	private boolean lastKmerPresent = false;
	private static int idxSubjectDebug = -2;
	private static int queryLengthDebug = -1;
	
	public KmerHitsCluster(CharSequence query, List<UngappedSearchHit> inputHits) {
		this.query = query;
		if(inputHits.size()==0) throw new RuntimeException("Invalid empty input hits for cluster"+inputHits.size());
		
		UngappedSearchHit firstHit = inputHits.get(0);
		sequenceIdx = firstHit.getSequenceIdx();
		sequenceName = firstHit.getSequenceName();
		sequenceLength = firstHit.getSequenceLength();
		if(sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("KmerHitsCluster. Clustering "+inputHits.size()+" hits. Subject idx: "+sequenceIdx);
		Map<Integer,List<UngappedSearchHit>> hitsMultiMap = new TreeMap<Integer, List<UngappedSearchHit>>();
		
		for(UngappedSearchHit hit:inputHits) {
			//if (sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" estq: "+estimateQueryStart(hit)+" - "+estimateQueryEnd(hit)+" estS: "+estimateSubjectStart(hit)+" - "+estimateSubjectEnd(hit));
			List<UngappedSearchHit> list = hitsMultiMap.computeIfAbsent(hit.getQueryIdx(), l -> new ArrayList<UngappedSearchHit>());
			list.add(hit);
		}
		if(sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("KmerHitsCluster. Num different kmers: "+hitsMultiMap.size());
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
		if(sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("KmerHitsCluster. Num unique: "+n);
		if(n<5) {
			//System.out.println("WARN. Few unique kmers for hits to "+sequenceIdx+" initial: "+hitsMultiMap.size()+" final: "+n);
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
		Distribution dist = new Distribution(0, 100, 1);
		for(int start:subjectStarts) {
			int distance = Math.abs(start-median);
			if (distance < stdev) dist.processDatapoint(distance);
		}
		
		if(sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("KmerHitsCluster. Num unique: "+n+" median: "+median+" variance: "+variance+" stdev: "+stdev+" distance avg: "+dist.getAverage()+" stdev "+Math.sqrt(dist.getVariance()));
		int maxDistance = (int) Math.min(dist.getAverage(), stdev);
		if(maxDistance < 100) maxDistance=100;
		//if(maxDistance<0.01*query.length()) maxDistance*=2;
		else if (maxDistance>0.05*query.length()) maxDistance/=2;
		
		subjectPredictedStart = -1;
		for(List<UngappedSearchHit> hits:hitsMultiMap.values()) {
			UngappedSearchHit hit = selectHit(hits,median, Math.min(query.length()/20, maxDistance));
			if(hit!=null) {
				if(subjectPredictedStart==-1) {
					subjectPredictedStart = estimateSubjectStart(hit);
					subjectPredictedEnd = estimateSubjectEnd(hit);
					subjectEvidenceStart = hit.getStart();
					subjectEvidenceEnd = subjectEvidenceStart+hit.getQuery().length();
					queryPredictedStart = estimateQueryStart(hit);
					queryPredictedEnd = estimateQueryEnd(hit);
					queryEvidenceStart = hit.getQueryIdx();
					queryEvidenceEnd = hit.getQueryIdx() + hit.getQuery().length();
					
				}
				addHit(hit);
			}
		}
		if(sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("Final hits: "+hitsMap.size()+" start: "+subjectPredictedStart+" end: "+subjectPredictedEnd);
	}
	private UngappedSearchHit selectHit(List<UngappedSearchHit> hits, int median, int maxDistance) {
		UngappedSearchHit answer = null;
		int minDistance = 0;
		for(UngappedSearchHit hit:hits) {
			int estStart = estimateSubjectStart(hit);
			int distance = Math.abs(estStart-median);
			//if (sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" estq: "+estimateQueryStart(hit)+" - "+estimateQueryEnd(hit)+" estS: "+estimateSubjectStart(hit)+" - "+estimateSubjectEnd(hit)+" distance: "+distance+ " max: "+maxDistance);
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
		sequenceLength = kmerHit.getSequenceLength();
		int kmerQueryStart = kmerHit.getQueryIdx();
		subjectPredictedStart = estimateSubjectStart(kmerHit);
		subjectPredictedEnd = estimateSubjectEnd(kmerHit);
		subjectEvidenceStart = kmerHit.getStart();
		subjectEvidenceEnd = subjectEvidenceStart+kmerHit.getQuery().length();
		queryPredictedStart = estimateQueryStart(kmerHit);
		queryPredictedEnd = estimateQueryEnd(kmerHit);
		queryEvidenceStart = kmerHit.getQueryIdx();
		queryEvidenceEnd = kmerHit.getQueryIdx() + kmerHit.getQuery().length();
		hitsMap.put(kmerQueryStart, kmerHit);
		numDifferentKmers = 1;
		firstKmerPresent = queryEvidenceStart == 0;
		lastKmerPresent = queryEvidenceEnd==query.length();
	}
	
	private void addHit(UngappedSearchHit hit) {
		hitsMap.put(hit.getQueryIdx(), hit);
		int estStart = estimateSubjectStart(hit);
		int estEnd = estimateSubjectEnd(hit);
		if(estStart!=subjectPredictedStart || estEnd!=subjectPredictedEnd) allConsistent = false;
		subjectPredictedStart = Math.min(subjectPredictedStart, estStart);
		subjectPredictedEnd = Math.max(subjectPredictedEnd, estEnd);
		subjectEvidenceStart = Math.min(subjectEvidenceStart, hit.getStart());
		subjectEvidenceEnd = Math.max(subjectEvidenceEnd, hit.getStart()+hit.getQuery().length());
		queryPredictedStart = Math.min(queryPredictedStart, estimateQueryStart(hit));
		queryPredictedEnd = Math.max(queryPredictedEnd, estimateQueryEnd(hit));
		queryEvidenceStart = Math.min(queryEvidenceStart, hit.getQueryIdx());
		queryEvidenceEnd = Math.max(queryEvidenceEnd, hit.getQueryIdx() + hit.getQuery().length());
		if(queryEvidenceStart==0) firstKmerPresent = true;
		if (queryEvidenceEnd==query.length()) lastKmerPresent = true;
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
	
	private int estimateQueryStart(UngappedSearchHit hit) {
		return hit.getQueryIdx() - hit.getStart();
	}

	private int estimateQueryEnd(UngappedSearchHit hit) {
		return hit.getQueryIdx()+(sequenceLength-hit.getStart());
	}

	
	public void summarize(double averageHitsQuery) {
		this.averageHitsQuery = averageHitsQuery;
		numDifferentKmers = hitsMap.size();
		weightedCount = 0;
		
		List<UngappedSearchHit> hits = new ArrayList<UngappedSearchHit>();
		hits.addAll(hitsMap.values());
		for(UngappedSearchHit hit: hits) {
			double n = hit.getTotalHitsQuery();
			if(n<=averageHitsQuery) weightedCount++;
			else weightedCount += averageHitsQuery/n;
		}
		//Hits are already sorted by query id
		predictQueryStart (hits);
		predictQueryEnd (hits);
		Collections.sort(hits, (h1,h2)->h1.getStart()-h2.getStart());
		predictSubjectStart (hits);
		predictSubjectEnd (hits);
		predictOverlap (hits);
	}

	private void predictOverlap(List<UngappedSearchHit> hits) {
		List<Integer> estimatedOverlaps = new ArrayList<Integer>();
		for(UngappedSearchHit hit:hits) {
			int overlap = estimateOverlap (hit);
			estimatedOverlaps.add(overlap);
		}
		predictedOverlap = estimatedOverlaps.get(estimatedOverlaps.size()/2);
		if (subjectPredictedStart>0 && subjectPredictedEnd>sequenceLength && queryPredictedEnd<query.length()) {
			//Average with estimation from subject start
			predictedOverlap = (predictedOverlap+(sequenceLength-subjectPredictedStart)+queryPredictedEnd)/3;
		} else if (subjectPredictedStart<0 && subjectPredictedEnd<sequenceLength && queryPredictedStart>0) {
			//Average with estimation from subject end
			predictedOverlap = (predictedOverlap+subjectPredictedEnd+(query.length()-queryPredictedStart))/3;
		}
	}

	private void predictSubjectStart(List<UngappedSearchHit> hits) {
		double totalSumSubject = 0;
		double weightSum = 0;
		int n = hits.size();
		for(int i=0;i<n && i < 50;i++) {
			UngappedSearchHit hit = hits.get(i);
			double weight = ((double)(n-i))/n;
			weightSum += weight;
			totalSumSubject += weight*estimateSubjectStart(hit);
		}
		subjectPredictedStart = (int) Math.round(totalSumSubject / weightSum);
	}
	private void predictQueryStart(List<UngappedSearchHit> hits) {
		double totalSumQuery = 0;
		double weightSum = 0;
		int n = hits.size();
		for(int i=0;i<n && i < 50;i++) {
			UngappedSearchHit hit = hits.get(i);
			double weight = ((double)(n-i))/n;
			weightSum += weight;
			totalSumQuery += weight*estimateQueryStart(hit);
			//if (sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" estqStart: "+estimateQueryStart(hit)+" weight: "+weight+" sum: "+totalSumQuery);
		}
		queryPredictedStart = (int) Math.round(totalSumQuery / weightSum);
		if(sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("KmerHitsCluster. evidenceStart "+queryEvidenceStart+" predicted start: "+queryPredictedStart+" sum: "+totalSumQuery+" weight sum: "+weightSum);
		//int d = queryPredictedStart-queryEvidenceStart;
		//if (n>100 && d>=50) System.out.println("WARN: Predicted start "+queryPredictedStart+" for cluster to subject: "+sequenceIdx+" much larger than evidence: "+queryEvidenceStart+" query length: "+query.length());
	}

	private void predictSubjectEnd(List<UngappedSearchHit> hits) {
		double totalSumSubject = 0;
		double weightSum = 0;
		int n = hits.size();
		for(int i=Math.max(0, n-50);i<n;i++) {
			UngappedSearchHit hit = hits.get(i);
			double weight = ((double)i)/n;
			weightSum += weight;
			totalSumSubject += weight*estimateSubjectEnd(hit);
		}
		subjectPredictedEnd = (int) Math.round(totalSumSubject / weightSum);
	}
	private void predictQueryEnd(List<UngappedSearchHit> hits) {
		double totalSumQuery = 0;
		double weightSum = 0;
		int n = hits.size();
		for(int i=Math.max(0, n-50);i<n;i++) {
			UngappedSearchHit hit = hits.get(i);
			double weight = ((double)i)/n;
			weightSum += weight;
			totalSumQuery += weight*estimateQueryEnd(hit);
			//if (sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" estqEnd: "+estimateQueryEnd(hit)+" weight: "+weight+" sum: "+totalSumQuery);
		}
		queryPredictedEnd = (int) Math.round(totalSumQuery / weightSum);
		if(sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("KmerHitsCluster. evidenceEnd "+queryEvidenceEnd+" predicted end: "+queryPredictedEnd+" sum: "+totalSumQuery+" weight sum: "+weightSum);
		//int d = queryEvidenceEnd-queryPredictedEnd;
		//if (n>100 && d>=50) System.out.println("WARN: Predicted end "+queryPredictedEnd+" for cluster to subject: "+sequenceIdx+" much smaller than evidence: "+queryEvidenceEnd+" query length: "+query.length());
	}
	
	private int estimateOverlap(UngappedSearchHit hit) {
		int start = estimateSubjectStart(hit);
		int end = estimateSubjectEnd(hit);
		int qS = estimateQueryStart(hit);
		int qE = estimateQueryEnd(hit);
		int overlap = Math.min(query.length(), sequenceLength);
		overlap = Math.min(overlap, sequenceLength-start);
		overlap = Math.min(overlap, end);
		overlap = Math.min(overlap, query.length()-qS);
		overlap = Math.min(overlap, qE);
		
		return overlap;
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
	public int getSubjectEvidenceStart() {
		return subjectEvidenceStart;
	}

	public int getSubjectEvidenceEnd() {
		return subjectEvidenceEnd;
	}
	public void setSubjectPredictedLimits (int predictedStart, int predictedEnd) {
		if(predictedEnd<=predictedStart) throw new IllegalArgumentException("Predicted end "+predictedEnd+" must be larger than predicted start: "+predictedStart);
		subjectPredictedStart = predictedStart;
		subjectPredictedEnd = predictedEnd;
	}
	
	public int getQueryPredictedStart() {
		return queryPredictedStart;
	}
	public int getQueryPredictedEnd() {
		return queryPredictedEnd;
	}
	public int getQueryEvidenceStart() {
		return queryEvidenceStart;
	}
	public int getQueryEvidenceEnd() {
		return queryEvidenceEnd;
	}
	
	public int getPredictedOverlap() {
		return predictedOverlap;
	}

	/**
	 * @return the query
	 */
	public CharSequence getQuery() {
		return query;
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
	
	public int getSelfHitsCountQuery() {
		return selfHitsCountQuery;
	}
	public void setSelfHitsCountQuery(int selfHitsCountQuery) {
		this.selfHitsCountQuery = selfHitsCountQuery;
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
	
	public void disposeHits () {
		hitsMap.clear();
	}
	
	public static List<KmerHitsCluster> clusterRegionKmerAlns(CharSequence query, List<UngappedSearchHit> sequenceHits, double minQueryCoverage) {
		List<KmerHitsCluster> answer = new ArrayList<>();
		Collections.sort(sequenceHits,(h1,h2) -> h1.getStart()-h2.getStart());
		KmerHitsCluster uniqueCluster = new KmerHitsCluster(query, sequenceHits);
		//if(querySequenceId==idxDebug) System.out.println("Hits to cluster: "+sequenceKmerHits.size()+" target: "+uniqueCluster.getSequenceIdx()+" first: "+uniqueCluster.getFirst()+" last: "+uniqueCluster.getLast()+" kmers: "+uniqueCluster.getNumDifferentKmers());
		if (uniqueCluster.getQueryEvidenceEnd()-uniqueCluster.getQueryEvidenceStart()<minQueryCoverage*query.length()) return answer;
		answer.add(uniqueCluster);
		if(uniqueCluster.getNumDifferentKmers()>0.8*sequenceHits.size()) return answer;
		//Cluster remaining hits
		List<UngappedSearchHit> remainingHits = new ArrayList<UngappedSearchHit>();
		for(UngappedSearchHit hit:sequenceHits) {
			if (hit!=uniqueCluster.getKmerHit(hit.getQueryIdx())) remainingHits.add(hit);
		}
		KmerHitsCluster cluster2 = new KmerHitsCluster(query, remainingHits);
		if(cluster2.getQueryEvidenceEnd()-cluster2.getQueryEvidenceStart()>=minQueryCoverage*query.length()) answer.add(cluster2);
		return answer;
	}
	
}