package ngsep.sequences;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.math.Distribution;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class KmerHitsCluster {

	private int queryLength;
	private int kmerLength;
	//Map indexed by query start
	private Map<Integer,UngappedSearchHit> hitsMap=new TreeMap<Integer, UngappedSearchHit>();
	private int subjectIdx;
	private String subjectName;
	private int subjectLength;
	private int subjectPredictedStart;
	private int subjectPredictedEnd;
	private int queryPredictedStart;
	private int queryPredictedEnd;
	private int queryEvidenceStart;
	private int queryEvidenceEnd;
	private int subjectEvidenceStart;
	private int subjectEvidenceEnd;
	private int predictedOverlap;
	private int averagePredictedOverlap;
	private int medianPredictedOverlap;
	private int fromLimitsPredictedOverlap;
	private double subjectStartSD=0;
	private double predictedOverlapSD=0;
	private int numDifferentKmers = 0;
	private double weightedCount=0;
	private boolean allConsistent = true;
	private boolean firstKmerPresent = false;
	private boolean lastKmerPresent = false;
	private static int idxSubjectDebug = 371;
	private static int queryLengthDebug = 25300;
	
	public KmerHitsCluster(int queryLength, int subjectLength, List<UngappedSearchHit> inputHits) {
		this.queryLength = queryLength;
		if(inputHits.size()==0) throw new RuntimeException("Invalid empty input hits for cluster"+inputHits.size());
		
		UngappedSearchHit firstHit = inputHits.get(0);
		subjectIdx = firstHit.getSequenceIdx();
		subjectName = firstHit.getSequenceName();
		kmerLength = firstHit.getQuery().length();
		this.subjectLength = subjectLength;
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("KmerHitsCluster. Clustering "+inputHits.size()+" hits. Subject idx: "+subjectIdx);
		
		//Index hits by query kmer start
		Map<Integer,List<UngappedSearchHit>> hitsMultiMap = new TreeMap<Integer, List<UngappedSearchHit>>();
		
		for(UngappedSearchHit hit:inputHits) {
			//if (sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" estq: "+estimateQueryStart(hit)+" - "+estimateQueryEnd(hit)+" estS: "+estimateSubjectStart(hit)+" - "+estimateSubjectEnd(hit));
			List<UngappedSearchHit> list = hitsMultiMap.computeIfAbsent(hit.getQueryIdx(), l -> new ArrayList<UngappedSearchHit>());
			list.add(hit);
		}
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("KmerHitsCluster. Num different kmers: "+hitsMultiMap.size());
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
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("KmerHitsCluster. Num unique: "+n);
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
		
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("KmerHitsCluster. Num unique: "+n+" median: "+median+" variance: "+variance+" stdev: "+stdev+" distance avg: "+dist.getAverage()+" stdev "+Math.sqrt(dist.getVariance()));
		int maxDistance = 5*queryLength; 
		if(subjectLength>maxDistance) {
			// This is only useful for mapping to a long reference subject
			maxDistance = (int) Math.min(dist.getAverage(), stdev);
			if(maxDistance < 100) maxDistance=100;
			//if(maxDistance<0.01*query.length()) maxDistance*=2;
			else if (maxDistance>0.05*queryLength) maxDistance/=2;
		}
		
		
		subjectPredictedStart = -1;
		for(List<UngappedSearchHit> hits:hitsMultiMap.values()) {
			UngappedSearchHit hit = selectHit(hits,median, Math.min(queryLength/20, maxDistance));
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
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Final hits: "+hitsMap.size()+" start: "+subjectPredictedStart+" end: "+subjectPredictedEnd);
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



	public KmerHitsCluster(int queryLength, int subjectLength, UngappedSearchHit kmerHit) {
		this.queryLength = queryLength;
		subjectIdx = kmerHit.getSequenceIdx();
		subjectName = kmerHit.getSequenceName();
		this.subjectLength = subjectLength;
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
		lastKmerPresent = queryEvidenceEnd==queryLength;
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
		if (queryEvidenceEnd==queryLength) lastKmerPresent = true;
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
		return hit.getStart()+(queryLength-hit.getQueryIdx());
	}
	
	private int estimateQueryStart(UngappedSearchHit hit) {
		return hit.getQueryIdx() - hit.getStart();
	}

	private int estimateQueryEnd(UngappedSearchHit hit) {
		return hit.getQueryIdx()+(subjectLength-hit.getStart());
	}

	
	public void summarize() {
		numDifferentKmers = hitsMap.size();
		weightedCount = 0;
		List<UngappedSearchHit> hits = new ArrayList<UngappedSearchHit>();
		hits.addAll(hitsMap.values());
		for(UngappedSearchHit hit: hits) {
			weightedCount += hit.getWeight();
		}
		//Hits are already sorted by query id
		predictQueryStart (hits);
		predictQueryEnd (hits);
		Collections.sort(hits, (h1,h2)->h1.getStart()-h2.getStart());
		predictSubjectStart (hits);
		predictSubjectEnd (hits);
		predictOverlap (hits);
		calculateSubjectStartSD(hits);
	}

	private void calculateSubjectStartSD(List<UngappedSearchHit> hits) {
		double sum=0;
		double sum2=0;
		int n = hits.size();
		for(UngappedSearchHit hit:hits) {
			int start = estimateSubjectStart(hit);
			sum+=start;
			sum2+=(start*start);
		}
		double startVariance = (sum2-sum*sum/n)/n-1;
		subjectStartSD = Math.sqrt(startVariance);
		
	}
	private void predictOverlap(List<UngappedSearchHit> hits) {
		List<Integer> estimatedOverlaps = new ArrayList<Integer>();
		double sum=0;
		double sum2=0;
		int n = hits.size();
		for(UngappedSearchHit hit:hits) {
			int overlap = estimateOverlap (hit);
			estimatedOverlaps.add(overlap);
			sum+=overlap;
			sum2+=overlap*overlap;
		}
		averagePredictedOverlap = (int)Math.round(sum/n);
		medianPredictedOverlap = estimatedOverlaps.get(estimatedOverlaps.size()/2);
		double predictedOverlapVariance = (sum2-sum*sum/n)/n-1;
		predictedOverlapSD = Math.sqrt(predictedOverlapVariance);
		predictedOverlap = medianPredictedOverlap;
		fromLimitsPredictedOverlap = 0;
		if (subjectPredictedStart>0 && subjectPredictedEnd>subjectLength && queryPredictedEnd<queryLength) {
			//Average with estimation from subject start
			fromLimitsPredictedOverlap = ((subjectLength-subjectPredictedStart)+queryPredictedEnd)/2;
			//predictedOverlap = (predictedOverlap+(sequenceLength-subjectPredictedStart)+queryPredictedEnd)/3;
		} else if (subjectPredictedStart<0 && subjectPredictedEnd<subjectLength && queryPredictedStart>0) {
			//Average with estimation from subject end
			fromLimitsPredictedOverlap = (subjectPredictedEnd+(queryLength-queryPredictedStart))/2;
			//predictedOverlap = (predictedOverlap+subjectPredictedEnd+(query.length()-queryPredictedStart))/3;
		}
		if(fromLimitsPredictedOverlap>0) predictedOverlap=fromLimitsPredictedOverlap;
		else predictedOverlap = averagePredictedOverlap;
		
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
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("KmerHitsCluster. evidenceStart "+queryEvidenceStart+" predicted start: "+queryPredictedStart+" sum: "+totalSumQuery+" weight sum: "+weightSum);
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
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("KmerHitsCluster. evidenceEnd "+queryEvidenceEnd+" predicted end: "+queryPredictedEnd+" sum: "+totalSumQuery+" weight sum: "+weightSum);
		//int d = queryEvidenceEnd-queryPredictedEnd;
		//if (n>100 && d>=50) System.out.println("WARN: Predicted end "+queryPredictedEnd+" for cluster to subject: "+sequenceIdx+" much smaller than evidence: "+queryEvidenceEnd+" query length: "+query.length());
	}
	
	private int estimateOverlap(UngappedSearchHit hit) {
		int start = estimateSubjectStart(hit);
		int end = estimateSubjectEnd(hit);
		int qS = estimateQueryStart(hit);
		int qE = estimateQueryEnd(hit);
		int overlap = Math.min(queryLength, subjectLength);
		overlap = Math.min(overlap, subjectLength-start);
		overlap = Math.min(overlap, end);
		overlap = Math.min(overlap, queryLength-qS);
		overlap = Math.min(overlap, qE);
		
		return overlap;
	}

	public String getSubjectName() {
		return subjectName;
	}
	
	/**
	 * @return the sequenceIdx
	 */
	public int getSubjectIdx() {
		return subjectIdx;
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
	
	public int getAveragePredictedOverlap() {
		return averagePredictedOverlap;
	}
	public int getMedianPredictedOverlap() {
		return medianPredictedOverlap;
	}
	public int getFromLimitsPredictedOverlap() {
		return fromLimitsPredictedOverlap;
	}
	public double getPredictedOverlapSD() {
		return predictedOverlapSD;
	}
	
	public double getSubjectStartSD() {
		return subjectStartSD;
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
	
	public void disposeHits () {
		hitsMap.clear();
	}
	
	public static List<KmerHitsCluster> clusterRegionKmerAlns(int queryLength, int subjectLength, List<UngappedSearchHit> sequenceHits, double minQueryCoverage) {
		List<KmerHitsCluster> answer = new ArrayList<>();
		KmerHitsCluster uniqueCluster = new KmerHitsCluster(queryLength, subjectLength, sequenceHits);
		//if(querySequenceId==idxDebug) System.out.println("Hits to cluster: "+sequenceKmerHits.size()+" target: "+uniqueCluster.getSequenceIdx()+" first: "+uniqueCluster.getFirst()+" last: "+uniqueCluster.getLast()+" kmers: "+uniqueCluster.getNumDifferentKmers());
		if (uniqueCluster.getQueryEvidenceEnd()-uniqueCluster.getQueryEvidenceStart()<minQueryCoverage*queryLength) return answer;
		answer.add(uniqueCluster);
		if(uniqueCluster.getNumDifferentKmers()>0.8*sequenceHits.size()) return answer;
		//Cluster remaining hits
		List<UngappedSearchHit> remainingHits = new ArrayList<UngappedSearchHit>();
		for(UngappedSearchHit hit:sequenceHits) {
			if (hit!=uniqueCluster.getKmerHit(hit.getQueryIdx())) remainingHits.add(hit);
		}
		KmerHitsCluster cluster2 = new KmerHitsCluster(queryLength, subjectLength, remainingHits);
		if(cluster2.getQueryEvidenceEnd()-cluster2.getQueryEvidenceStart()>=minQueryCoverage*queryLength) answer.add(cluster2);
		return answer;
	}
	public void completeMissingHits(String subjectSequence, Map<Integer, Long> queryCodes) {
		//Map<Long,Integer> subjectUniqueCodes = KmersExtractor.extractLocallyUniqueKmerCodes(subjectSequence, kmerLength, subjectEvidenceStart, subjectEvidenceEnd);
		Map<Integer,Long> subjectCodes = KmersExtractor.extractDNAKmerCodes(subjectSequence, kmerLength, subjectEvidenceStart, subjectEvidenceEnd);
		
		Map<Long,List<Integer>> subjectCodesPos = new HashMap<Long, List<Integer>>();
		for(int i:subjectCodes.keySet()) {
			long code = subjectCodes.get(i);
			List<Integer> posList = subjectCodesPos.computeIfAbsent(code, v->new ArrayList<Integer>());
			posList.add(i);
		}
		int lastQueryStart = -1;
		int lastSubjectStart = -1;
		UngappedSearchHit lastHit = null;
		ArrayList<Integer> queryStarts = new ArrayList<Integer>(hitsMap.size());
		queryStarts.addAll(hitsMap.keySet());
		for(int queryStart:queryStarts) {
			UngappedSearchHit hit = hitsMap.get(queryStart);
			int subjectStart = hit.getStart();
			int diffQ = queryStart-lastQueryStart;
			int diffS = subjectStart - lastSubjectStart;
			int diffC = Math.abs(diffQ-diffS);
			if(lastHit!=null && diffC<10 && diffQ<100 && diffS>0) {
				for(int i=lastQueryStart+1;i<queryStart;i++) {
					Long code = queryCodes.get(i);
					if(code == null) continue;
					List<Integer> posList = subjectCodesPos.get(code);
					if(posList==null) continue;
					Integer selectedPos = null;
					for(int subjectPos:posList) {
						if(lastSubjectStart<subjectPos && subjectPos<subjectStart) {
							selectedPos = subjectPos;
						}
					}
					if(selectedPos==null) continue;
					CharSequence kmer = new String(AbstractLimitedSequence.getSequence(code, kmerLength, DNASequence.EMPTY_DNA_SEQUENCE));
					UngappedSearchHit rescuedHit = new UngappedSearchHit(kmer, subjectIdx, selectedPos);
					rescuedHit.setQueryIdx(i);
					double weight = (lastHit.getWeight()+hit.getWeight())/2;
					rescuedHit.setWeight(weight);
					addHit(rescuedHit);
				}
			}
			lastQueryStart = queryStart;
			lastSubjectStart = subjectStart;
			lastHit = hit;
		}
	}
	
}