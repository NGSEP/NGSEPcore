package ngsep.sequences;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
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
	private double weightedCount=0;
	private int rawKmerHits=0;
	private double rawKmerHitsSubjectStartSD = 0;
	private boolean allConsistent = true;
	private boolean firstKmerPresent = false;
	private boolean lastKmerPresent = false;
	private static int idxSubjectDebug = -1;
	private static int queryLengthDebug = -1;
	
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
		rawKmerHits = inputHits.size();
		for(UngappedSearchHit hit:inputHits) {
			//if (sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" estq: "+estimateQueryStart(hit)+" - "+estimateQueryEnd(hit)+" estS: "+estimateSubjectStart(hit)+" - "+estimateSubjectEnd(hit));
			List<UngappedSearchHit> list = hitsMultiMap.computeIfAbsent(hit.getQueryIdx(), l -> new ArrayList<UngappedSearchHit>());
			list.add(hit);
		}
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("KmerHitsCluster. Num different kmers: "+hitsMultiMap.size());
		List<Integer> subjectStarts = new ArrayList<Integer>();
		double sum = 0;
		double sum2 = 0;
		double n = 0;
		for(List<UngappedSearchHit> hits:hitsMultiMap.values()) {
			for(UngappedSearchHit hit:hits) {
				int estStart = estimateSubjectStart(hit);
				subjectStarts.add(estStart);
				sum+=1.0*estStart;
				sum2+=(1.0*estStart*estStart);
				n++;
			}
		}
		
		Collections.sort(subjectStarts);
		int median = subjectStarts.get(subjectStarts.size()/2);
		//System.out.println("Sum: "+sum+" sum2: "+sum2);
		double variance = (sum2-sum*sum/n)/(n-1);
		rawKmerHitsSubjectStartSD = (variance>0)?Math.sqrt(variance):0;
		Distribution dist = new Distribution(0, 100, 1);
		for(int start:subjectStarts) {
			int distance = Math.abs(start-median);
			if (distance < rawKmerHitsSubjectStartSD) dist.processDatapoint(distance);
		}
		//System.out.println(subjectStarts);
		
		int maxDistance = 5*queryLength; 
		if(subjectLength>maxDistance) {
			// This is only useful for mapping to a long reference subject
			maxDistance = (int) Math.max(dist.getAverage(), rawKmerHitsSubjectStartSD);
			maxDistance *=3;
			if(maxDistance < 100) maxDistance=100;
			//if(maxDistance<0.01*query.length()) maxDistance*=2;
			else if (maxDistance>0.05*queryLength) maxDistance/=2;
		} else if (queryLength>2000) {
			maxDistance = Math.min(queryLength/20, maxDistance);
		}
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("KmerHitsCluster. Num hits: "+n+" median: "+median+" variance: "+variance+" stdev: "+rawKmerHitsSubjectStartSD+" distance avg: "+dist.getAverage()+" stdev "+Math.sqrt(dist.getVariance())+" max distance: "+maxDistance);
		
		List<UngappedSearchHit> selectedHits = new ArrayList<UngappedSearchHit>(hitsMultiMap.size());
		int minSubjectStart = firstHit.getStart();
		int maxSubjectStart = firstHit.getStart();
		for(List<UngappedSearchHit> hits:hitsMultiMap.values()) {
			UngappedSearchHit hit = selectHit(hits, queryLength, median, maxDistance);
			if(hit!=null) {
				if (hit.getSequenceIdx()==idxSubjectDebug && queryLength == queryLengthDebug /*&& hit.getQueryIdx()<1000*/) System.out.println("Selected hits. Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" estq: "+estimateQueryStart(hit)+" - "+estimateQueryEnd(hit)+" estS: "+estimateSubjectStart(hit)+" - "+estimateSubjectEnd(hit));
				selectedHits.add(hit);
				minSubjectStart = Math.min(minSubjectStart, hit.getStart());
				maxSubjectStart = Math.max(maxSubjectStart, hit.getStart());
			}
		}
		if(selectedHits.size()<1) System.err.println("WARN. Empty list of selected hits for subject: "+subjectIdx+" "+subjectName);
		List<UngappedSearchHit> filteredHits = removeDisorganized (selectedHits);
		if(filteredHits.size()<1) {
			System.err.println("WARN. Empty list of sorted hits for subject: "+subjectIdx+" "+subjectName+" selected hits: "+selectedHits.size()+" query length: "+queryLength);
			return;
		}
		//Create cluster with selected hits
		boolean initialized = false;
		for(UngappedSearchHit hit:filteredHits) {
			int distance = Math.abs(estimateSubjectStart(hit)-median);
			if (hit.getSequenceIdx()==idxSubjectDebug && queryLength == queryLengthDebug /*&& hit.getQueryIdx()<1000*/) System.out.println("Filtered hits. Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" estq: "+estimateQueryStart(hit)+" - "+estimateQueryEnd(hit)+" estS: "+estimateSubjectStart(hit)+" - "+estimateSubjectEnd(hit)+" distance: "+distance+ " max: "+maxDistance);
			if(!initialized) {
				subjectPredictedStart = estimateSubjectStart(hit);
				subjectPredictedEnd = estimateSubjectEnd(hit);
				subjectEvidenceStart = hit.getStart();
				subjectEvidenceEnd = subjectEvidenceStart+hit.getQuery().length();
				queryPredictedStart = estimateQueryStart(hit);
				queryPredictedEnd = estimateQueryEnd(hit);
				queryEvidenceStart = hit.getQueryIdx();
				queryEvidenceEnd = hit.getQueryIdx() + hit.getQuery().length();
				initialized = true;
			}
			addHit(hit);
		}
		if(filteredHits.size()>0) summarize();
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Final hits: "+hitsMap.size()+" start: "+subjectPredictedStart+" end: "+subjectPredictedEnd);
	}
	
	private UngappedSearchHit selectHit(List<UngappedSearchHit> hits, int queryLength, int median, int maxDistance) {
		UngappedSearchHit answer = null;
		int minDistance = 0;
		for(UngappedSearchHit hit:hits) {
			int estStart = estimateSubjectStart(hit);
			int distance = Math.abs(estStart-median);
			if(distance <= maxDistance) {
				if(answer == null || minDistance>distance) {
					answer = hit;
					minDistance = distance;
				}
			}
		}
		return answer;
	}

	private List<UngappedSearchHit> removeDisorganized(List<UngappedSearchHit> selectedHits) {
		int n = selectedHits.size();
		if(n<=30) return selectedHits;
		boolean [] vicinityConsistent = new boolean [n];
		Arrays.fill(vicinityConsistent, false);
		for(int i=15;i<n-15;i++) {
			UngappedSearchHit hit = selectedHits.get(i);
			int start = hit.getStart();
			boolean consistent = true;
			for(int j=i-1;j>=i-15 && consistent;j--) {
				UngappedSearchHit hit2 = selectedHits.get(j);
				if(hit2.getStart()>=start) consistent = false; 
			}
			for(int j=i+1;j<=i+15 && consistent;j++) {
				UngappedSearchHit hit2 = selectedHits.get(j);
				if(hit2.getStart()<=start) consistent = false; 
			}
			vicinityConsistent[i] = consistent;
		}
		List<UngappedSearchHit> dpSortedHits = new ArrayList<UngappedSearchHit>();
		double sum = 0;
		double sum2 = 0;
		double nS = 0;
		int i=0;
		while(i<n) {
			UngappedSearchHit hit = selectedHits.get(i);
			if(vicinityConsistent[i]) {
				dpSortedHits.add(hit);
				double s = estimateSubjectStart(hit);
				sum+=s;
				sum2+=s*s;
				nS++;
				i++;
			} else {
				int j;
				for(j=i+1;j<n && !vicinityConsistent[j];j++);
				dpSortedHits.addAll(selectSorted(selectedHits,i,j-1));
				i=j;
			}
		}
		//Remove final outliers
		int n2= dpSortedHits.size();
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Initial hits: "+n+" consistent: "+nS+" dpSorted: "+n2);
		if(n2<30 || nS < 20) return dpSortedHits;
		int averageStart = (int) (sum/nS);
		double variance = Math.max(1, (sum2-sum*sum/nS)/(nS-1));
		double stdev = Math.sqrt(variance);
		double maxDistance = Math.max(30, 3*stdev);
		int countOutliers =0;
		for(UngappedSearchHit hit:dpSortedHits) {
			double distance = Math.abs(estimateSubjectStart(hit)-averageStart);
			if(distance > maxDistance) countOutliers++;
		}
		
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Average start: "+averageStart+" stdev: "+stdev+" max distance: "+maxDistance+" outliers: "+countOutliers+" n "+n2+" nconsistent: "+nS);
		List<UngappedSearchHit> answer = new ArrayList<UngappedSearchHit>(n2);
		if(countOutliers<0.02*n2) {
			for(UngappedSearchHit hit:dpSortedHits) {
				double distance = Math.abs(estimateSubjectStart(hit)-averageStart);
				if(distance <= maxDistance) answer.add(hit);
				else if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Removing outlier: "+hit.getQueryIdx());
			}
			return answer;
		}
		Collections.sort(dpSortedHits,(h1,h2)->Math.abs(estimateSubjectStart(h2)-averageStart)-Math.abs(estimateSubjectStart(h1)-averageStart));
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Sorted by diestance");
		for (int k=0;k<dpSortedHits.size();k++) {
			UngappedSearchHit hit = dpSortedHits.get(k);
			double distance = Math.abs(estimateSubjectStart(hit)-averageStart);
			if(distance <= maxDistance || k>=0.02*n2) answer.add(hit);
			else if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Removing outlier after sorting: "+hit.getQueryIdx());
		}
		Collections.sort(answer,(h1,h2)->h1.getQueryIdx()-h2.getQueryIdx());
		return answer;
	}


	private List<UngappedSearchHit> selectSorted(List<UngappedSearchHit> selectedHits, int first, int last) {
		//Check sorted first
		List<UngappedSearchHit> regionHits = new ArrayList<UngappedSearchHit>(last-first+1);
		int lastStartSubject = -1;
		int lastStartQuery = -1;
		boolean sortedSubject = true;
		boolean sortedQuery = true;
		for(int i=first;i<=last;i++) {
			UngappedSearchHit hit = selectedHits.get(i);
			regionHits.add(hit);
			int startS = hit.getStart();
			if(lastStartSubject>=startS) sortedSubject = false;
			lastStartSubject = startS;
			int startQ = hit.getQueryIdx();
			if(lastStartQuery>=startQ) sortedQuery = false;
			lastStartQuery = startQ;
		}
		if(!sortedQuery) {
			System.err.println("WARN. Disorganized query hits list");
			Collections.sort(regionHits,(h1,h2)->h1.getQueryIdx()-h2.getQueryIdx());
			
		} else if (sortedSubject) return regionHits;
		List<UngappedSearchHit> regionHitsBySubject = new ArrayList<UngappedSearchHit>();
		regionHitsBySubject.addAll(regionHits);
		Collections.sort(regionHitsBySubject,(h1,h2)->h1.getStart()-h2.getStart());
		int n = regionHits.size();
		//Dynamic programming to select the best organized subset
		int [][] scores = new int [n+1][n+1];
		for(int i=0;i<=n;i++) {
			for(int j=0;j<=n;j++) {
				if(i==0 || j==0) scores[i][j] = 0;
				else {
					UngappedSearchHit hit1 = regionHits.get(i-1);
					UngappedSearchHit hit2 = regionHitsBySubject.get(j-1);
					if(hit1==hit2) scores[i][j] = 1 + scores[i-1][j-1];
					else scores[i][j] = Math.max(scores[i-1][j], scores[i][j-1]);
				}
			}
			 
		}
		List<UngappedSearchHit> answer = new LinkedList<UngappedSearchHit>();
		int i= n;
		int j = n;
		while(i>0 && j>0) {
			UngappedSearchHit hit1 = regionHits.get(i-1);
			UngappedSearchHit hit2 = regionHitsBySubject.get(j-1);
			if(hit1==hit2) {
				answer.add(0,hit1);
				i--;
				j--;
			} else if (scores[i][j] == scores[i-1][j]) i--;
			else j--;
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
	
	private static int estimateSubjectStart(UngappedSearchHit hit) {
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
		predictedOverlap = averagePredictedOverlap;
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
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Predicting overlap. Avg: "+averagePredictedOverlap+" median: "+medianPredictedOverlap+" from limits: "+fromLimitsPredictedOverlap);
		//if(fromLimitsPredictedOverlap>0) predictedOverlap=fromLimitsPredictedOverlap;
		//else predictedOverlap = averagePredictedOverlap;
		//predictedOverlap = (averagePredictedOverlap+medianPredictedOverlap)/2;
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
	
	public int getRawKmerHits() {
		return rawKmerHits;
	}
	public double getRawKmerHitsSubjectStartSD() {
		return rawKmerHitsSubjectStartSD;
	}
	/**
	 * 
	 * @return the kmerNumbers
	 */
	public int getNumDifferentKmers() {
		return hitsMap.size();
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
		Map<Integer,Integer> countsByQueryIdx = new HashMap<Integer, Integer>();
		double minHits = Math.min(20,0.01*queryLength);
		if(sequenceHits.size()<minHits) return new ArrayList<KmerHitsCluster>();
		for(UngappedSearchHit hit:sequenceHits) {
			//if (sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" estq: "+estimateQueryStart(hit)+" - "+estimateQueryEnd(hit)+" estS: "+estimateSubjectStart(hit)+" - "+estimateSubjectEnd(hit));
			countsByQueryIdx.compute(hit.getQueryIdx(), (k,v)->v==null?1:v+1);
		}
		//double estimatedClusters = 0.5*sequenceHits.size()/queryLength;
		double avg = 0;
		for(int count:countsByQueryIdx.values()) avg+=count;
		if(avg>0) avg/=countsByQueryIdx.size();
		double estimatedClusters = avg;
		UngappedSearchHit firstHit = sequenceHits.get(0);
		int subjectIdx = firstHit.getSequenceIdx();
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Clustering hits: "+sequenceHits.size()+" estimatedCLusters: "+estimatedClusters);
		//TODO: Use good clustering also for mapping
		if(minQueryCoverage==0) return clusterRegionKmerAlnsMultiple(queryLength, subjectLength, sequenceHits, estimatedClusters);
		//if(estimatedClusters>1.5) return clusterRegionKmerAlnsMultiple(queryLength, subjectLength, sequenceHits, estimatedClusters);
		List<KmerHitsCluster> answer = new ArrayList<>();
		KmerHitsCluster uniqueCluster = new KmerHitsCluster(queryLength, subjectLength, sequenceHits);
		//if(querySequenceId==idxDebug) System.out.println("Hits to cluster: "+sequenceKmerHits.size()+" target: "+uniqueCluster.getSequenceIdx()+" first: "+uniqueCluster.getFirst()+" last: "+uniqueCluster.getLast()+" kmers: "+uniqueCluster.getNumDifferentKmers());
		if (uniqueCluster.getNumDifferentKmers()<minHits || uniqueCluster.getQueryEvidenceEnd()-uniqueCluster.getQueryEvidenceStart()<minQueryCoverage*queryLength) return answer;
		answer.add(uniqueCluster);
		if(uniqueCluster.getNumDifferentKmers()>0.8*sequenceHits.size()) return answer;
		//Cluster remaining hits
		List<UngappedSearchHit> remainingHits = new ArrayList<UngappedSearchHit>();
		for(UngappedSearchHit hit:sequenceHits) {
			if (hit!=uniqueCluster.getKmerHit(hit.getQueryIdx())) remainingHits.add(hit);
		}
		KmerHitsCluster cluster2 = new KmerHitsCluster(queryLength, subjectLength, remainingHits);
		if(cluster2.getNumDifferentKmers()>=minHits && cluster2.getQueryEvidenceEnd()-cluster2.getQueryEvidenceStart()>=minQueryCoverage*queryLength) answer.add(cluster2);
		return answer;
	}
	private static List<KmerHitsCluster> clusterRegionKmerAlnsMultiple(int queryLength, int subjectLength, List<UngappedSearchHit> sequenceHits, double estimatedClusters) {
		List<KmerHitsCluster> answer = new ArrayList<>();
		UngappedSearchHit firstHit = sequenceHits.get(0);
		int subjectIdx = firstHit.getSequenceIdx();
		//Initial clustering
		Map<Integer,List<UngappedSearchHit>> hitsByBin = new HashMap<Integer, List<UngappedSearchHit>>();
		for(UngappedSearchHit hit:sequenceHits) {
			int estStart = estimateSubjectStart(hit);
			int bin = estStart/1000;
			if(estStart<0) bin--;
			List<UngappedSearchHit> hitsBin = hitsByBin.computeIfAbsent(bin, (v)-> new ArrayList<UngappedSearchHit>());
			hitsBin.add(hit);
		}
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Clustering kmer hits. Initial bin starts: "+hitsByBin.keySet()+" estimated number of clusters: "+estimatedClusters);
		//Second cluster centered in best averages
		List<List<UngappedSearchHit>> sortedClusters = new ArrayList<List<UngappedSearchHit>>(hitsByBin.size());
		sortedClusters.addAll(hitsByBin.values());
		Collections.sort(sortedClusters,(c1,c2)->c2.size()-c1.size());
		hitsByBin.clear();
		List<Integer> clusterAverages = new ArrayList<Integer>((int) (2*estimatedClusters+1));
		for(int i=0;i<sortedClusters.size() && i<=2*(estimatedClusters);i++) {
			List<UngappedSearchHit> cluster = sortedClusters.get(i);
			if(cluster.size()<20) break;
			int average = getAverageEstimatedSubjectStart(cluster);
			if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Clustering kmer hits. Average predicted start next cluster: "+average+" size: "+cluster.size());
			clusterAverages.add(average);
		}
		if (clusterAverages.size()==0) return answer;
		Collections.sort(clusterAverages, (c1,c2)->c1-c2);
		int next = clusterAverages.get(0);
		for(int average:clusterAverages) {
			if(average-next<500) next = (next+average)/2;
			else {
				hitsByBin.put(next, new ArrayList<UngappedSearchHit>());
				next = average;
			}
		}
		hitsByBin.put(next, new ArrayList<UngappedSearchHit>());
		for(UngappedSearchHit hit:sequenceHits) {
			int estStart = estimateSubjectStart(hit);
			int minS = 0;
			int minD = 500;
			for(int start:hitsByBin.keySet()) {
				int d = Math.abs(estStart-start);
				if(d<minD) {
					minS = start;
					minD = d;
				}
			}
			if(minD<500) hitsByBin.get(minS).add(hit);
		}
		if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Clustering kmer hits. Final bin starts: "+hitsByBin.keySet());
		for(List<UngappedSearchHit> hits:hitsByBin.values()) {
			List<List<UngappedSearchHit>> subclusters = breakByQueryStart(hits);
			for(List<UngappedSearchHit> subcluster:subclusters) {
				if(subcluster.size()<20) continue;
				KmerHitsCluster cluster = new KmerHitsCluster(queryLength, subjectLength, subcluster);
				if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Next cluster subject predicted coords: "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" subject evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" query evidence: "+cluster.getQueryEvidenceStart()+" "+cluster.getQueryEvidenceEnd()+" unique kmers: "+cluster.getNumDifferentKmers());
				if(cluster.getNumDifferentKmers()>=20) answer.add(cluster);
			}
		}
		return answer;
	}
	private static List<List<UngappedSearchHit>> breakByQueryStart(List<UngappedSearchHit> hits) {
		List<List<UngappedSearchHit>> answer = new ArrayList<List<UngappedSearchHit>>();
		Collections.sort(hits, (h1,h2)->h1.getQueryIdx()-h2.getQueryIdx());
		List<UngappedSearchHit> nextSubcluster = new ArrayList<UngappedSearchHit>();
		int nextQueryStart = 0;
		for(UngappedSearchHit hit:hits) {
			//TODO: define better this parameter
			if(hit.getQueryIdx()>nextQueryStart+30000) {
				if(nextSubcluster.size()>0) answer.add(nextSubcluster); 
				nextSubcluster = new ArrayList<UngappedSearchHit>();
			}
			nextSubcluster.add(hit);
			nextQueryStart = hit.getQueryIdx()+hit.getQuery().length();
		}
		if(nextSubcluster.size()>0) answer.add(nextSubcluster);
		return answer;
	}
	private static int getAverageEstimatedSubjectStart(List<UngappedSearchHit> clusterHits) {
		int average = 0;
		for(UngappedSearchHit hit:clusterHits) {
			average+= estimateSubjectStart(hit);
			
		}
		return average/clusterHits.size();
	}
	public void completeMissingHits(Map<Integer,Long> subjectCodes, Map<Integer, Long> queryCodes) {
		
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
				int j=lastSubjectStart+1;
				for(int i=lastQueryStart+1;i<queryStart;i++) {
					Long code = queryCodes.get(i);
					if(code == null) continue;
					List<Integer> posList = subjectCodesPos.get(code);
					if(posList==null) continue;
					Integer selectedPos = null;
					if(posList.size()<10) {
						for(int k:posList) {
							if(k>=j && k<j+10) {
								selectedPos = k;
								break;
							}
						}
					} else {
						for(int k=j;k<j+10;k++) {
							Long codeS = subjectCodes.get(k);
							if(codeS!=null && codeS.longValue() == code.longValue()) {
								selectedPos = k;
								break;
							}
						}
					}
					if(selectedPos==null) {
						j++;
					} else {
						CharSequence kmer = new String(AbstractLimitedSequence.getSequence(code, kmerLength, DNASequence.EMPTY_DNA_SEQUENCE));
						UngappedSearchHit rescuedHit = new UngappedSearchHit(kmer, subjectIdx, selectedPos);
						rescuedHit.setQueryIdx(i);
						double weight = (lastHit.getWeight()+hit.getWeight())/2;
						rescuedHit.setWeight(weight);
						addHit(rescuedHit);
						j = selectedPos+1;
					}
				}
			}
			lastQueryStart = queryStart;
			lastSubjectStart = subjectStart;
			lastHit = hit;
		}
	}
	
}