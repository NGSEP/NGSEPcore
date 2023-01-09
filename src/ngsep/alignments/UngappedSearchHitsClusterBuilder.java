package ngsep.alignments;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.math.Distribution;
import ngsep.sequences.UngappedSearchHit;

public class UngappedSearchHitsClusterBuilder {
	private int idxSubjectDebug = -1;
	private int queryLengthDebug = -1;
	private boolean debug = false;
	public List<UngappedSearchHitsCluster> clusterRegionKmerAlns(int queryLength, int subjectIdx, int subjectLength, List<UngappedSearchHit> sequenceHits, double minQueryCoverage) {
		Map<Integer,Integer> countsByQueryIdx = new HashMap<Integer, Integer>();
		double minHits = Math.min(20,0.01*queryLength);
		debug = subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug;
		if(sequenceHits.size()<minHits) return new ArrayList<UngappedSearchHitsCluster>();
		for(UngappedSearchHit hit:sequenceHits) {
			//if (debug && hit.getSubjectStart()<10000) System.out.println("Next qpos "+hit.getQueryStart()+" hit: "+hit.getSubjectStart());
			countsByQueryIdx.compute(hit.getQueryStart(), (k,v)->v==null?1:v+1);
		}
		//double estimatedClusters = 0.5*sequenceHits.size()/queryLength;
		double avg = 0;
		for(int count:countsByQueryIdx.values()) avg+=count;
		if(avg>0) avg/=countsByQueryIdx.size();
		double estimatedClusters = avg;
		 
		if(debug) System.out.println("Clustering hits: "+sequenceHits.size()+" estimatedCLusters: "+estimatedClusters);
		
		//if(minQueryCoverage==0) 
		return clusterRegionKmerAlnsMultiple(queryLength, subjectIdx, subjectLength, sequenceHits, estimatedClusters);
		//if(estimatedClusters>1.5) return clusterRegionKmerAlnsMultiple(queryLength, subjectLength, sequenceHits, estimatedClusters);
		/*List<UngappedSearchHitsCluster> answer = new ArrayList<>();
		List<UngappedSearchHit> selectedHits = collapseAndSelectSortedHits(queryLength, subjectLength, sequenceHits);
		UngappedSearchHitsCluster uniqueCluster = new UngappedSearchHitsCluster(queryLength, subjectLength, selectedHits);
		uniqueCluster.setRawKmerHits(sequenceHits.size());
		uniqueCluster.setRawKmerHitsSubjectStartSD(1);
		//if(querySequenceId==idxDebug) System.out.println("Hits to cluster: "+sequenceKmerHits.size()+" target: "+uniqueCluster.getSequenceIdx()+" first: "+uniqueCluster.getFirst()+" last: "+uniqueCluster.getLast()+" kmers: "+uniqueCluster.getNumDifferentKmers());
		if (uniqueCluster.getNumDifferentKmers()<minHits || uniqueCluster.getQueryEvidenceEnd()-uniqueCluster.getQueryEvidenceStart()<minQueryCoverage*queryLength) return answer;
		answer.add(uniqueCluster);
		if(uniqueCluster.getNumDifferentKmers()>0.8*sequenceHits.size()) return answer;
		//Cluster remaining hits
		List<UngappedSearchHit> remainingHits = new ArrayList<UngappedSearchHit>();
		for(UngappedSearchHit hit:sequenceHits) {
			if (hit!=uniqueCluster.getKmerHit(hit.getQueryIdx())) remainingHits.add(hit);
		}
		selectedHits = collapseAndSelectSortedHits(queryLength, subjectLength, remainingHits);
		UngappedSearchHitsCluster cluster2 = new UngappedSearchHitsCluster(queryLength, subjectLength, selectedHits);
		uniqueCluster.setRawKmerHits(remainingHits.size());
		uniqueCluster.setRawKmerHitsSubjectStartSD(1);
		if(cluster2.getNumDifferentKmers()>=minHits && cluster2.getQueryEvidenceEnd()-cluster2.getQueryEvidenceStart()>=minQueryCoverage*queryLength) answer.add(cluster2);
		return answer;*/
	}
	private List<UngappedSearchHitsCluster> clusterRegionKmerAlnsMultiple(int queryLength, int subjectIdx, int subjectLength, List<UngappedSearchHit> sequenceHits, double estimatedClusters) {
		double minHits = Math.min(20,0.01*queryLength);
		List<UngappedSearchHitsCluster> answer = new ArrayList<>();
		//UngappedSearchHit firstHit = sequenceHits.get(0);
		//int subjectIdx = firstHit.getSequenceIdx();
		//Initial clustering
		Map<Integer,List<UngappedSearchHit>> hitsByBin = new HashMap<Integer, List<UngappedSearchHit>>();
		for(UngappedSearchHit hit:sequenceHits) {
			int estStart = estimateSubjectStart(hit);
			int bin = estStart/1000;
			if(estStart<0) bin--;
			//if(hit.getQueryIdx()==10851) System.out.println("Initial hit clustering. Next hit start: "+hit.getStart()+" estimated: "+estStart+" cluster bin: "+bin);
			List<UngappedSearchHit> hitsBin = hitsByBin.computeIfAbsent(bin, (v)-> new ArrayList<UngappedSearchHit>());
			hitsBin.add(hit);
		}
		if(debug) System.out.println("Clustering kmer hits. Initial bin starts: "+hitsByBin.keySet()+" estimated number of clusters: "+estimatedClusters);
		//Second cluster centered in best averages
		List<List<UngappedSearchHit>> sortedClusters = new ArrayList<List<UngappedSearchHit>>(hitsByBin.size());
		sortedClusters.addAll(hitsByBin.values());
		Collections.sort(sortedClusters,(c1,c2)->c2.size()-c1.size());
		hitsByBin.clear();
		List<Integer> clusterAverages = new ArrayList<Integer>((int) (2*estimatedClusters+1));
		for(int i=0;i<sortedClusters.size() && i<=2*(estimatedClusters);i++) {
			List<UngappedSearchHit> cluster = sortedClusters.get(i);
			//if(cluster.size()<minHits) break;
			int average = getAverageEstimatedSubjectStart(cluster);
			if(debug) System.out.println("Clustering kmer hits. Estimated start First hit: "+estimateSubjectStart(cluster.get(0))+" Average predicted start next cluster: "+average+" size: "+cluster.size());
			clusterAverages.add(average);
		}
		if (clusterAverages.size()==0) return answer;
		Collections.sort(clusterAverages, (c1,c2)->c1-c2);
		int next = clusterAverages.get(0);
		for(int average:clusterAverages) {
			if(average-next<500) next = (next+average)/2;
			else {
				hitsByBin.put(next, new ArrayList<UngappedSearchHit>());
				if(debug) System.out.println("Clustering kmer hits. Next average start: "+next);
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
			//if(hit.getQueryIdx()==10851) System.out.println("Clustering kmer hits. Next hit: "+hit.getQueryIdx()+" estimated start "+estStart+" min s: "+minS+" minD: "+minD);
			//if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) 
			if(minD<500) hitsByBin.get(minS).add(hit);
		}
		if(debug) System.out.println("Clustering kmer hits. Final bin starts: "+hitsByBin.keySet());
		List<List<UngappedSearchHit>> hitsClusters=new ArrayList<List<UngappedSearchHit>>(hitsByBin.values());
		Collections.sort(hitsClusters,(l1,l2)->l2.size()-l1.size());
		for(List<UngappedSearchHit> hits:hitsClusters) {
			if(debug) System.out.println("Next candidate cluster size: "+hits.size()+" current limit: "+minHits);
			if(hits.size()<minHits) break;
			List<List<UngappedSearchHit>> subclusters = breakByQueryStart(hits);
			for(List<UngappedSearchHit> subcluster:subclusters) {
				if(subcluster.size()<minHits) continue;
				List<UngappedSearchHit> selectedHits = collapseAndSelectSortedHits(queryLength, subjectIdx, subjectLength, subcluster);
				if(selectedHits.size()<minHits) continue;
				UngappedSearchHitsCluster cluster = new UngappedSearchHitsCluster(queryLength, subjectIdx, subjectLength, selectedHits);
				cluster.setRawKmerHits(sequenceHits.size());
				cluster.setRawKmerHitsSubjectStartSD(1);
				if(debug) System.out.println("Next cluster subject predicted coords: "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" subject evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" query evidence: "+cluster.getQueryEvidenceStart()+" "+cluster.getQueryEvidenceEnd()+" unique kmers: "+cluster.getNumDifferentKmers());
				if(cluster.getNumDifferentKmers()>=minHits) answer.add(cluster);
				minHits = Math.max(minHits, cluster.getNumDifferentKmers()/2);
			}
		}
		return answer;
	}
	private List<List<UngappedSearchHit>> breakByQueryStart(List<UngappedSearchHit> hits) {
		List<List<UngappedSearchHit>> answer = new ArrayList<List<UngappedSearchHit>>();
		Collections.sort(hits, (h1,h2)->h1.getQueryStart()-h2.getQueryStart());
		List<UngappedSearchHit> nextSubcluster = new ArrayList<UngappedSearchHit>();
		int lastQueryStart = 0;
		for(UngappedSearchHit hit:hits) {
			//TODO: define better this parameter
			if(hit.getQueryStart()>lastQueryStart+30000) {
				if(nextSubcluster.size()>0) answer.add(nextSubcluster); 
				nextSubcluster = new ArrayList<UngappedSearchHit>();
			}
			nextSubcluster.add(hit);
			lastQueryStart = hit.getQueryStart();
		}
		if(nextSubcluster.size()>0) answer.add(nextSubcluster);
		return answer;
	}
	private int getAverageEstimatedSubjectStart(List<UngappedSearchHit> clusterHits) {
		long average = 0;
		for(UngappedSearchHit hit:clusterHits) {
			average+= estimateSubjectStart(hit);
			
		}
		return (int) (average/clusterHits.size());
	}
	private int estimateSubjectStart(UngappedSearchHit hit) {
		return hit.getSubjectStart() - hit.getQueryStart();
	}
	private int estimateQueryStart(UngappedSearchHit hit) {
		return hit.getQueryStart() - hit.getSubjectStart();
	}
	/**
	 * 
	 * @param queryLength
	 * @param subjectLength
	 * @param hits sorted by queryidx
	 * @return
	 */
	private List<UngappedSearchHit> collapseAndSelectSortedHits(int queryLength, int subjectIdx, int subjectLength, List<UngappedSearchHit> inputHits) {
		
		if(debug) System.out.println("KmerHitsCluster. Clustering "+inputHits.size()+" hits. Subject idx: "+subjectIdx);
		
		//Index hits by query kmer start
		Map<Integer,List<UngappedSearchHit>> hitsMultiMap = new TreeMap<Integer, List<UngappedSearchHit>>();
		for(UngappedSearchHit hit:inputHits) {
			//if (subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" kmer: "+hit.getQuery()+ " estq: "+estimateQueryStart(hit)+" - "+estimateQueryEnd(hit)+" estS: "+estimateSubjectStart(hit)+" - "+estimateSubjectEnd(hit));
			List<UngappedSearchHit> list = hitsMultiMap.computeIfAbsent(hit.getQueryStart(), l -> new ArrayList<UngappedSearchHit>());
			list.add(hit);
		}
		if(debug) System.out.println("KmerHitsCluster. Num different kmers: "+hitsMultiMap.size());
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
		double rawKmerHitsSubjectStartSD = (variance>0)?Math.sqrt(variance):1;
		Distribution dist = new Distribution(-400, 400, 50);
		Distribution distAbs = new Distribution(0, 500, 1);
		for(int start:subjectStarts) {
			int distance = start-median;
			int distanceAbs = Math.abs(distance);
			if (distanceAbs < 2*rawKmerHitsSubjectStartSD) {
				dist.processDatapoint(distance);
				distAbs.processDatapoint(distanceAbs);
			}
		}
		int modeDist = (int) Math.round(dist.getLocalMode(-300, 300));
		if(Math.abs(modeDist)>=100) median+=(modeDist/2);
		
		//System.out.println(subjectStarts);
		if(debug) dist.printDistributionInt(System.out);
		int maxDistance = (int) Math.max(distAbs.getAverage(), Math.sqrt(distAbs.getVariance()+1));
		maxDistance *=5;
		maxDistance=Math.min(queryLength/20,maxDistance);
		maxDistance=Math.max(queryLength/40,maxDistance);
		maxDistance=Math.max(10,maxDistance);
		//if(maxDistance<0.01*query.length()) maxDistance*=2;
		if(debug) System.out.println("KmerHitsCluster. Num hits: "+n+" median: "+median+" average: "+(sum/n)+" variance: "+variance+" stdev: "+rawKmerHitsSubjectStartSD+" abs distance avg: "+distAbs.getAverage()+" stdev "+Math.sqrt(distAbs.getVariance())+" max distance: "+maxDistance);
		
		List<UngappedSearchHit> selectedHits = selectHitsByDistanceWithMedian(hitsMultiMap, median, maxDistance);
		if(selectedHits.size()<1) { 
			if(debug) System.err.println("WARN. Empty list of selected hits for subject: "+subjectIdx);
			return selectedHits;
		}
		replaceHitsByLocalAgreement(selectedHits, hitsMultiMap, median, maxDistance, queryLength);
		List<UngappedSearchHit> filteredHits = removeDisorganized (selectedHits, median);
		if(filteredHits.size()<1) {
			if(debug) System.err.println("WARN. Empty list of sorted hits for subject: "+subjectIdx+" selected hits: "+selectedHits.size()+" query length: "+queryLength);
			
		}
		return filteredHits;
	}
	private List<UngappedSearchHit> selectHitsByDistanceWithMedian(Map<Integer, List<UngappedSearchHit>> hitsMultiMap, int median, int maxDistance) {
		List<UngappedSearchHit> selectedHits = new ArrayList<UngappedSearchHit>(hitsMultiMap.size());
		for(List<UngappedSearchHit> hits:hitsMultiMap.values()) {
			UngappedSearchHit hit = selectHit(hits, median, maxDistance);
			if(hit!=null) {
				if (debug) System.out.println("Selected hits. Next qpos "+hit.getQueryStart()+" hit: "+hit.getSubjectStart()+" estq: "+estimateQueryStart(hit)+" estS: "+estimateSubjectStart(hit)+" all starts: "+calculateHitStarts(hits));
				selectedHits.add(hit);
			}
		}
		return selectedHits;
	}

	private List<Integer> calculateHitStarts(List<UngappedSearchHit> hits) {
		List<Integer> hitStarts = new ArrayList<Integer>();
		for(UngappedSearchHit hit2:hits) hitStarts.add(estimateSubjectStart(hit2));
		return hitStarts;
	}
	
	private UngappedSearchHit selectHit(List<UngappedSearchHit> hits, int median, int maxDistance) {
		UngappedSearchHit answer = null;
		int minCost = 0;
		
		for(UngappedSearchHit hit:hits) {
			int estStart = estimateSubjectStart(hit);
			int distance = Math.abs(estStart-median);
			int cost = distance;
			//if(debug && hit.getQueryIdx()==5314) System.out.println("Selecting hit. Next hit: "+hit.getStart()+" distance: "+distance+" max "+maxDistance+" cost: "+cost);
			if(distance <= maxDistance) {
				if(answer == null || minCost>cost) {
					answer = hit;
					minCost = cost;
				}
			}
		}
		return answer;
	}
	
	private void replaceHitsByLocalAgreement(List<UngappedSearchHit> selectedHits, Map<Integer, List<UngappedSearchHit>> hitsMultiMap, int median, int maxDistance, int queryLength) {
		//Find trustable site
		int minHitPos = -1;
		int minCost = -1;
		int midPoint = queryLength/2;
		for(int i=0;i<selectedHits.size();i++) {
			UngappedSearchHit hit = selectedHits.get(i);
			int estStart = estimateSubjectStart(hit);
			int distance = Math.abs(estStart-median);
			int d2 = Math.abs(hit.getQueryStart()-midPoint);
			int cost = distance + d2/10;
			if(minHitPos==-1 || cost<minCost) {
				minHitPos = i;
				minCost = cost;
			}
		}
		if(minHitPos==-1) return;
		UngappedSearchHit minHit = selectedHits.get(minHitPos);
		if (debug) System.out.println("Hit closest to median. qpos "+minHit.getQueryStart()+" hit: "+minHit.getSubjectStart()+" cost: "+minCost+" estq: "+estimateQueryStart(minHit)+" estS: "+estimateSubjectStart(minHit));
		int vicinityEstStart = estimateSubjectStart(minHit);
		LinkedList<Integer> localValues = new LinkedList<Integer>();
		localValues.add(vicinityEstStart);
		for(int i=minHitPos-1;i>=0;i--) {
			UngappedSearchHit nextHit = selectedHits.get(i);
			int estStart = estimateSubjectStart(nextHit);
			if(Math.abs(estStart-vicinityEstStart)> 5) {
				List<UngappedSearchHit> hits = hitsMultiMap.get(nextHit.getQueryStart());
				vicinityEstStart = getMedian(localValues);
				UngappedSearchHit updatedHit = selectHit(hits, vicinityEstStart, 10);
				if (debug) System.out.println("Trying to replace outlier at qpos "+nextHit.getQueryStart()+" hit: "+nextHit.getSubjectStart()+" ests: "+estStart+" local median: "+vicinityEstStart+" all starts: "+calculateHitStarts(hits));
				if(updatedHit!=null && updatedHit!=nextHit) {
					if (debug) System.out.println("Replacing hit. qpos "+updatedHit.getQueryStart()+" hit: "+updatedHit.getSubjectStart()+" estq: "+estimateQueryStart(updatedHit)+" estS: "+estimateSubjectStart(updatedHit));
					selectedHits.set(i, updatedHit);
					estStart = estimateSubjectStart(updatedHit);
				}
			} else {
				vicinityEstStart = estStart;
			}
			localValues.add(estStart);
			if(localValues.size()>10) localValues.removeFirst();
		}
		vicinityEstStart = estimateSubjectStart(minHit);
		localValues.clear();
		localValues.add(vicinityEstStart);
		for(int i=minHitPos+1;i<selectedHits.size();i++) {
			UngappedSearchHit nextHit = selectedHits.get(i);
			int estStart = estimateSubjectStart(nextHit);
			if(Math.abs(estStart-vicinityEstStart)> 5) {
				List<UngappedSearchHit> hits = hitsMultiMap.get(nextHit.getQueryStart());
				vicinityEstStart = getMedian(localValues);
				UngappedSearchHit updatedHit = selectHit(hits, vicinityEstStart, 10);
				if (debug) System.out.println("Trying to replace outlier at qpos "+nextHit.getQueryStart()+" hit: "+nextHit.getSubjectStart()+" ests: "+estStart+" local median: "+vicinityEstStart+" all starts: "+calculateHitStarts(hits));
				if(updatedHit!=null && updatedHit!=nextHit) {
					if (debug) System.out.println("Replacing hit. qpos "+updatedHit.getQueryStart()+" hit: "+updatedHit.getSubjectStart()+" estq: "+estimateQueryStart(updatedHit)+" estS: "+estimateSubjectStart(updatedHit));
					selectedHits.set(i, updatedHit);
					estStart = estimateSubjectStart(updatedHit);
				}
			} else {
				vicinityEstStart = estStart;
			}
			localValues.add(estStart);
			if(localValues.size()>10) localValues.removeFirst();
		}
		
	}

	private int getMedian(LinkedList<Integer> localValues) {
		ArrayList<Integer> copy = new ArrayList<Integer>(localValues);
		Collections.sort(copy);
		return copy.get(copy.size()/2);
	}

	private List<UngappedSearchHit> removeDisorganized(List<UngappedSearchHit> selectedHits, int median) {
		int n = selectedHits.size();
		//if(n<=30) return selectSorted(selectedHits, 0, n-1, median);
		boolean [] vicinityConsistent = new boolean [n];
		Arrays.fill(vicinityConsistent, false);
		for(int i=10;i<n-10;i++) {
			UngappedSearchHit hit = selectedHits.get(i);
			int start = hit.getSubjectStart();
			boolean consistent = true;
			for(int j=i-1;j>=i-10 && consistent;j--) {
				UngappedSearchHit hit2 = selectedHits.get(j);
				if(hit2.getSubjectStart()>=start) consistent = false; 
			}
			for(int j=i+1;j<=i+10 && consistent;j++) {
				UngappedSearchHit hit2 = selectedHits.get(j);
				if(hit2.getSubjectStart()<=start) consistent = false; 
			}
			vicinityConsistent[i] = consistent;
			if(!consistent) {
				for(int j=i-1;j>=i-10;j--) vicinityConsistent[j] = false;
				for(int j=i+1;j<=i+10 ;j++) vicinityConsistent[j] = false;
				i+=10;
			}
		}
		List<UngappedSearchHit> dpSortedHits = new ArrayList<UngappedSearchHit>();
		double sum = 0;
		double sum2 = 0;
		double nS = 0;
		int i=0;
		while(i<n) {
			UngappedSearchHit hit = selectedHits.get(i);
			if(vicinityConsistent[i]) {
				//if(debug) System.out.println("Adding consistent hit at "+hit.getQueryIdx());
				dpSortedHits.add(hit);
				double s = estimateSubjectStart(hit);
				sum+=s;
				sum2+=s*s;
				nS++;
				i++;
			} else {
				int j;
				for(j=i+1;j<n && !vicinityConsistent[j];j++);
				dpSortedHits.addAll(selectSorted(selectedHits,i,j-1,median));
				i=j;
			}
		}
		//Remove final outliers
		int n2= dpSortedHits.size();
		if(debug) System.out.println("Initial hits: "+n+" consistent: "+nS+" dpSorted: "+n2);
		if(n2<30 || nS < 20) return dpSortedHits;
		int averageStart = (int) (sum/nS);
		double variance = Math.max(1, (sum2-sum*sum/nS)/(nS-1));
		double stdev = Math.sqrt(variance);
		double maxDistance = Math.max(100, 3*stdev);
		int countOutliers =0;
		for(UngappedSearchHit hit:dpSortedHits) {
			double distance = Math.abs(estimateSubjectStart(hit)-averageStart);
			if(distance > maxDistance) countOutliers++;
		}
		
		if(debug) System.out.println("Average start: "+averageStart+" stdev: "+stdev+" max distance: "+maxDistance+" outliers: "+countOutliers+" n "+n2+" nconsistent: "+nS);
		List<UngappedSearchHit> answer = new ArrayList<UngappedSearchHit>(n2);
		if(countOutliers<0.02*n2) {
			for(UngappedSearchHit hit:dpSortedHits) {
				double distance = Math.abs(estimateSubjectStart(hit)-averageStart);
				if(distance <= maxDistance) answer.add(hit);
				else if(debug) System.out.println("Removing outlier: "+hit.getQueryStart()+" subject: "+hit.getSubjectStart()+" estimated start: "+estimateSubjectStart(hit));
			}
			return answer;
		}
		Collections.sort(dpSortedHits,(h1,h2)->Math.abs(estimateSubjectStart(h2)-averageStart)-Math.abs(estimateSubjectStart(h1)-averageStart));
		//if(debug) System.out.println("Sorted by distance");
		for (int k=0;k<dpSortedHits.size();k++) {
			UngappedSearchHit hit = dpSortedHits.get(k);
			double distance = Math.abs(estimateSubjectStart(hit)-averageStart);
			if(distance <= maxDistance || k>=0.02*n2) answer.add(hit);
			else if(debug) System.out.println("Removing outlier after sorting: "+hit.getQueryStart()+" subject: "+hit.getSubjectStart()+" estimated start: "+estimateSubjectStart(hit));
		}
		Collections.sort(answer,(h1,h2)->h1.getQueryStart()-h2.getQueryStart());
		return answer;
	}


	private List<UngappedSearchHit> selectSorted(List<UngappedSearchHit> selectedHits, int first, int last, int median) {
		if(debug) System.out.println("DP Sorting hits from "+first+" to "+last+" median start: "+median);
		//Check sorted first
		List<UngappedSearchHit> regionHits = new ArrayList<UngappedSearchHit>(last-first+1);
		int lastStartSubject = -1;
		int lastStartQuery = -1;
		boolean sortedSubject = true;
		boolean sortedQuery = true;
		for(int i=first;i<=last;i++) {
			UngappedSearchHit hit = selectedHits.get(i);
			regionHits.add(hit);
			int startS = hit.getSubjectStart();
			if(lastStartSubject>=startS) sortedSubject = false;
			lastStartSubject = startS;
			int startQ = hit.getQueryStart();
			if(lastStartQuery>=startQ) sortedQuery = false;
			lastStartQuery = startQ;
		}
		if(!sortedQuery) {
			System.err.println("WARN. Disorganized query hits list");
			Collections.sort(regionHits,(h1,h2)->h1.getQueryStart()-h2.getQueryStart());
			
		} else if (sortedSubject) return regionHits;
		List<UngappedSearchHit> regionHitsBySubject = new ArrayList<UngappedSearchHit>();
		regionHitsBySubject.addAll(regionHits);
		Collections.sort(regionHitsBySubject,(h1,h2)->h1.getSubjectStart()-h2.getSubjectStart());
		int n = regionHits.size();
		//Dynamic programming to select the best organized subset
		int [][] scores = new int [n+1][n+1];
		for(int i=0;i<=n;i++) {
			for(int j=0;j<=n;j++) {
				if(i==0 || j==0) scores[i][j] = 0;
				else {
					UngappedSearchHit hit1 = regionHits.get(i-1);
					UngappedSearchHit hit2 = regionHitsBySubject.get(j-1);
					if(hit1==hit2) scores[i][j] = calculateScore(hit1,median) + scores[i-1][j-1];
					else scores[i][j] = Math.max(scores[i-1][j], scores[i][j-1]);
				}
			}
			 
		}
		//if(debug) printMatrix(scores);
		List<UngappedSearchHit> answer = new ArrayList<UngappedSearchHit>();
		int i= n;
		int j = n;
		while(i>0 && j>0) {
			UngappedSearchHit hit1 = regionHits.get(i-1);
			UngappedSearchHit hit2 = regionHitsBySubject.get(j-1);
			if(hit1==hit2) {
				//if(debug) System.out.println("Added DP sorted hit at "+hit1.getQueryIdx());
				answer.add(hit1);
				i--;
				j--;
			} else if (scores[i][j] == scores[i-1][j]) i--;
			else j--;
		}
		Collections.reverse(answer);
		return answer;
	}
	public void printMatrix(int[][] scores) {
		for(int i=0;i<scores.length;i++) {
			for(int j=0;j<scores[i].length;j++) {
				if(j>0) System.out.print("\t");
				System.out.print(""+scores[i][j]);
			}
			System.out.println();
		}
		
	}
	private int calculateScore(UngappedSearchHit hit, int median) {
		/*int start = estimateSubjectStart(hit);
		int distance = Math.abs(start-median);
		int score = (int) (1000*Math.max(0, 1-((double)distance)/1000.0));
		if(debug) System.out.println("Score for hit: "+hit.getQueryIdx()+" "+score+" start. "+start+" distance: "+distance);
		return score;*/
		return 1;
	}
	
}
