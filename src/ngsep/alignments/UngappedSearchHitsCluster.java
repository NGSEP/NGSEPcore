/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.alignments;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.sequences.UngappedSearchHit;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class UngappedSearchHitsCluster {

	private int queryLength;
	//private int kmerLength;
	//Map indexed by query start
	private Map<Integer,UngappedSearchHit> hitsMap=new TreeMap<Integer, UngappedSearchHit>();
	private int subjectIdx;
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
	
	public UngappedSearchHitsCluster(int queryLength, int subjectIdx, int subjectLength, List<UngappedSearchHit> hits) {
		
		
		if(hits.size()==0) throw new RuntimeException("Invalid empty input hits for cluster"+hits.size());
		
		this.queryLength = queryLength;
		this.subjectIdx = subjectIdx;
		this.subjectLength = subjectLength;
		//Create cluster with selected hits
		boolean initialized = false;
		for(UngappedSearchHit hit:hits) {
			if(!initialized) {
				subjectPredictedStart = estimateSubjectStart(hit);
				subjectPredictedEnd = estimateSubjectEnd(hit);
				subjectEvidenceStart = hit.getSubjectStart();
				subjectEvidenceEnd = subjectEvidenceStart+hit.getHitLength();
				queryPredictedStart = estimateQueryStart(hit);
				queryPredictedEnd = estimateQueryEnd(hit);
				queryEvidenceStart = hit.getQueryStart();
				queryEvidenceEnd = hit.getQueryStart() + hit.getHitLength();
				initialized = true;
			}
			addHit(hit);
			//if(hits.size()==534) System.out.println("UngappedSeaarchHitsCluster. Added hit: "+hit.getQueryStart()+" "+hit.getSubjectStart()+" subjectPredicted: "+estimateSubjectStart(hit)+" "+estimateSubjectEnd(hit));
		}
		summarize();
	}

	

	public UngappedSearchHitsCluster(int queryLength, int subjectIdx, int subjectLength, UngappedSearchHit kmerHit) {
		this.queryLength = queryLength;
		this.subjectIdx = subjectIdx;
		this.subjectLength = subjectLength;
		int kmerQueryStart = kmerHit.getQueryStart();
		subjectPredictedStart = estimateSubjectStart(kmerHit);
		subjectPredictedEnd = estimateSubjectEnd(kmerHit);
		subjectEvidenceStart = kmerHit.getSubjectStart();
		subjectEvidenceEnd = subjectEvidenceStart+kmerHit.getHitLength();
		queryPredictedStart = estimateQueryStart(kmerHit);
		queryPredictedEnd = estimateQueryEnd(kmerHit);
		queryEvidenceStart = kmerQueryStart;
		queryEvidenceEnd = kmerQueryStart + kmerHit.getHitLength();
		hitsMap.put(kmerQueryStart, kmerHit);
		firstKmerPresent = queryEvidenceStart == 0;
		lastKmerPresent = queryEvidenceEnd==queryLength;
	}
	
	private void addHit(UngappedSearchHit hit) {
		hitsMap.put(hit.getQueryStart(), hit);
		int estStart = estimateSubjectStart(hit);
		int estEnd = estimateSubjectEnd(hit);
		if(estStart!=subjectPredictedStart || estEnd!=subjectPredictedEnd) allConsistent = false;
		subjectPredictedStart = Math.min(subjectPredictedStart, estStart);
		subjectPredictedEnd = Math.max(subjectPredictedEnd, estEnd);
		subjectEvidenceStart = Math.min(subjectEvidenceStart, hit.getSubjectStart());
		subjectEvidenceEnd = Math.max(subjectEvidenceEnd, hit.getSubjectStart()+hit.getHitLength());
		queryPredictedStart = Math.min(queryPredictedStart, estimateQueryStart(hit));
		queryPredictedEnd = Math.max(queryPredictedEnd, estimateQueryEnd(hit));
		queryEvidenceStart = Math.min(queryEvidenceStart, hit.getQueryStart());
		queryEvidenceEnd = Math.max(queryEvidenceEnd, hit.getQueryStart() + hit.getHitLength());
		if(queryEvidenceStart==0) firstKmerPresent = true;
		if (queryEvidenceEnd==queryLength) lastKmerPresent = true;
	}
	
	public boolean addKmerHit(UngappedSearchHit kmerHit, int toleranceChange) {
		int estStart = estimateSubjectStart(kmerHit);
		int estEnd = estimateSubjectEnd(kmerHit);
		//if (kmerHit.getSubjectStart()<100000) System.out.println("Hit with idx: "+kmerHit.getQueryStart()+" Previous coords: "+subjectPredictedStart+"-"+subjectPredictedEnd+" next cords: "+estStart+"-"+estEnd);
		if(subjectPredictedStart > estEnd || subjectPredictedEnd < estStart) return false;
		if(toleranceChange>0 && Math.abs(subjectPredictedStart-estStart)>toleranceChange) return false;
		if(toleranceChange>0 && Math.abs(subjectPredictedEnd-estEnd)>toleranceChange) return false;
		addHit(kmerHit);
		return true;
	}
	
	private int estimateSubjectStart(UngappedSearchHit hit) {
		return hit.getSubjectStart() - hit.getQueryStart();
	}

	private int estimateSubjectEnd(UngappedSearchHit hit) {
		return hit.getSubjectStart()+(queryLength-hit.getQueryStart());
	}
	
	private int estimateQueryStart(UngappedSearchHit hit) {
		return hit.getQueryStart() - hit.getSubjectStart();
	}

	private int estimateQueryEnd(UngappedSearchHit hit) {
		return hit.getQueryStart()+(subjectLength-hit.getSubjectStart());
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
		Collections.sort(hits, (h1,h2)->h1.getSubjectStart()-h2.getSubjectStart());
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
		//if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Predicting overlap. Avg: "+averagePredictedOverlap+" median: "+medianPredictedOverlap+" from limits: "+fromLimitsPredictedOverlap);
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
		//if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("KmerHitsCluster. evidenceStart "+queryEvidenceStart+" predicted start: "+queryPredictedStart+" sum: "+totalSumQuery+" weight sum: "+weightSum);
		//int d = queryPredictedStart-queryEvidenceStart;
		//if (n>100 && d>=50) System.out.println("WARN: Predicted start "+queryPredictedStart+" for cluster to subject: "+sequenceIdx+" much larger than evidence: "+queryEvidenceStart+" query length: "+query.length());
	}

	private void predictSubjectEnd(List<UngappedSearchHit> hits) {
		double totalSumSubject = 0;
		double weightSum = 0;
		int n = hits.size();
		for(int i=Math.max(0, n-50);i<n;i++) {
			UngappedSearchHit hit = hits.get(i);
			double weight = ((double)(i+1))/n;
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
			double weight = ((double)(i+1))/n;
			weightSum += weight;
			totalSumQuery += weight*estimateQueryEnd(hit);
			//if (sequenceIdx==idxSubjectDebug && query.length() == queryLengthDebug) System.out.println("Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" estqEnd: "+estimateQueryEnd(hit)+" weight: "+weight+" sum: "+totalSumQuery);
		}
		queryPredictedEnd = (int) Math.round(totalSumQuery / weightSum);
		//if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("KmerHitsCluster. evidenceEnd "+queryEvidenceEnd+" predicted end: "+queryPredictedEnd+" sum: "+totalSumQuery+" weight sum: "+weightSum);
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
	
	public void setRawKmerHits(int rawKmerHits) {
		this.rawKmerHits = rawKmerHits;
	}
	public double getRawKmerHitsSubjectStartSD() {
		return rawKmerHitsSubjectStartSD;
	}
	public void setRawKmerHitsSubjectStartSD(double rawKmerHitsSubjectStartSD) {
		this.rawKmerHitsSubjectStartSD = rawKmerHitsSubjectStartSD;
	}
	/**
	 * 
	 * @return the kmerNumbers
	 */
	public int getCountKmerHitsCluster() {
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
}