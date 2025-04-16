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
package ngsep.sequences;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.main.ThreadPoolManager;
import ngsep.math.Distribution;
import ngsep.math.EntropyCalculator;
import ngsep.math.ShannonEntropyCalculator;

/**
 * @author Jorge Duitama
 */
public class ShortKmerCodesTable {
	
	private EntropyCalculator entropyCalculator = new ShannonEntropyCalculator(4);
	private static final long [] EMPTY_LONG_ARRAY = new long[0];
	
	private Logger log = Logger.getLogger(ShortKmerCodesTable.class.getName());
	
	private int maxHitsKmerCode = 0;
	
	private ShortKmerCodesSampler codesSampler;
	
	private int mode=1;
	private int kmerDistModeLocalSD=5;
	
	//Structures to implement the minimizers like codes hash table
	//Map with code as key and row of the sequencesByCodeTable as value
	private Map<Long,Integer> matrixRowMap;
	//Table with encoded entries. Each row corresponds to a single code
	private long [][] sequencesByCodeTable;
	//Actual lengths of the lists within the table
	private short [] sequencesByCodeTableColumnLengths;
	
	private long totalEntries = 0;
	private ThreadPoolManager poolAddKmerCodes;
	
	
	public ShortKmerCodesTable(ShortKmerCodesSampler sampler) {
		this(sampler,1000, false);
	}
	public ShortKmerCodesTable(ShortKmerCodesSampler sampler, int capacity, boolean useThreadToAddCodes) {
		this.codesSampler = sampler;
		initializeTable(capacity, useThreadToAddCodes);
	}
	private void initializeTable(int capacity, boolean useThreadToAddCodes) {
		matrixRowMap = new HashMap<Long, Integer>(capacity);
		sequencesByCodeTable = new long [capacity][1];
		for(int i=0;i<sequencesByCodeTable.length;i++) Arrays.fill(sequencesByCodeTable[i], 0);
		sequencesByCodeTableColumnLengths = new short [capacity];
		Arrays.fill(sequencesByCodeTableColumnLengths, (short)0);
		if(useThreadToAddCodes) {
			poolAddKmerCodes = new ThreadPoolManager(1, 100);
			//poolAddKmerCodes.setSecondsPerTask(5);
		}
	}

	public int getMode() {
		return mode;
	}
	public void setMode(int mode) {
		this.mode = mode;
	}
	public int getKmerDistModeLocalSD() {
		return kmerDistModeLocalSD;
	}
	public void setKmerDistModeLocalSD(int kmerDistModeLocalSD) {
		this.kmerDistModeLocalSD = kmerDistModeLocalSD;
	}
	
	//Hash table management methods
	private long[] lookupHits(long code) {
		Integer row = matrixRowMap.get(code);
		if(row==null) return EMPTY_LONG_ARRAY;
		//Select hits up to the real number of hits
		return Arrays.copyOf(sequencesByCodeTable[row], sequencesByCodeTableColumnLengths[row]);
	}
	public int size() {
		return matrixRowMap.size();
	}
	
	/**
	 * Calculates the number of times that the code has been observed
	 * @param code number to query
	 * @return int times that the given code has been observed
	 */
	public int getTotalHits(long code) {
		Integer row = matrixRowMap.get(code);
		if(row==null) return 0;
		return sequencesByCodeTableColumnLengths[row];
	}
	
	private void addCodeSequence (long code, List<Long> entries) {
		Integer row = matrixRowMap.get(code);
		if(row==null) {
			row = size();
			if(row ==Integer.MAX_VALUE) {
				log.warning("Reached maximum number of minimizers that can be saved "+row);
				return;
			}
			matrixRowMap.put(code, row);
			if(row==sequencesByCodeTable.length) resizeTable();
		}
		int currentCount = sequencesByCodeTableColumnLengths[row];
		int newCount = currentCount+entries.size(); 
		if (newCount<Short.MAX_VALUE && (maxHitsKmerCode==0 || newCount<maxHitsKmerCode)) {
			for (long entry:entries) addToTable(row, entry);
			totalEntries+=entries.size();
		} 
		//else System.out.println("Rejected codes for kmer: "+new String(DNASequence.getDNASequence(code, kmerLength))+" current count: "+currentCount+" new entries: "+entries.size());
	}
	private void resizeTable() {
		log.info("Resizing codes table. Current number of codes: "+size()+" current capacity: "+sequencesByCodeTable.length);
		int newCapacity =  2*sequencesByCodeTable.length;
		if(newCapacity<0) newCapacity = Integer.MAX_VALUE;
		sequencesByCodeTableColumnLengths = Arrays.copyOf(sequencesByCodeTableColumnLengths, newCapacity);
		long [][] newTable = new long [newCapacity][0];
		for(int i=0;i<newCapacity;i++) {
			if(i<sequencesByCodeTable.length) newTable[i] = sequencesByCodeTable[i];
			else newTable[i] = new long[1];
		}
		sequencesByCodeTable = newTable;
		log.info("Resized codes table. New capacity: "+sequencesByCodeTable.length);
	}
	private void addToTable(int row, long value) {
		long [] codeEntries = sequencesByCodeTable[row];
		int column = sequencesByCodeTableColumnLengths[row];
		if(column == codeEntries.length ) {
			//Resize entries
			sequencesByCodeTable[row] = Arrays.copyOf(codeEntries, 2*codeEntries.length);
		}
		sequencesByCodeTable[row][column] = value;
		sequencesByCodeTableColumnLengths[row]++;
	}
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	
	public int getMaxHitsKmerCode() {
		return maxHitsKmerCode;
	}
	public void setMaxHitsKmerCode(int maxHitsKmerCode) {
		this.maxHitsKmerCode = maxHitsKmerCode;
	}

	/**
	 * Adds selected codes of the given sequence to the table
	 * @param sequenceId Id of the sequence to add
	 * @param sequence to add
	 */
	public void addSequence (int sequenceId, CharSequence sequence) {
		long startTime = System.currentTimeMillis();
		int n = sequence.length();
		if(n>1000000) log.info("Starting Sequence "+sequenceId);
		String sequenceStr = sequence.toString();
		int step = 500000;
		//Map<Long, List<Long>> codesSeq = new HashMap<Long, List<Long>>();
		List<KmerCodesTableEntry> codesSeq = new ArrayList<KmerCodesTableEntry>();
		int totalMinimizers = 0;
		for (int start = 0;start < n;start+=step) {
			List<KmerCodesTableEntry> codeEntriesList = codesSampler.computeSequenceCodes(sequenceId, sequenceStr, start, Math.min(n, start+step));
			//if(sequenceId==0) System.out.println("Seq id: "+sequenceId+" start: "+start+" raw codes: "+codeEntriesList.size()+" time: "+(System.currentTimeMillis()-startTime)/1000);
			codesSeq.addAll(codeEntriesList);
			if(codesSeq.size()>=1000000) {
				addCodesSequenceTask(sequenceId, n, codesSeq);
				//addCodesSequence(sequenceId, n, codesSeq);
				totalMinimizers+=codesSeq.size();
				codesSeq.clear();
			}
		}
		addCodesSequenceTask(sequenceId, n, codesSeq);
		//addCodesSequence(sequenceId, n, codesSeq);
		totalMinimizers+=codesSeq.size();
		long time = (System.currentTimeMillis()-startTime)/1000;
		if(n>1000000) log.info("Sequence "+sequenceId+" length: "+n+" total minimizers: "+totalMinimizers+" pct: "+(100*totalMinimizers/sequence.length())+" time(s): "+time);
	}
	private synchronized void addCodesSequenceTask(int sequenceId, int n, List<KmerCodesTableEntry> codesSeq) {
		final List<KmerCodesTableEntry> codesSeqToAdd = new ArrayList<KmerCodesTableEntry>(codesSeq);
		if(poolAddKmerCodes!=null) {
			try {
				poolAddKmerCodes.queueTask(()->addCodesSequence(sequenceId, n, codesSeqToAdd));
			} catch (InterruptedException e) {
				e.printStackTrace();
				throw new RuntimeException("Concurrence error creating minimizers table",e);
			}
		} else addCodesSequence(sequenceId, n, codesSeqToAdd);
	}
	private void addCodesSequence(int sequenceId, int seqLen, List<KmerCodesTableEntry> codesSeq) {
		//long startTime = System.currentTimeMillis();
		for(KmerCodesTableEntry entry:codesSeq) {
			List<Long> a = new ArrayList<Long>(1);
			a.add(entry.encode());
			addCodeSequence (entry.getKmerCode(), a);
		}
		//long time = (System.currentTimeMillis()-startTime)/1000;
		//if(seqLen>1000000) log.info("Sequence "+sequenceId+" length: "+seqLen+" entries added: "+codesSeq.size()+" time(s): "+time);
		codesSeq.clear();
	}
	public void endAddingSequences() {
		try {
			poolAddKmerCodes.terminatePool();
		} catch (InterruptedException e) {
			e.printStackTrace();
			throw new RuntimeException("Concurrence error creating minimizers table",e);
		}
	}

	/**
	 * Calculates the hits of the given query
	 * @param queryIdx Id of the query
	 * @param query sequence
	 * @return Map<Integer,List<MinimizersTableEntry>> Sequences matching kmers of the given query indexed by subject
	 */
	public Map<Integer,List<UngappedSearchHit>> match (int queryIdx, CharSequence query) {
		KmerSearchResultsCompressedTable hitsCompressed = matchCompressed(queryIdx, query, -1);
		return hitsCompressed.getAllHits();
	}
	/**
	 * Calculates the hits of the given query
	 * @param queryIdx Id of the query
	 * @param query sequence
	 * @param maxSubjectIdx Maximum id of a subject to process. -1 to disable this filter
	 * @return KmerSearchResultsCompressedTable Object encoding the results of the query 
	 */
	public KmerSearchResultsCompressedTable matchCompressed (int queryIdx, CharSequence query, int maxSubjectIdx) {
		//long [] segmentCodes = KmersExtractor.extractDNAKmerCodes(query.toString(), codesSampler.getKmerLength(), 0, query.length());
		//List<KmerCodesTableEntry> selectedCodeEntries = codesSampler.computeSequenceCodes(queryIdx, 0, segmentCodes);
		List<KmerCodesTableEntry> selectedCodeEntries = codesSampler.computeSequenceCodes(queryIdx, query.toString(), 0, query.length());
		Map<Integer, Long> selectedCodes = new LinkedHashMap<>();
		for(KmerCodesTableEntry entry:selectedCodeEntries) {
			long code = entry.getKmerCode();
			int start = entry.getStart();
			selectedCodes.put(start,code);
		}
		return matchCompressed(queryIdx, query.length(), selectedCodes, maxSubjectIdx);
		
	}
	/**
	 * Calculates the hits of the query represented by the given codes
	 * @param queryIdx Id of the query
	 * @param queryLength Length of the query
	 * @param codes representing the sequence. The codes length must coincide with this kmerLength
	 * @param maxSubjectIdx
	 * @return
	 */
	public KmerSearchResultsCompressedTable matchCompressed (int queryIdx, int queryLength, Map<Integer, Long> codes, int maxSubjectIdx) {
		int idxDebug = -2;
		int kmerLength = codesSampler.getKmerLength();
		int limitHits = 200;
		if(mode>1) limitHits = Math.max(limitHits, 5*mode);
		if (queryIdx == idxDebug) System.out.println("ShortKmerCodesTable. Aligning a total of "+codes.size()+" codes. Mode: "+mode+" kmer length: "+kmerLength+" limit hits: "+limitHits);
		
		Set<Integer> preselectedSubjectIds = null;
		if(queryLength>20000) {
			preselectedSubjectIds = preselectSubjectIds(queryLength, codes);
		}
		
		KmerSearchResultsCompressedTable result = new KmerSearchResultsCompressedTable(codes, kmerLength, 10);
		Map<Long,Integer> internalMultiHitCodes = calculateInternalMultihitKmers(codes);
		
		int numUsedCodes = 0;
		int notFoundCodes = 0;
		int multiSequenceCodes = 0;
		int multihitCodes = 0;
		int [] normalizedCountsDist = new int [10];
		Set<Long> usedCodes = new HashSet<>();
		Set<Long> internalMultiUsedCodes = new HashSet<>();
		Set<Long> highDepthUsedCodes = new HashSet<>();
		int selfSequenceCount = 0;
		for(Map.Entry<Integer, Long> entry:codes.entrySet()) {
			int startQuery = entry.getKey();
			long kmerCode = entry.getValue();
			//if (queryIdx == idxDebug ) System.out.println("Minimizers table. For pos "+startQuery+" code: "+kmerCode+" kmer: "+new String (DNASequence.getDNASequence(kmerCode, kmerLength)));
			
			long [] codesMatching = lookupHits(kmerCode);
			int numHits = codesMatching.length;
			int normalized = (int)Math.round((double)numHits/mode);
			if(normalized<normalizedCountsDist.length) normalizedCountsDist[normalized]++;
			else normalizedCountsDist[normalizedCountsDist.length-1]++;
			int internalCount = internalMultiHitCodes.getOrDefault(kmerCode,1);
			if(normalized>=2) {
				highDepthUsedCodes.add(kmerCode);
				if(internalCount>1) internalMultiUsedCodes.add(kmerCode);
			}
			if (queryIdx == idxDebug) System.out.println("Minimizers table. For pos "+startQuery+" kmer: "+new String (DNASequence.getDNASequence(kmerCode, kmerLength))+" codes matching: "+codesMatching.length+" internal multihit: "+(internalMultiHitCodes.getOrDefault(kmerCode,0)));
			if(numHits>limitHits) {
				multihitCodes++;
				if(internalCount>1) {
					limitHits = Math.max(200,10*internalCount);
					if(mode>1) limitHits = Math.max(limitHits,2*internalCount*mode);
					limitHits = Math.min(1000,limitHits);
					if(codesMatching.length>limitHits) codesMatching = subsampleCodes(codesMatching,limitHits);
				}
				else codesMatching = subsampleCodes(codesMatching,limitHits);
			} 
			if(numHits>0) {
				usedCodes.add(kmerCode);	 
				numUsedCodes++;
			} else notFoundCodes++;
			
			for(long entryCode:codesMatching) {
				int [] dec = KmerCodesTableEntry.decode(entryCode);
				int subjectIdx = dec[0];
				if (queryIdx == idxDebug && subjectIdx==0) System.err.println("KmerCodesTable. For pos "+startQuery+" kmer: "+new String (DNASequence.getDNASequence(kmerCode, kmerLength))+" next match: "+subjectIdx+" start: "+dec[1]);
				if (subjectIdx < 0) {
					System.err.println("Invalid subject "+subjectIdx+" query code: "+kmerCode+" matching code: "+entryCode+" start: "+dec[1]);
					continue;
				}
				if(subjectIdx==queryIdx) selfSequenceCount++;
				if(preselectedSubjectIds==null || preselectedSubjectIds.contains(subjectIdx)) {
					if(maxSubjectIdx==-1 || subjectIdx<=maxSubjectIdx) result.addKmerHit(startQuery,subjectIdx,dec[1]);
				}
				
			}
		}
		result.setUsedCodesCount(numUsedCodes);
		result.setMultisequenceCodesCount(multiSequenceCodes);
		result.setMultihitCodesCount(multihitCodes);
		result.setNotFoundCodesCount(notFoundCodes);
		result.setDistinctUsedCodesCount(usedCodes.size());
		result.setNormalizedCountsDist(normalizedCountsDist);
		result.setKmerWeights(calculateCodeWeights(codes, kmerLength));
		result.setHighDepthUsedKmerCodes(highDepthUsedCodes);
		result.setInternalMultihitUsedKmerCodes(internalMultiUsedCodes);
		if (queryIdx == idxDebug) System.out.println("ShortKmerCodesTable. Total codes used: "+numUsedCodes+" not found: "+notFoundCodes+" self sequence count: "+selfSequenceCount+" multihit used: "+multihitCodes+" internal multihit total: "+internalMultiHitCodes.size()+" used: "+internalMultiUsedCodes.size()+" total hits recorded: "+result.getTotalHits());
		return result;
		
	}
	
	private long[] subsampleCodes(long[] codesMatching, int limitHits) {
		if(codesMatching.length<=limitHits) return codesMatching;
		int factor = codesMatching.length/limitHits + 1;
		int outLength = codesMatching.length/factor;
		List<Long> outputCodes = new ArrayList<Long>(outLength);
		for(int i=0;i<codesMatching.length;i++) {
			if(i%factor==0) outputCodes.add(codesMatching[i]);
		}
		int n = outputCodes.size();
		long [] answer = new long[n];
		for(int i=0;i<n;i++) answer[i] = outputCodes.get(i);
		return answer;
	}
	private Set<Integer> preselectSubjectIds(int queryLength, Map<Integer, Long> codes) {
		int minHits = queryLength/200;
		
		Map<Integer,Integer> subjectHitCounts = new HashMap<>();
		for(Map.Entry<Integer, Long> entry:codes.entrySet()) {
			long kmerCode = entry.getValue();
			//int count = codesLocalCounts.getOrDefault(kmerCode, 0);
			long [] codesMatching = lookupHits(kmerCode);
			for(long entryCode:codesMatching) {
				int [] dec = KmerCodesTableEntry.decode(entryCode);
				int subjectIdx = dec[0];
				subjectHitCounts.compute(subjectIdx, (k,v)-> ((v!=null)?v+1:1));
			}
		}
		Set<Integer> answer = new HashSet<>();
		for(Map.Entry<Integer, Integer> entry:subjectHitCounts.entrySet()) {
			if(entry.getValue()>=minHits) answer.add(entry.getKey());
		}
		return answer;
	}
	private Map<Long,Integer> calculateInternalMultihitKmers(Map<Integer, Long> searchCodes) {
		Map<Long,Integer> answer = new HashMap<>();
		Map<Long,Integer> countsMap = new HashMap<>();
		for(Map.Entry<Integer,Long> entry:searchCodes.entrySet()) {
			countsMap.compute(entry.getValue(), (k,v)->v!=null?v+1:1);
		}
		for(Map.Entry<Long,Integer> entry:countsMap.entrySet()) {
			if(entry.getValue()>2) {
				answer.put(entry.getKey(),entry.getValue());
			}
		}
		return answer;
	}
	public Map<Long,Double> calculateCodeWeights(Map<Integer, Long> codes, int length) {
		Map<Long,Double> answer = new HashMap<>();
		for(long code:codes.values()) {
			answer.computeIfAbsent(code, v->calculateWeight(code, length));
		}
		return answer;
	}

	public double calculateWeight(long code, int length) {
		if(mode > 1) {
			//if(kmersMap==null) return 1;
			//int countDifferent = getCountDifferentSequences(code);
			int totalCount = getTotalHits(code);
			/*int diff1 = countDifferent-mode;
			int diff2 = totalCount/countQuery-mode;
			if(diff1<=kmerDistModeLocalSD && diff2<=kmerDistModeLocalSD) return 1;
			int diff3=diff1+diff2-2*kmerDistModeLocalSD;
			if(diff3<1) diff3=1;*/
			int diff1 = totalCount-mode;
			if(diff1<=kmerDistModeLocalSD) return 1;
			int diff3=diff1-kmerDistModeLocalSD;
			return 1.0*mode/(mode+diff3);
		} 
		else {
			CharSequence sequence = new String (DNASequence.getDNASequence(code, length));
			double entropy = entropyCalculator.calculateEntropy(sequence);
			return entropyCalculator.normalizeEntropy(entropy);
		}
	}
	public Distribution calculateDistributionHits() {
		Distribution dist = new Distribution(1, 300, 1);
		int numCodes = size();
		for(int i=0;i<numCodes;i++) {
			dist.processDatapoint(sequencesByCodeTableColumnLengths[i]);	
		}
		return dist;
	}

	public long getTotalEntries() {
		return totalEntries;
	}
	
	
}
