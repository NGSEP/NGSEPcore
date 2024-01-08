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
import ngsep.math.PrimeNumbers;

public class ShortKmerCodesTable {
	
	private static final long [] EMPTY_LONG_ARRAY = new long[0];
	
	private Logger log = Logger.getLogger(ShortKmerCodesTable.class.getName());
	
	private int kmerLength;
	private int windowLength;
	private int maxHitsKmerCode = 0;
	
	
	private KmersMapAnalyzer kmersAnalyzer;
	private KmersMap kmersMap;
	private int mode=1;
	private int kmerDistModeLocalSD=5;
	private int limitHitsPerSequence = 10;
	
	//Structures to implement the minimizers like codes hash table
	//Map with code as key and row of the sequencesByCodeTable as value
	private Map<Long,Integer> matrixRowMap;
	//Table with encoded entries. Each row corresponds to a single code
	private long [][] sequencesByCodeTable;
	//Actual lengths of the lists within the table
	private short [] sequencesByCodeTableColumnLengths;
	//Count of different sequences reporting each code
	private short [] codeCountDifferentSequences;
	
	private long totalEntries = 0;
	private PrimeNumbers primeNumbersHelper;
	private ThreadPoolManager poolAddKmerCodes;
	
	
	public ShortKmerCodesTable(int kmerLength, int windowLength) {
		this(kmerLength,windowLength,10, false);
	}
	public ShortKmerCodesTable(int kmerLength, int windowLength, int capacity, boolean useThreadToAddCodes) {
		this.kmerLength = kmerLength;
		this.windowLength = windowLength;
		this.primeNumbersHelper = new PrimeNumbers();
		initializeTable(capacity, useThreadToAddCodes);
	}
	public ShortKmerCodesTable(KmersMapAnalyzer kmersAnalyzer, int kmerLength, int windowLength) {
		this.kmersAnalyzer = kmersAnalyzer;
		this.kmersMap = kmersAnalyzer.getKmersMap();
		this.kmerLength = kmerLength;
		this.windowLength = windowLength;
		this.mode = Math.max(1,kmersAnalyzer.getMode());
		this.kmerDistModeLocalSD = kmersAnalyzer.getModeLocalSD();
		//TODO: This is ok for assemblies but it is not ok for reference genomes
		int capacityPrimes = (int) Math.min(100000000, kmersAnalyzer.getNumKmers(this.mode));
		this.primeNumbersHelper = new PrimeNumbers(capacityPrimes);
		
		//Create the structures with appropriate initial capacity
		int capacity = kmersMap.size()/10;
		initializeTable(capacity,false);
	}
	private void initializeTable(int capacity, boolean useThreadToAddCodes) {
		matrixRowMap = new HashMap<Long, Integer>(capacity);
		sequencesByCodeTable = new long [capacity][1];
		for(int i=0;i<sequencesByCodeTable.length;i++) Arrays.fill(sequencesByCodeTable[i], 0);
		sequencesByCodeTableColumnLengths = new short [capacity];
		Arrays.fill(sequencesByCodeTableColumnLengths, (short)0);
		codeCountDifferentSequences = new short [capacity];
		Arrays.fill(codeCountDifferentSequences, (short)0);
		if(useThreadToAddCodes) poolAddKmerCodes = new ThreadPoolManager(1, 100);
	}
	
	//Hash table management methods
	private long[] lookupHits(long code) {
		Integer row = matrixRowMap.get(code);
		if(row==null) return EMPTY_LONG_ARRAY;
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
	
	/**
	 * Calculates the number of sequences where the minimizer has been observed
	 * @param minimizer number to query
	 * @return int number of different sequences where the minimizer has been observed
	 */
	public int getCountDifferentSequences(long code) {
		Integer row = matrixRowMap.get(code);
		if(row==null) return 0;
		return codeCountDifferentSequences[row];
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
			codeCountDifferentSequences[row]++;
		} 
		//else System.out.println("Rejected codes for kmer: "+new String(DNASequence.getDNASequence(code, kmerLength))+" current count: "+currentCount+" new entries: "+entries.size());
	}
	private void resizeTable() {
		log.info("Resizing codes table. Current number of codes: "+size()+" current capacity: "+sequencesByCodeTable.length);
		int newCapacity =  2*sequencesByCodeTable.length;
		if(newCapacity<0) newCapacity = Integer.MAX_VALUE;
		sequencesByCodeTableColumnLengths = Arrays.copyOf(sequencesByCodeTableColumnLengths, newCapacity);
		codeCountDifferentSequences = Arrays.copyOf(codeCountDifferentSequences, newCapacity);
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
	
	public int getKmerLength() {
		return kmerLength;
	}
	public KmersMap getKmersMap() {
		return kmersMap;
	}
	
	public void setKmersMap(KmersMap kmersMap) {
		this.kmersMap = kmersMap;
	}
	
	
	public int getMaxHitsKmerCode() {
		return maxHitsKmerCode;
	}
	public void setMaxHitsKmerCode(int maxHitsKmerCode) {
		this.maxHitsKmerCode = maxHitsKmerCode;
	}

	public int getLimitHitsPerSequence() {
		return limitHitsPerSequence;
	}
	public void setLimitHitsPerSequence(int limitHitsPerSequence) {
		this.limitHitsPerSequence = limitHitsPerSequence;
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
			List<KmerCodesTableEntry> codeEntriesList = computeSequenceCodes(sequenceId, sequenceStr, start, Math.min(n, start+step));
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
	private void addCodesSequenceTask(int sequenceId, int n, List<KmerCodesTableEntry> codesSeq) {
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
	private synchronized void addCodesSequence(int sequenceId, int seqLen, List<KmerCodesTableEntry> codesSeq) {
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
	 * Calculates the selected codes of the given sequence following the same algorithm used for minimizers but saving the codes instead of the hashes
	 * @param sequenceId Id of the sequence to calculate
	 * @param sequence characters of the sequence to calculate
	 * @return List<MinimizersTableEntry> Codes selected for the given sequence.
	 */
	public Map<Integer,Long> computeSequenceCodesAsMap(String sequence,int start,int end) {
		List<KmerCodesTableEntry> entries = computeSequenceCodes(-1, sequence, start , end);
		Map<Integer,Long> answer = new LinkedHashMap<>();
		for(KmerCodesTableEntry entry:entries) {
			answer.put(entry.getStart(), entry.getKmerCode());
		}
		return answer;
	}
	public Map<Integer,Long> computeSequenceCodesAsMap(String sequence,int start,int end, long [] kmerCodes) {
		List<KmerCodesTableEntry> entries = computeSequenceCodes(-1, start, Math.min(end, sequence.length()-kmerLength-windowLength), kmerCodes);
		Map<Integer,Long> answer = new LinkedHashMap<>();
		for(KmerCodesTableEntry entry:entries) {
			answer.put(entry.getStart(), entry.getKmerCode());
		}
		return answer;
	}
	/**
	 * Calculates the selected codes of the given sequence following the same algorithm used for minimizers but saving the codes instead of the hashes
	 * @param sequenceId Id of the sequence to calculate
	 * @param sequence characters of the sequence to calculate
	 * @return List<MinimizersTableEntry> Codes selected for the given sequence.
	 */
	public List<KmerCodesTableEntry> computeSequenceCodes(int sequenceId, String sequence,int start,int end) {
		//new PrimeNumbers(1000000);
		long [] segmentCodes = KmersExtractor.extractDNAKmerCodes(sequence, kmerLength, start, Math.min(sequence.length(),end+windowLength+kmerLength));
		//log.info("Extracted codes for sequence "+sequenceId+" from "+start+" to "+end+" Number of codes: "+codes.size());
		//return new ArrayList<KmerCodesTableEntry>();
		return computeSequenceCodes(sequenceId, start, Math.min(end, sequence.length()-kmerLength-windowLength), segmentCodes);
	}
	/**
	 * Calculates the selected codes of the given codes following the same algorithm used for minimizers but saving the codes instead of the hashes
	 * @param sequenceId Id of the sequence to calculate
	 * @param start of the sequence to consider
	 * @param end of the sequence to consider
	 * @param kmerCodes Input codes to be selected
	 * @return List<KmerCodesTableEntry> selected codes
	 */
	private List<KmerCodesTableEntry> computeSequenceCodes(int sequenceId, int start, int end, long [] segmentCodes) {
		int debugIdx = -2;
		List<KmerCodesTableEntry> answer = new ArrayList<KmerCodesTableEntry>();
		Integer [] hashcodes = new Integer [segmentCodes.length];
		Arrays.fill(hashcodes, null);
		for(int i=0;i<segmentCodes.length;i++) {
			long code = segmentCodes[i];
			if(code<0) continue;
			hashcodes[i]=getHash(code);
		}
		if(sequenceId==debugIdx) System.err.println("Calculated hash codes for sequence "+sequenceId+" from "+start+" to "+end);
		Integer previousMinimizer = null;
		int previousMinimizerPos = -1;
		for(int i=start;i<end;i++) {
			int posIArray = i-start;
			int posLastKmerWindow = posIArray+windowLength-1;
			if(posLastKmerWindow>=segmentCodes.length) break;
			Integer minimizerI = null;
			int minPos = -1;
			Integer newHash = hashcodes[posLastKmerWindow];
			boolean lastInRange = previousMinimizer!=null && previousMinimizerPos>=i;
			if(lastInRange && (newHash==null || previousMinimizer < newHash)) {
				minimizerI = previousMinimizer;
				minPos = previousMinimizerPos;
			} else if (newHash!=null && (previousMinimizer==null || newHash <= previousMinimizer)) {
				minimizerI = newHash;
				minPos = i+windowLength-1;
			}
			if(minimizerI == null) {
				for(int j=0;j<windowLength;j++) {
					Integer hash = hashcodes[posIArray+j];
					if (hash!=null && (minimizerI==null || hash <= minimizerI)) {
						minimizerI = hash;
						minPos = i+j;
					}
				}
				//if(sequenceId==debugIdx && i>0 && i<3000) System.err.println("Minimizer calculated with cycle. Start: "+i+" New pos: "+minPos+" new minimizer: "+minimizerI+" previous: "+previousMinimizer+" kmer code: "+kmerCodes.get(minPos)+" total: "+minimizersSeq.size());
			}
			if (minimizerI==previousMinimizer) continue;
			if(minimizerI != null) {
				long originalCode = segmentCodes[minPos-start];
				KmerCodesTableEntry entry = new KmerCodesTableEntry(originalCode, sequenceId, minPos);
				answer.add(entry);
			}
			previousMinimizer = minimizerI;
			previousMinimizerPos = minPos;
		}
		if(sequenceId==debugIdx) System.err.println("Selected codes for sequence "+sequenceId+" from "+start+" to "+end+" Codes: "+answer.size());
		return answer;
	}

	private int getHash(long dnaHash) {

		long range = (long)(Math.pow(2, Long.toBinaryString(dnaHash).length())-1);
		int prime = 1073676287;
		//if(kmersAnalyzer!=null) {
		if(kmersAnalyzer==null) {

			long answer = (~dnaHash +(dnaHash << 21)) & range;
			answer = (answer ^ answer >> 24);
			answer = (answer + (answer << 3) + (answer << 8)) & range;
			answer = (answer ^ answer >> 14);
			answer = (answer + (answer << 2) + (answer << 4)) & range;
			answer = (answer ^ answer >> 28);
			answer = (answer + (answer << 31)) & range;

			return (int) answer;
		}
		int count;
		if(kmersMap instanceof ShortArrayDNAKmersMapImpl) {
			//Select the 15-mer suffix
			count = ((ShortArrayDNAKmersMapImpl)kmersMap).getCount(dnaHash & 0x3FFFFFFF);
			//count = ((ShortArrayDNAKmersMapImpl)kmersMap).getCount(dnaHash);
			//count = 20;
		} else {
			String kmer = new String(AbstractLimitedSequence.getSequence(dnaHash, kmerLength, DNASequence.EMPTY_DNA_SEQUENCE));
			count = kmersMap.getCount(kmer);
		}
		long hash = Integer.MAX_VALUE;
		if(count > 0) {
			long rankingStart=kmersAnalyzer.getRanking(count);
			long kmersWithCount = kmersAnalyzer.getNumKmers(count);
			if(kmersWithCount<primeNumbersHelper.getCapacity()) prime = primeNumbersHelper.getNextPrime((int) kmersWithCount);
			hash = rankingStart+(dnaHash%prime);
			if(hash>=Integer.MAX_VALUE) hash = Integer.MAX_VALUE-1;
		}
		return (int)hash;
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
		long [] segmentCodes = KmersExtractor.extractDNAKmerCodes(query.toString(), kmerLength, 0, query.length());
		Map<Integer, Long> selectedCodes = new LinkedHashMap<>();
		List<KmerCodesTableEntry> selectedCodeEntries = computeSequenceCodes(queryIdx, 0, query.length(), segmentCodes);
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
		//int idxDebug = -1;
		
		if (queryIdx == idxDebug) System.out.println("ShortKmerCodesTable. Aligning a total of "+codes.size()+" codes. Mode: "+mode+" kmer length: "+kmerLength);
		//int limitSequences = Math.max(sequenceLengths.size()/10, 4*mode);
		int limitSequences = Math.max(100, 4*mode);
		
		Set<Integer> preselectedSubjectIds = null;
		if(queryLength>100000) preselectedSubjectIds = preselectSubjectIds(queryLength, limitSequences, limitHitsPerSequence, codes);
		KmerSearchResultsCompressedTable result = new KmerSearchResultsCompressedTable(codes, kmerLength, 10);
		Set<Long> internalMultiHitCodes = calculateInternalMultihitKmers(codes);
		/*Map<Long,Integer> codesLocalCounts = new HashMap<Long, Integer>();
		for(Long code:codes.values()) {
			codesLocalCounts.compute(code, (k,v)->(v==null?1:v+1));
		}
		if (queryIdx == idxDebug) System.out.println("Minimizers table. Counting hits for query. Codes: "+codes.size()+" unique: "+codesLocalCounts.size());
		*/
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
			//int count = codesLocalCounts.getOrDefault(kmerCode, 0);
			int countSeqs = getCountDifferentSequences(kmerCode);
			//if (queryIdx == idxDebug && (countSeqs>10 || startQuery==0)) System.out.println("Minimizers table. For pos "+startQuery+" kmer: "+new String (DNASequence.getDNASequence(kmerCode, kmerLength))+" count sequences: "+countSeqs+" limit "+limitSequences);
			if (queryIdx == idxDebug ) System.out.println("Minimizers table. For pos "+startQuery+" kmer: "+new String (DNASequence.getDNASequence(kmerCode, kmerLength))+" count sequences: "+countSeqs+" limit "+limitSequences);
			if (countSeqs>limitSequences) {
				multiSequenceCodes++;
				continue;
			}
			
			long [] codesMatching = lookupHits(kmerCode);
			int numHits = codesMatching.length;
			int normalized = (int)Math.round((double)numHits/mode);
			if(normalized<normalizedCountsDist.length) normalizedCountsDist[normalized]++;
			else normalizedCountsDist[normalizedCountsDist.length-1]++;
			if (queryIdx == idxDebug && startQuery==0) System.out.println("Minimizers table. For pos "+startQuery+" kmer: "+new String (DNASequence.getDNASequence(kmerCode, kmerLength))+" codes matching: "+codesMatching.length+" limit: "+(limitHitsPerSequence*countSeqs));
			if(numHits>limitHitsPerSequence*countSeqs) {
				multihitCodes++;
				continue;
			}
			else if(numHits>0) {
				usedCodes.add(kmerCode);
				 
				numUsedCodes++;
			}
			else notFoundCodes++;
			if(normalized>=2) {
				highDepthUsedCodes.add(kmerCode);
				if(internalMultiHitCodes.contains(kmerCode)) internalMultiUsedCodes.add(kmerCode);
			}
			for(long entryCode:codesMatching) {
				int [] dec = KmerCodesTableEntry.decode(entryCode);
				int subjectIdx = dec[0];
				if (queryIdx == idxDebug && startQuery==0) System.out.println("Minimizers table. For pos "+startQuery+" kmer: "+new String (DNASequence.getDNASequence(kmerCode, kmerLength))+" next match: "+subjectIdx+" start: "+dec[1]);
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
		result.setKmerWeights(calculateCodeWeights(codes));
		result.setHighDepthUsedKmerCodes(highDepthUsedCodes);
		result.setInternalMultihitUsedKmerCodes(internalMultiUsedCodes);
		if (queryIdx == idxDebug) System.out.println("ShortKmerCodesTable. Total codes used: "+numUsedCodes+" not found: "+notFoundCodes+" self sequence count: "+selfSequenceCount+" codes with hits in multiple sequences: "+multiSequenceCodes);
		return result;
		
	}
	
	private Set<Integer> preselectSubjectIds(int queryLength, int limitSequences, int limitHitsPerSequence, Map<Integer, Long> codes) {
		int minHits = queryLength/100;
		
		Map<Integer,Integer> subjectHitCounts = new HashMap<>();
		for(Map.Entry<Integer, Long> entry:codes.entrySet()) {
			long kmerCode = entry.getValue();
			//int count = codesLocalCounts.getOrDefault(kmerCode, 0);
			int countSeqs = getCountDifferentSequences(kmerCode);
			if (countSeqs>limitSequences) continue;
			
			long [] codesMatching = lookupHits(kmerCode);
			if(codesMatching.length>limitHitsPerSequence*countSeqs) continue;
			
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
	private Set<Long> calculateInternalMultihitKmers(Map<Integer, Long> searchCodes) {
		Set<Long> answer = new HashSet<>();
		Map<Long,Integer> countsMap = new HashMap<>();
		for(Map.Entry<Integer,Long> entry:searchCodes.entrySet()) {
			countsMap.compute(entry.getValue(), (k,v)->v!=null?v+1:1);
		}
		for(Map.Entry<Long,Integer> entry:countsMap.entrySet()) {
			if(entry.getValue()>2) {
				answer.add(entry.getKey());
			}
		}
		return answer;
	}
	public Map<Long,Double> calculateCodeWeights(Map<Integer, Long> codes) {
		Map<Long,Double> answer = new HashMap<>();
		for(long code:codes.values()) {
			answer.computeIfAbsent(code, v->calculateWeight(code));
		}
		return answer;
	}

	public double calculateWeight(long code) {
		//if(kmersMap==null) return 1;
		int countDifferent = getCountDifferentSequences(code);
		/*int totalCount = getTotalHits(minimizer);
		int diff1 = countDifferent-mode;
		int diff2 = totalCount/countQuery-mode;
		if(diff1<=kmerDistModeLocalSD && diff2<=kmerDistModeLocalSD) return 1;
		int diff3=diff1+diff2-2*kmerDistModeLocalSD;
		if(diff3<1) diff3=1;*/
		int modeMinimizers = Math.max(1, mode/2);
		int diff1 = countDifferent-modeMinimizers;
		if(diff1<=kmerDistModeLocalSD) return 1;
		int diff3=diff1-kmerDistModeLocalSD;
		return 1.0*modeMinimizers/(modeMinimizers+diff3);
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
