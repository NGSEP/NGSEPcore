package ngsep.sequences;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.math.Distribution;
import ngsep.math.PrimeNumbers;

public class MinimizersTable {
	
	private static final long [] EMPTY_LONG_ARRAY = new long[0];
	
	private Logger log = Logger.getLogger(MinimizersTable.class.getName());
	
	private int kmerLength;
	private int windowLength;
	
	
	private KmersMapAnalyzer kmersAnalyzer;
	private KmersMap kmersMap;
	private Map<Long,Integer> explicitKmerHashCodes = new HashMap<Long, Integer>();
	private int mode=1;
	private int kmerDistModeLocalSD=5;
	
	//Structures to implement the minimizers hash table
	//Map with minimizer as key and row of the sequencesByMinimizerTable as value
	private Map<Integer,Integer> matrixRowMap;
	//Table with encoded entries. Each row corresponds to a single minimizer
	private long [][] sequencesByMinimizerTable;
	//Actual lengths of the lists within the table
	private short [] sequencesByMinimizerTableColumnLengths;
	//Count of different sequences reporting each minimizer
	private short [] minimizerCountDifferentSequences;
	
	private Map<Integer,Integer> sequenceLengths = new HashMap<Integer, Integer>();
	private long totalEntries = 0;
	private PrimeNumbers primeNumbersHelper;
	
	

	public MinimizersTable(int kmerLength, int windowLength) {
		this.kmerLength = kmerLength;
		this.windowLength = windowLength;
		this.primeNumbersHelper = new PrimeNumbers();
		initializeTable(100000);
	}
	public MinimizersTable(KmersMapAnalyzer kmersAnalyzer, int kmerLength, int windowLength) {
		this.kmersAnalyzer = kmersAnalyzer;
		this.kmersMap = kmersAnalyzer.getKmersMap();
		this.kmerLength = kmerLength;
		this.windowLength = windowLength;
		this.mode = kmersAnalyzer.getMode();
		this.kmerDistModeLocalSD = kmersAnalyzer.getModeLocalSD();
		//TODO: This is ok for assemblies but it is not ok for reference genomes
		int capacityPrimes = (int) Math.min(100000000, kmersAnalyzer.getNumKmers(this.mode));
		this.primeNumbersHelper = new PrimeNumbers(capacityPrimes);
		
		/*if(!kmersAnalyzer.isAssembly()) {
			long [] codesUniqueZone = kmersAnalyzer.extractKmerCodesInLocalSDZone();
			for(int i=0;i<codesUniqueZone.length && codesUniqueZone[i]>=0;i++) {
				explicitKmerHashCodes.put(codesUniqueZone[i],i);
			}
			System.out.println("Indexed "+explicitKmerHashCodes.size()+" kmer unique codes.");
		}*/
		
		//Create the structures with appropriate initial capacity
		int capacity = kmersMap.size()/10;
		initializeTable(capacity);
	}
	public void initializeTable(int capacity) {
		matrixRowMap = new HashMap<Integer, Integer>(capacity);
		sequencesByMinimizerTable = new long [capacity][1];
		for(int i=0;i<sequencesByMinimizerTable.length;i++) Arrays.fill(sequencesByMinimizerTable[i], 0);
		sequencesByMinimizerTableColumnLengths = new short [capacity];
		Arrays.fill(sequencesByMinimizerTableColumnLengths, (short)0);
		minimizerCountDifferentSequences = new short [capacity];
		Arrays.fill(minimizerCountDifferentSequences, (short)0);
	}
	
	//Hash table management methods
	private long[] lookupHits(int minimizer) {
		Integer row = matrixRowMap.get(minimizer);
		if(row==null) return EMPTY_LONG_ARRAY;
		return Arrays.copyOf(sequencesByMinimizerTable[row], sequencesByMinimizerTableColumnLengths[row]);
	}
	public int size() {
		return matrixRowMap.size();
	}
	
	/**
	 * Calculates the number of times that the minimizer has been observed
	 * @param minimizer number to query
	 * @return int times that the given minimizer has been observed
	 */
	public int getTotalHits(int minimizer) {
		Integer row = matrixRowMap.get(minimizer);
		if(row==null) return 0;
		return sequencesByMinimizerTableColumnLengths[row];
	}
	
	/**
	 * Calculates the number of sequences where the minimizer has been observed
	 * @param minimizer number to query
	 * @return int number of different sequences where the minimizer has been observed
	 */
	public int getCountDifferentSequences(int minimizer) {
		Integer row = matrixRowMap.get(minimizer);
		if(row==null) return 0;
		return minimizerCountDifferentSequences[row];
	}
	
	private void addMinimizerSequence (int minimizer, List<KmerCodesTableEntry> entries) {
		Integer row = matrixRowMap.get(minimizer);
		if(row==null) {
			row = size();
			if(row ==Integer.MAX_VALUE) {
				log.warning("Reached maximum number of minimizers that can be saved "+row);
				return;
			}
			matrixRowMap.put(minimizer, row);
			if(row==sequencesByMinimizerTable.length) resizeTable();
		}
		int currentCount = sequencesByMinimizerTableColumnLengths[row];
		if (currentCount+entries.size()<Short.MAX_VALUE) {
			for (KmerCodesTableEntry entry:entries) addToTable(row, entry.encode());
			totalEntries+=entries.size();
			minimizerCountDifferentSequences[row]++;
		}
	}
	private void resizeTable() {
		log.info("Resizing minimizers table. Current number of minimizers: "+size()+" current capacity: "+sequencesByMinimizerTable.length);
		int newCapacity =  2*sequencesByMinimizerTable.length;
		if(newCapacity<0) newCapacity = Integer.MAX_VALUE;
		sequencesByMinimizerTableColumnLengths = Arrays.copyOf(sequencesByMinimizerTableColumnLengths, newCapacity);
		minimizerCountDifferentSequences = Arrays.copyOf(minimizerCountDifferentSequences, newCapacity);
		long [][] newTable = new long [newCapacity][0];
		for(int i=0;i<newCapacity;i++) {
			if(i<sequencesByMinimizerTable.length) newTable[i] = sequencesByMinimizerTable[i];
			else newTable[i] = new long[1];
		}
		sequencesByMinimizerTable = newTable;
		log.info("Resized minimizers table. New capacity: "+sequencesByMinimizerTable.length);
	}
	private void addToTable(int row, long value) {
		long [] minimizerEntries = sequencesByMinimizerTable[row];
		int column = sequencesByMinimizerTableColumnLengths[row];
		if(column == minimizerEntries.length ) {
			//Resize entries
			sequencesByMinimizerTable[row] = Arrays.copyOf(minimizerEntries, 2*minimizerEntries.length);
		}
		sequencesByMinimizerTable[row][column] = value;
		sequencesByMinimizerTableColumnLengths[row]++;
	}
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	
	public KmersMap getKmersMap() {
		return kmersMap;
	}
	
	public void setKmersMap(KmersMap kmersMap) {
		this.kmersMap = kmersMap;
	}
	/**
	 * Adds the minimizers of the given sequence to the table
	 * @param sequenceId Id of the sequence to add
	 * @param sequence to add
	 */
	public void addSequence (int sequenceId, CharSequence sequence) {
		int n = sequence.length();
		String sequenceStr = sequence.toString();
		int step = 500000;
		Map<Integer, List<KmerCodesTableEntry>> minimizersSeq = new HashMap<Integer, List<KmerCodesTableEntry>>();
		for (int start = 0;start < n;start+=step) {
			List<KmerCodesTableEntry> minEntriesList = computeSequenceMinimizers(sequenceId, sequenceStr, start, Math.min(n, start+step));
			for(KmerCodesTableEntry entry:minEntriesList) {
				List<KmerCodesTableEntry> minList = minimizersSeq.computeIfAbsent((int)entry.getKmerCode(), l->new ArrayList<KmerCodesTableEntry>());
				minList.add(entry);
			}
		}
		//log.info("Sequence "+sequenceId+" number of minimizers: "+minimizersSeq.size());
		
		synchronized (sequenceLengths) {
			for(Map.Entry<Integer, List<KmerCodesTableEntry>> minEntry:minimizersSeq.entrySet()) {
				if (minEntry.getValue().size()== 0) continue;
				addMinimizerSequence (minEntry.getKey(), minEntry.getValue());
			}
			sequenceLengths.put(sequenceId, n);
		}	
	}

	/**
	 * Calculates the minimizers of the given sequence
	 * @param sequenceId Id of the sequence to calculate
	 * @param sequence characters of the sequence to calculate
	 * @return Map<Integer, List<MinimizersTableEntry>> Minimizers calculated for the given sequence indexed by the minimizer
	 */
	public List<KmerCodesTableEntry> computeSequenceMinimizers(int sequenceId, String sequence,int start,int end) {
		//Map<Integer, String> kmers = KmersExtractor.extractKmersAsMap(sequence.toString(), kmerLength, 1, start, Math.min(sequence.length(),end+windowLength+kmerLength), false, true, true);
		Map<Integer, Long> codes = KmersExtractor.extractDNAKmerCodes(sequence, kmerLength, start, Math.min(sequence.length(),end+windowLength+kmerLength));
		//log.info("Extracted codes for sequence "+sequenceId+" from "+start+" to "+end+" Nuber of codes: "+codes.size());
		return computeSequenceMinimizers(sequenceId, start, Math.min(end, sequence.length()-kmerLength-windowLength), codes);
	}
	/**
	 * Calculates the minimizers of the sequence represented by the given kmers
	 * @param sequenceId Id of the sequence to calculate
	 * @param start of the sequence to consider
	 * @param end of the sequence to consider
	 * @param sequenceKmers Kmers considered to build minimizers
	 * @return Map<Integer, List<MinimizersTableEntry>> Minimizers calculated for the given sequence indexed by the minimizer
	 */
	private List<KmerCodesTableEntry> computeSequenceMinimizers(int sequenceId, int start, int end, Map<Integer, Long> kmerCodes) {
		int debugIdx = -2;
		List<KmerCodesTableEntry> minimizersSeq = new ArrayList<KmerCodesTableEntry>();
		Map<Integer, Integer> hashcodesForward = new HashMap<Integer, Integer>();
		for(Map.Entry<Integer, Long> entry: kmerCodes.entrySet()) {
			int i = entry.getKey();
			long code = entry.getValue();
			int hash = getHash(code);
			hashcodesForward.put(i, hash);
		}
		if(sequenceId==debugIdx) System.err.println("Filtered codes for sequence "+sequenceId+" from "+start+" to "+end+" Filtered codes: "+hashcodesForward.size());
		Integer previousMinimizer = null;
		int previousMinimizerPos = -1;
		for(int i=start;i<end;i++) {
			Integer minimizerI = null;
			int minPos = -1;
			Integer newHash = hashcodesForward.get(i+windowLength-1);
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
					Integer hashForward = hashcodesForward.get(i+j);
					if (hashForward!=null && (minimizerI==null || hashForward <= minimizerI)) {
						minimizerI = hashForward;
						minPos = i+j;
					}
				}
				//if(sequenceId==debugIdx && i>0 && i<3000) System.err.println("Minimizer calculated with cycle. Start: "+i+" New pos: "+minPos+" new minimizer: "+minimizerI+" previous: "+previousMinimizer+" kmer code: "+kmerCodes.get(minPos)+" total: "+minimizersSeq.size());
			}
			if (minimizerI==previousMinimizer) continue;
			if(minimizerI != null) {
				KmerCodesTableEntry entry = new KmerCodesTableEntry(minimizerI, sequenceId, minPos);
				minimizersSeq.add(entry);
			}
			previousMinimizer = minimizerI;
			previousMinimizerPos = minPos;
		}
		if(sequenceId==debugIdx) System.err.println("Calculated minimizers for sequence "+sequenceId+" from "+start+" to "+end+" Minimizers: "+minimizersSeq.size());
		return minimizersSeq;
	}

	private int getHash(long dnaHash) {
		Integer code = explicitKmerHashCodes.get(dnaHash);
		if(code!=null) return code;
		int prime = 1073676287;
		//if(kmersAnalyzer!=null) {
		if(kmersAnalyzer==null) {
			long answer = (dnaHash+1)%prime;
			return (int) answer;
		}
		int count;
		if(kmersMap instanceof ShortArrayDNAKmersMapImpl) {
			count = ((ShortArrayDNAKmersMapImpl)kmersMap).getCount(dnaHash);
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
		/*int distance = count-mode;
		long hash;
		if(distance>0) {
			hash = distance << 24;
			hash = Math.max(hash, explicitKmerHashCodes.size());
			hash+= (dnaHash & 0xFFFFFFF);
			if(hash>Integer.MAX_VALUE) hash = Integer.MAX_VALUE;
		} else {
			hash = (-distance) << 25;
			hash+= (dnaHash & 0xFFFFFFF);
			if(hash>Integer.MAX_VALUE) hash = Integer.MAX_VALUE;
		}*/
		//explicitKmerHashCodes.put(dnaHash, (int)hash);
		return (int)hash;
	}
	
	
	/**
	 * Calculates the hits of the given query
	 * @param query sequence
	 * @return Map<Integer,List<MinimizersTableEntry>> Sequences matching kmers of the given query indexed by subject and sorted by subject start position
	 */
	public Map<Integer,List<UngappedSearchHit>> match (int queryIdx, CharSequence query) {
		Map<Integer, Long> codes = KmersExtractor.extractDNAKmerCodes(query.toString(), kmerLength, 0, query.length());
		//computeSequenceMinimizers(-1, 0, query.length(), codes);
		return match(queryIdx, query.length(), codes);
		//return new HashMap<Integer, List<UngappedSearchHit>>();
		/*Random r = new Random();
		Map<Integer, Long> selectedCodes = new HashMap<Integer, Long>();
		for(Map.Entry<Integer, Long> entry:codes.entrySet()) {
			if(r.nextDouble()<0.2) selectedCodes.put(entry.getKey(), entry.getValue());
		}
		return match(queryIdx, query.length(), selectedCodes);
		*/
	}
	/**
	 * Calculates the hits of the given query
	 * @param query sequence
	 * @return Map<Integer,List<MinimizersTableEntry>> Sequences matching kmers of the given query indexed by subject and sorted by subject start position
	 */
	public Map<Integer,List<UngappedSearchHit>> match (int queryIdx, int queryLength, Map<Integer, Long> codes) {
		int idxDebug = -2;
		//int idxDebug = 1;
		//int limitSequences = Math.max(sequenceLengths.size()/10, 4*mode);
		int limitSequences = Math.max(100, 4*mode);
		List<KmerCodesTableEntry> minimizersQueryList = computeSequenceMinimizers(-1, 0, queryLength, codes);
		
		Map<Integer,Integer> minimizersLocalCounts = new HashMap<Integer, Integer>();
		for(KmerCodesTableEntry entry:minimizersQueryList) {
			minimizersLocalCounts.compute((int)entry.getKmerCode(), (k,v)->(v==null?1:v+1));
		}
		if (queryIdx == idxDebug) {
			Set<Long> uniqueCodes = new HashSet<Long>();
			uniqueCodes.addAll(codes.values());
			System.out.println("Minimizers table. Counting hits for query. Codes: "+codes.size()+" unique: "+uniqueCodes.size()+" minimizer counts. total: "+minimizersQueryList.size()+" unique: "+minimizersLocalCounts.size());
		}
		int numUsedMinimizers = 0;
		int multihitMinimizers = 0;
		int withoutkmerMinimizers = 0;
		int selfSequenceCount = 0;
		Map<Integer,List<UngappedSearchHit>> answer = new HashMap<Integer, List<UngappedSearchHit>>();
		//Set<Integer> usedMinimizers = new HashSet<Integer>();
		for(KmerCodesTableEntry entry:minimizersQueryList) {
			int minimizer = (int)entry.getKmerCode();
			//if(usedMinimizers.contains(minimizer)) continue;
			//usedMinimizers.add(minimizer);
			int count = minimizersLocalCounts.getOrDefault(minimizer, 0);
			int countSeqs = getCountDifferentSequences(minimizer);
			//if (queryIdx == idxDebug && count>1) System.out.println("Minimizers table. For minimizer: "+minimizer+" query entries: "+count+" count sequences: "+countSeqs+" mode "+mode);
			if (countSeqs>limitSequences) {
				multihitMinimizers++;
				continue;
			}
			Long kmerCode = codes.get(entry.getStart());
			if(kmerCode == null) {
				//Kmers that are not a minimizers are not considered
				withoutkmerMinimizers++;
				continue;
			}
			CharSequence kmer = new String(AbstractLimitedSequence.getSequence(kmerCode, kmerLength, DNASequence.EMPTY_DNA_SEQUENCE));
			long [] codesMatching = lookupHits(minimizer);
			if(codesMatching.length>0) numUsedMinimizers++;
			for(long entryCode:codesMatching) {
				KmerCodesTableEntry matchingEntry = new KmerCodesTableEntry(minimizer, entryCode);
				int subjectIdx = matchingEntry.getSequenceId();
				if (subjectIdx < 0) {
					System.err.println("Invalid subject "+subjectIdx+" minimizer: "+minimizer+" matching code: "+entryCode+" start: "+matchingEntry.getStart());
					continue;
				}
				UngappedSearchHit hit = new UngappedSearchHit(kmer, subjectIdx, matchingEntry.getStart());
				hit.setQueryIdx(entry.getStart());
				hit.setWeight(calculateWeight(minimizer, count));
				List<UngappedSearchHit> targetHits = answer.computeIfAbsent(subjectIdx,l -> new ArrayList<UngappedSearchHit>());
				targetHits.add(hit);
				if(subjectIdx==queryIdx) selfSequenceCount++;
			}
		}
		if (queryIdx == idxDebug) System.out.println("Minimizers table. Total minimizers used: "+numUsedMinimizers+" self sequence count: "+selfSequenceCount+" minimizers with multiple hits: "+multihitMinimizers+" kmer not found: "+withoutkmerMinimizers);
		return answer;
		
	}

	private double calculateWeight(int minimizer, int countQuery) {
		//if(kmersMap==null) return 1;
		int countDifferent = getCountDifferentSequences(minimizer);
		/*int totalCount = getTotalHits(minimizer);
		int diff1 = countDifferent-mode;
		int diff2 = totalCount/countQuery-mode;
		if(diff1<=kmerDistModeLocalSD && diff2<=kmerDistModeLocalSD) return 1;
		int diff3=diff1+diff2-2*kmerDistModeLocalSD;
		if(diff3<1) diff3=1;*/
		int modeMinimizers = mode/2;
		int diff1 = countDifferent-modeMinimizers;
		if(diff1<=kmerDistModeLocalSD) return 1;
		int diff3=diff1-kmerDistModeLocalSD;
		return 1.0*modeMinimizers/(modeMinimizers+diff3);
	}
	public Distribution calculateDistributionHits() {
		Distribution dist = new Distribution(1, 300, 1);
		int numMinimizers = size();
		for(int i=0;i<numMinimizers;i++) {
			dist.processDatapoint(sequencesByMinimizerTableColumnLengths[i]);	
		}
		return dist;
	}

	public long getTotalEntries() {
		return totalEntries;
	}
	
	
}
