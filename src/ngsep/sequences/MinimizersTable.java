package ngsep.sequences;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.math.Distribution;

public class MinimizersTable {
	
	private static final long [] EMPTY_LONG_ARRAY = new long[0];
	
	private Logger log = Logger.getLogger(MinimizersTable.class.getName());
	
	private int kmerLength;
	private int windowLength;
	private int maxAbundanceMinimizer = 0;
	private boolean saveRepeatedMinimizersWithinSequence = false;
	private boolean keepSingletons = false;
	
	
	
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
	
	private Map<Integer,Integer> sequenceLengths = new HashMap<Integer, Integer>();
	private long totalEntries = 0;
	
	

	public MinimizersTable(int kmerLength, int windowLength) {
		this.kmerLength = kmerLength;
		this.windowLength = windowLength;
		initializeTable(100000);
	}
	public MinimizersTable(KmersMapAnalyzer kmersAnalyzer, int kmerLength, int windowLength) {
		this.kmersMap = kmersAnalyzer.getKmersMap();
		this.kmerLength = kmerLength;
		this.windowLength = windowLength;
		this.mode = kmersAnalyzer.getMode();
		this.kmerDistModeLocalSD = kmersAnalyzer.getModeLocalSD();
		long [] codesUniqueZone = kmersAnalyzer.extractKmerCodesInLocalSDZone();
		for(int i=0;i<codesUniqueZone.length && codesUniqueZone[i]>=0;i++) {
			explicitKmerHashCodes.put(codesUniqueZone[i],i);
		}
		System.out.println("Indexed "+explicitKmerHashCodes.size()+" kmer unique codes.");
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
	
	private synchronized void add (int minimizer, List<MinimizersTableEntry> entries) {
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
		if (currentCount+entries.size()<Short.MAX_VALUE && (maxAbundanceMinimizer==0 || currentCount < maxAbundanceMinimizer)) {
			if(saveRepeatedMinimizersWithinSequence) {
				for (MinimizersTableEntry entry:entries) addToTable(row, entry.encode());
				totalEntries+=entries.size();
			} else {
				addToTable(row, entries.get(0).encode());
				totalEntries++;
			}
		}
	}
	private void resizeTable() {
		log.info("Resizing minimizers table. Current number of minimizers: "+size()+" current capacity: "+sequencesByMinimizerTable.length);
		int newCapacity =  2*sequencesByMinimizerTable.length;
		if(newCapacity<0) newCapacity = Integer.MAX_VALUE;
		sequencesByMinimizerTableColumnLengths = Arrays.copyOf(sequencesByMinimizerTableColumnLengths, newCapacity);
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
	public int getMaxAbundanceMinimizer() {
		return maxAbundanceMinimizer;
	}
	public void setMaxAbundanceMinimizer(int maxAbundanceMinimizer) {
		this.maxAbundanceMinimizer = maxAbundanceMinimizer;
	}
	
	public boolean isSaveRepeatedMinimizersWithinSequence() {
		return saveRepeatedMinimizersWithinSequence;
	}
	public void setSaveRepeatedMinimizersWithinSequence(boolean saveRepeatedMinimizersWithinSequence) {
		this.saveRepeatedMinimizersWithinSequence = saveRepeatedMinimizersWithinSequence;
	}
	
	public KmersMap getKmersMap() {
		return kmersMap;
	}
	
	public boolean isKeepSingletons() {
		return keepSingletons;
	}

	public void setKeepSingletons(boolean keepSingletons) {
		this.keepSingletons = keepSingletons;
	}

	public void addSequences (QualifiedSequenceList sequences) {
		for(int i=0;i<sequences.size();i++) {
			addSequence(i, sequences.get(i).getCharacters());
		}
	}
	/**
	 * Adds the minimizers of the given sequence to the table
	 * @param sequenceId Id of the sequence to add
	 * @param sequence to add
	 */
	public void addSequence (int sequenceId, CharSequence sequence) {
		int n = sequence.length();
		String sequenceStr = sequence.toString();
		int step = 10000000;
		Map<Integer, List<MinimizersTableEntry>> minimizersSeq = new HashMap<Integer, List<MinimizersTableEntry>>();
		for (int start = 0;start < n;start+=step) {
			List<MinimizersTableEntry> minEntriesList = computeSequenceMinimizers(sequenceId, sequenceStr, start, Math.min(n, start+step));
			for(MinimizersTableEntry entry:minEntriesList) {
				List<MinimizersTableEntry> minList = minimizersSeq.computeIfAbsent(entry.getMinimizer(), l->new ArrayList<MinimizersTableEntry>());
				minList.add(entry);
			}
		}
		
		for(int minimizer:minimizersSeq.keySet()) {
			List<MinimizersTableEntry> entries = minimizersSeq.get(minimizer);
			if (entries.size()== 0) continue; 
			if(!saveRepeatedMinimizersWithinSequence && !overlapping(entries)) continue;
			add (minimizer, entries);
		}
		synchronized (sequenceLengths) {
			sequenceLengths.put(sequenceId, n);
		}	
	}

	private boolean overlapping(List<MinimizersTableEntry> entries) {
		if(entries.size()<2) return true;
		int firstStart = entries.get(0).getStart();
		for(MinimizersTableEntry entry:entries) {
			if(Math.abs(firstStart-entry.getStart())>2*windowLength) return false;
		}
		return true;
	}

	/**
	 * Calculates the minimizers of the given sequence
	 * @param sequenceId Id of the sequence to calculate
	 * @param sequence characters of the sequence to calculate
	 * @return Map<Integer, List<MinimizersTableEntry>> Minimizers calculated for the given sequence indexed by the minimizer
	 */
	public List<MinimizersTableEntry> computeSequenceMinimizers(int sequenceId, String sequence,int start,int end) {
		//Map<Integer, String> kmers = KmersExtractor.extractKmersAsMap(sequence.toString(), kmerLength, 1, start, Math.min(sequence.length(),end+windowLength+kmerLength), false, true, true);
		Map<Integer, Long> codes = KmersExtractor.extractDNAKmerCodes(sequence, kmerLength, start, Math.min(sequence.length(),end+windowLength+kmerLength));
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
	private List<MinimizersTableEntry> computeSequenceMinimizers(int sequenceId, int start, int end, Map<Integer, Long> kmerCodes) {
		List<MinimizersTableEntry> minimizersSeq = new ArrayList<MinimizersTableEntry>();
		Map<Integer, Integer> hashcodesForward = new HashMap<Integer, Integer>();
		for(int i: kmerCodes.keySet()) {
			long code = kmerCodes.get(i);
			if(kmersMap!=null && kmersMap instanceof ShortArrayDNAKmersMapImpl) {
				int count = ((ShortArrayDNAKmersMapImpl)kmersMap).getCount(code);
				if(!keepSingletons && count == 1) continue;
				if(maxAbundanceMinimizer >0 && count>maxAbundanceMinimizer) continue;
			}
			hashcodesForward.put(i, getHash(code));
		}
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
			}
			if (minimizerI==previousMinimizer) continue;
			if(minimizerI != null) {
				MinimizersTableEntry entry = new MinimizersTableEntry(minimizerI, sequenceId, minPos);
				minimizersSeq.add(entry);
			}
			previousMinimizer = minimizerI;
			previousMinimizerPos = minPos;
		}
		return minimizersSeq;
	}

	private int getHash(long dnaHash) {
		Integer code = explicitKmerHashCodes.get(dnaHash);
		if(code!=null) return code;
		if(kmersMap==null) {
			if(dnaHash<=Integer.MAX_VALUE) return (int)dnaHash;
			return (int) dnaHash %100000000;
		}
		int count;
		if(kmersMap instanceof ShortArrayDNAKmersMapImpl) {
			count = ((ShortArrayDNAKmersMapImpl)kmersMap).getCount(dnaHash);
		} else {
			String kmer = new String(AbstractLimitedSequence.getSequence(dnaHash, kmerLength, DNASequence.EMPTY_DNA_SEQUENCE));
			count = kmersMap.getCount(kmer);
		}
		int distance = count-mode;
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
		}
		
		return (int)hash;
	}
	
	
	/**
	 * Calculates the hits of the given query
	 * @param query sequence
	 * @return Map<Integer,List<MinimizersTableEntry>> Sequences matching kmers of the given query indexed by subject and sorted by subject start position
	 */
	public Map<Integer,List<UngappedSearchHit>> match (CharSequence query) {
		//Map<Integer, String> kmers = KmersExtractor.extractKmersAsMap(query.toString(), kmerLength, 1, 0, query.length(), false, true, true);
		Map<Integer, Long> codes = KmersExtractor.extractDNAKmerCodes(query.toString(), kmerLength, 0, query.length());
		List<MinimizersTableEntry> minimizersQueryList = computeSequenceMinimizers(-1, 0, query.length(), codes);
		//if (querySequenceId == idxDebug) log.info("Minimizers table. Counting hits for query: "+querySequenceId+" "+queryRC+" queryCount: "+queryCount+" min count: "+minCount);
		Map<Integer,List<MinimizersTableEntry>> minimizersQuery = new HashMap<Integer, List<MinimizersTableEntry>>();
		for(MinimizersTableEntry entry:minimizersQueryList) {
			List<MinimizersTableEntry> minList = minimizersQuery.computeIfAbsent(entry.getMinimizer(), l->new ArrayList<MinimizersTableEntry>());
			minList.add(entry);
		}
		Map<Integer,List<UngappedSearchHit>> answer = new HashMap<Integer, List<UngappedSearchHit>>();
		for(int minimizer:minimizersQuery.keySet()) {
			List<MinimizersTableEntry> queryHits = minimizersQuery.get(minimizer);
			if (queryHits == null || queryHits.size() != 1) continue;
			MinimizersTableEntry queryEntry = queryHits.get(0);
			Long kmerCode = codes.get(queryEntry.getStart());
			if(kmerCode == null) {
				//Kmers that are not a minimizers are not considered
				continue;
			}
			CharSequence kmer = new String(AbstractLimitedSequence.getSequence(kmerCode, kmerLength, DNASequence.EMPTY_DNA_SEQUENCE));
			long [] codesMatching = lookupHits(minimizer);
			for(long entryCode:codesMatching) {
				MinimizersTableEntry matchingEntry = new MinimizersTableEntry(minimizer, entryCode);
				int subjectIdx = matchingEntry.getSequenceId();
				if (subjectIdx <0) {
					System.err.println("Invalid subject "+subjectIdx+" minimizer: "+minimizer+" matching code: "+entryCode+" start: "+matchingEntry.getStart());
					continue;
				}
				UngappedSearchHit hit = new UngappedSearchHit(kmer, subjectIdx, matchingEntry.getStart());
				hit.setQueryIdx(queryEntry.getStart());
				hit.setWeight(calculateWeight(kmer));
				List<UngappedSearchHit> targetHits = answer.computeIfAbsent(subjectIdx,l -> new ArrayList<UngappedSearchHit>());
				targetHits.add(hit);
			}
		}
		return answer;
		
	}

	private double calculateWeight(CharSequence kmer) {
		if(kmersMap==null) return 1;
		int count = kmersMap.getCount(kmer);
		int diff = Math.abs(mode-count);
		if(diff<kmerDistModeLocalSD) return 1;
		int diff2=diff-kmerDistModeLocalSD;
		return 1.0*mode/(mode+2*diff2);
	}
	public Distribution calculateDistributionHits() {
		Distribution dist = new Distribution(1, Math.max(100, maxAbundanceMinimizer), 1);
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
