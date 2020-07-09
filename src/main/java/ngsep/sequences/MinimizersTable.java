package ngsep.sequences;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.math.Distribution;

public class MinimizersTable {
	
	private Logger log = Logger.getLogger(MinimizersTable.class.getName());
	private int kmerLength;
	private int windowLength;
	private int maxAbundanceMinimizer = 0;
	private boolean saveRepeatedMinimizersWithinSequence = false;
	
	
	private Map<Integer, List<Long>> sequencesByMinimizer = new HashMap<Integer, List<Long>>();
	private Map<Integer,Integer> sequenceLengths = new HashMap<Integer, Integer>(); 
	private int numSequences = 0;
	
	

	public MinimizersTable(int kmerLength, int windowLength) {
		this.kmerLength = kmerLength;
		this.windowLength = windowLength;
		 
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
		sequenceLengths.put(sequenceId, n);
		int step = 50000000;
		Map<Integer, List<MinimizersTableEntry>> minimizersSeq = new HashMap<Integer, List<MinimizersTableEntry>>();
		for (int start = 0;start < n;start+=step) {
			List<MinimizersTableEntry> minEntriesList = computeSequenceMinimizers(sequenceId, sequence, start, Math.min(n, start+step));
			for(MinimizersTableEntry entry:minEntriesList) {
				List<MinimizersTableEntry> minList = minimizersSeq.computeIfAbsent(entry.getMinimizer(), l->new ArrayList<MinimizersTableEntry>());
				minList.add(entry);
			}
		}
		
		for(int minimizer:minimizersSeq.keySet()) {
			List<MinimizersTableEntry> entries = minimizersSeq.get(minimizer);
			if (entries.size()== 0) continue; 
			if(!saveRepeatedMinimizersWithinSequence && !overlapping(entries)) continue;
			synchronized (sequencesByMinimizer) {
				List<Long> minList = sequencesByMinimizer.computeIfAbsent(minimizer, l -> new ArrayList<Long>());
				if (maxAbundanceMinimizer==0 || minList.size()<=maxAbundanceMinimizer) {
					if(saveRepeatedMinimizersWithinSequence) {
						for (MinimizersTableEntry entry:entries) minList.add(entry.encode());
					} else {
						 minList.add(entries.get(0).encode());
					}
				}
			}
		}
		numSequences++;
		if (numSequences%100==0) System.out.println("Added "+numSequences+" sequences. Total minimizers:"+sequencesByMinimizer.size());
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
	public List<MinimizersTableEntry> computeSequenceMinimizers(int sequenceId, CharSequence sequence,int start,int end) {
		Map<Integer, CharSequence> kmers = KmersExtractor.extractKmersAsMap(sequence.toString(), kmerLength, 1, start, Math.min(sequence.length(),end+windowLength+kmerLength), false, true, true);
		return computeSequenceMinimizers(sequenceId, start, Math.min(end, sequence.length()-kmerLength-windowLength), kmers);
	}
	/**
	 * Calculates the minimizers of the sequence represented by the given kmers
	 * @param sequenceId Id of the sequence to calculate
	 * @param start of the sequence to consider
	 * @param end of the sequence to consider
	 * @param sequenceKmers Kmers considered to build minimizers
	 * @return Map<Integer, List<MinimizersTableEntry>> Minimizers calculated for the given sequence indexed by the minimizer
	 */
	public List<MinimizersTableEntry> computeSequenceMinimizers(int sequenceId, int start, int end, Map<Integer, CharSequence> sequenceKmers) {
		List<MinimizersTableEntry> minimizersSeq = new ArrayList<MinimizersTableEntry>();
		Map<Integer, Integer> hashcodesForward = new HashMap<Integer, Integer>();
		for(int i: sequenceKmers.keySet()) {
			CharSequence kmer = sequenceKmers.get(i);
			if(kmer == null) continue;
			hashcodesForward.put(i, getHash(kmer));
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

	private int getHash(CharSequence kmer) {
		//return kmer.toString().hashCode();
		return kmer.hashCode();
	}
	
	/**
	 * Calculates the number of times that the minimizer has been observed
	 * @param minimizer number to query
	 * @return int times that the given minimizer has been observed
	 */
	public int getTotalHits(int minimizer) {
		List<Long> hits = sequencesByMinimizer.get(minimizer);
		return (hits!=null)?hits.size():0;
	}
	
	/**
	 * Calculates the hits of the given query
	 * @param query sequence
	 * @return Map<Integer,List<MinimizersTableEntry>> Sequences matching kmers of the given query indexed by subject and sorted by subject start position
	 */
	public Map<Integer,List<UngappedSearchHit>> match (CharSequence query) {
		Map<Integer, CharSequence> kmers = KmersExtractor.extractKmersAsMap(query, kmerLength, 1, 0, query.length(), false, true, true);
		List<MinimizersTableEntry> minimizersQueryList = computeSequenceMinimizers(-1, 0, query.length(), kmers);
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
			CharSequence kmer = kmers.get(queryEntry.getStart());
			if(kmer == null) {
				//Neighbor kmers normally share minimizers
				continue;
			}
			List<Long> codesMatching = sequencesByMinimizer.get(minimizer);
			if (codesMatching==null) continue;
			for(long entryCode:codesMatching) {
				MinimizersTableEntry matchingEntry = new MinimizersTableEntry(minimizer, entryCode);
				int subjectIdx = matchingEntry.getSequenceId();
				if (subjectIdx <0) {
					System.err.println("Invalid subject "+subjectIdx+" minimizer: "+minimizer+" matching code: "+entryCode+" start: "+matchingEntry.getStart());
					continue;
				}
				UngappedSearchHit hit = new UngappedSearchHit(kmer, subjectIdx, matchingEntry.getStart());
				hit.setQueryIdx(queryEntry.getStart());
				hit.setSequenceLength(sequenceLengths.get(subjectIdx));
				hit.setTotalHitsQuery(codesMatching.size());
				List<UngappedSearchHit> targetHits = answer.computeIfAbsent(subjectIdx,l -> new ArrayList<UngappedSearchHit>());
				targetHits.add(hit);
			}
		}
		return answer;
		
	}

	public Distribution calculateDistributionHits() {
		Distribution dist = new Distribution(1, 100, 1);
		for(List<Long> list:sequencesByMinimizer.values()) {
			dist.processDatapoint(list.size());	
		}
		return dist;
	}

	/**
	 * Removes minimizers observed only one time 
	 */
	public void clearOverrepresentedMinimizers() {
		List<Integer> minimizers = new ArrayList<Integer>();
		minimizers.addAll(sequencesByMinimizer.keySet());
		for(int minimizer:minimizers) {
			int totalHits = getTotalHits(minimizer);
			if (maxAbundanceMinimizer>0 && totalHits>maxAbundanceMinimizer) sequencesByMinimizer.remove(minimizer); 
		}
	}

	public int getTotalMinimizers() {
		return sequencesByMinimizer.size();
	}
}
