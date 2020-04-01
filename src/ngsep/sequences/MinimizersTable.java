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
		Map<Integer, List<MinimizersTableEntry>> minimizersSeq = computeSequenceMinimizers(sequenceId, sequence);
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
		if ((sequenceId+1)%100==0) System.out.println("Added "+numSequences+" sequences. Total minimizers:"+sequencesByMinimizer.size());
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
	public Map<Integer, List<MinimizersTableEntry>> computeSequenceMinimizers(int sequenceId, CharSequence sequence) {
		Map<Integer, CharSequence> kmers = KmersExtractor.extractKmersAsMap(sequence.toString(), kmerLength, 1, 0, sequence.length(), false, true, true);
		return computeSequenceMinimizers(sequenceId, sequence.length(), kmers);
	}
	/**
	 * Calculates the minimizers of the sequence represented by the given kmers
	 * @param sequenceId Id of the sequence to calculate
	 * @param sequenceLength Length of the sequence
	 * @param sequenceKmers Kmers considered to build minimizers
	 * @return Map<Integer, List<MinimizersTableEntry>> Minimizers calculated for the given sequence indexed by the minimizer
	 */
	public Map<Integer, List<MinimizersTableEntry>> computeSequenceMinimizers(int sequenceId, int sequenceLength, Map<Integer, CharSequence> sequenceKmers) {
		Map<Integer, List<MinimizersTableEntry>> minimizersSeq = new HashMap<Integer, List<MinimizersTableEntry>>();
		Map<Integer, Integer> hashcodesForward = new HashMap<Integer, Integer>();
		for(int i: sequenceKmers.keySet()) {
			CharSequence kmer = sequenceKmers.get(i);
			if(kmer == null) continue;
			hashcodesForward.put(i, getHash(kmer));
		}
		Integer previousMinimizer = null;
		int end = sequenceLength - kmerLength - windowLength;
		for(int i=0;i<end;i++) {
			int minimizerI = 0;
			int minJ = -1;
			for(int j=0;j<windowLength;j++) {
				Integer hashForward = hashcodesForward.get(i+j);
				if (hashForward!=null && (minJ==-1 || hashForward <= minimizerI)) {
					minimizerI = hashForward;
					minJ = j;
				}
			}
			if(minJ == -1) {
				previousMinimizer = null;
				continue;
			}
			if (previousMinimizer!=null && minimizerI==previousMinimizer) continue;
			List<MinimizersTableEntry> minList = minimizersSeq.computeIfAbsent(minimizerI, l -> new ArrayList<MinimizersTableEntry>());
			MinimizersTableEntry entry = new MinimizersTableEntry(minimizerI, sequenceId, i+minJ);
			minList.add(entry);
			previousMinimizer = minimizerI;
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
		Map<Integer, List<MinimizersTableEntry>> minimizersQuery = computeSequenceMinimizers(-1, query.length(), kmers);
		//if (querySequenceId == idxDebug) log.info("Minimizers table. Counting hits for query: "+querySequenceId+" "+queryRC+" queryCount: "+queryCount+" min count: "+minCount);
		
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
				
				UngappedSearchHit hit = new UngappedSearchHit(kmer, subjectIdx, matchingEntry.getStart());
				hit.setQueryIdx(queryEntry.getStart());
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
