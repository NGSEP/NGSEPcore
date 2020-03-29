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

	public void addSequence (int sequenceId, CharSequence sequence) {
		Map<Integer, CharSequence> kmers = KmersExtractor.extractKmersAsMap(sequence.toString(), kmerLength, 1, 0, sequence.length(), false, true, true);
		Map<Integer, List<MinimizersTableEntry>> minimizersSeq = computeSequenceMinimizers(sequenceId, sequence.length(), kmers);
		for(int minimizer:minimizersSeq.keySet()) {
			List<MinimizersTableEntry> entries = minimizersSeq.get(minimizer);
			if (entries.size()== 0 || !overlapping(entries)) continue; 
			synchronized (sequencesByMinimizer) {
				List<Long> minList = sequencesByMinimizer.computeIfAbsent(minimizer, l -> new ArrayList<Long>());
				if (maxAbundanceMinimizer==0 || minList.size()<=maxAbundanceMinimizer) minList.add(entries.get(0).encode());
				/*for(MinimizersTableEntry entry:entries ) {
					minList.add(entry.encode());
				}*/
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

	public Map<Integer, List<MinimizersTableEntry>> computeSequenceMinimizers(int sequenceId, int sequenceLength, Map<Integer, CharSequence> sequenceKmers) {
		Map<Integer, List<MinimizersTableEntry>> minimizersSeq = new HashMap<Integer, List<MinimizersTableEntry>>();
		Map<Integer, Integer> hashcodesForward = new HashMap<Integer, Integer>();
		//Map<Integer, Integer> hashcodesReverse = new HashMap<Integer, Integer>();
		for(int i: sequenceKmers.keySet()) {
			CharSequence kmer = sequenceKmers.get(i);
			if(kmer == null) continue;
			hashcodesForward.put(i, getHash(kmer));
			//CharSequence kmerR = DNAMaskedSequence.getReverseComplement(kmer);
			//hashcodesReverse.put(i, getHash(kmerR));
		}
		Integer previousMinimizer = null;
		int end = sequenceLength - kmerLength - windowLength;
		for(int i=0;i<end;i++) {
			int minimizerI = 0;
			int minJ = -1;
			for(int j=0;j<windowLength;j++) {
				Integer hashForward = hashcodesForward.get(i+j);
				//Integer hashReverse = hashcodesReverse.get(i+j);
				//if(hashcodesForward==hashcodesReverse) continue;
				if (hashForward!=null && (minJ==-1 || hashForward <= minimizerI)) {
					minimizerI = hashForward;
					minJ = j;
				}
				/*if (hashReverse!=null && (minJ==-1 || hashReverse < minimizerI)) {
					minimizerI = hashReverse;
					minJ = j;
				}*/
			}
			if(minJ == -1) {
				previousMinimizer = null;
				continue;
			}
			if (previousMinimizer!=null && minimizerI==previousMinimizer) continue;
			List<MinimizersTableEntry> minList = minimizersSeq.computeIfAbsent(minimizerI, l -> new ArrayList<MinimizersTableEntry>());
			MinimizersTableEntry entry = new MinimizersTableEntry(minimizerI, sequenceId, i+minJ);
			minList.add(entry);
			
			/*for(int j=0;j<windowLength;j++) {
				Integer hashForward = hashcodesForward.get(i+j);
				if(hashForward!=null && hashForward==minimizerI) {
					
					break;
				}
				Integer hashReverse = hashcodesReverse.get(i+j);
				if(hashReverse!=null && hashReverse==minimizerI) {
					MinimizersTableEntry entry = new MinimizersTableEntry(nextSequenceId, i, true);
					minList.add(entry);
				}
			}*/
			previousMinimizer = minimizerI;
		}
		return minimizersSeq;
	}

	private int getHash(CharSequence kmer) {
		//return kmer.toString().hashCode();
		return kmer.hashCode();
	}
	
	public int getTotalHits(int minimizer) {
		List<Long> hits = sequencesByMinimizer.get(minimizer);
		return (hits!=null)?hits.size():0;
	}
	
	public Map<Integer,List<MinimizersTableEntry>> calculateMinimizerHits(int queryIdx, Map<Integer, List<MinimizersTableEntry>> minimizersQuery) {
		//System.out.println("Computed: "+minimizersSeq.size()+" minimizers for sequence "+queryIdx);
		Map<Integer,List<MinimizersTableEntry>> answer = new HashMap<Integer, List<MinimizersTableEntry>>();
		for(int minimizer:minimizersQuery.keySet()) {
			List<Long> codesMatching = sequencesByMinimizer.get(minimizer);
			if (codesMatching==null) {
				//System.out.println("Warning. Zero matches for minimizer: "+minimizer+" from sequence "+queryIdx+" matches query: "+minimizersSeq.get(minimizer).size());
				continue;
			}
			
			for(long entryCode:codesMatching) {
				MinimizersTableEntry matchingEntry = new MinimizersTableEntry(minimizer, entryCode);
				int targetId = matchingEntry.getSequenceId();
				if (matchingEntry.getSequenceId()==queryIdx) continue;
				List<MinimizersTableEntry> targetEntries = answer.computeIfAbsent(targetId,l -> new ArrayList<MinimizersTableEntry>());
				targetEntries.add(matchingEntry); 
			}
		}
		return answer;
	}

	public Distribution calculateDistributionHits() {
		Distribution dist = new Distribution(1, 50, 1);
		for(List<Long> list:sequencesByMinimizer.values()) {
			dist.processDatapoint(list.size());	
		}
		return dist;
	}

	/**
	 * Removes minimizers observed only one time 
	 */
	public void clearSingletonAndOverrepresentedMinimizers() {
		List<Integer> minimizers = new ArrayList<Integer>();
		minimizers.addAll(sequencesByMinimizer.keySet());
		for(int minimizer:minimizers) {
			int totalHits = getTotalHits(minimizer); 
			if(totalHits==1 ) sequencesByMinimizer.remove(minimizer);
			else if (maxAbundanceMinimizer>0 && totalHits>maxAbundanceMinimizer) sequencesByMinimizer.remove(minimizer); 
		}
	}

	public int getTotalMinimizers() {
		return sequencesByMinimizer.size();
	}
}
