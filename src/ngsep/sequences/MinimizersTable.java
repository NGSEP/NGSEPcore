package ngsep.sequences;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.math.Distribution;

public class MinimizersTable {
	private int kmerLength;
	private int windowLength;
	
	
	private Map<Integer, List<Integer>> sequencesByMinimizer = new HashMap<Integer, List<Integer>>();
	private List<Integer> minimizerCountsBySequence = new ArrayList<Integer>();
	
	
	

	public MinimizersTable(int kmerLength, int windowLength) {
		this.kmerLength = kmerLength;
		this.windowLength = windowLength;
		 
	}

	public void addSequence (CharSequence sequence) {
		Map<Integer, CharSequence> kmers = KmersExtractor.extractKmersAsMap(sequence, kmerLength, 1, 0, sequence.length(), false, true, true);
		Map<Integer, List<MinimizerEntry>> minimizersSeq = computeSequenceMinimizers(kmers, sequence.length());
		for(int minimizer:minimizersSeq.keySet()) {
			//List<MinimizerEntry> entries = minimizersSeq.get(minimizer);
			synchronized (sequencesByMinimizer) {
				List<Integer> minList = sequencesByMinimizer.computeIfAbsent(minimizer, l -> new ArrayList<Integer>());
				minList.add(minimizerCountsBySequence.size());
			}
		}
		minimizerCountsBySequence.add(minimizersSeq.size());
		int seqId = minimizerCountsBySequence.size();
		if (seqId%100==0) System.out.println("Added "+seqId+" sequences. Total minimizers:"+sequencesByMinimizer.size());
	}

	private Map<Integer, List<MinimizerEntry>> computeSequenceMinimizers(Map<Integer, CharSequence> kmers, int n) {
		Map<Integer, List<MinimizerEntry>> minimizersSeq = new HashMap<Integer, List<MinimizerEntry>>();
		Map<Integer, Integer> hashcodesForward = new HashMap<Integer, Integer>();
		Map<Integer, Integer> hashcodesReverse = new HashMap<Integer, Integer>();
		for(int i: kmers.keySet()) {
			CharSequence kmer = kmers.get(i);
			if(kmer == null) continue;
			hashcodesForward.put(i, getHash(kmer));
			CharSequence kmerR = DNAMaskedSequence.getReverseComplement(kmer);
			hashcodesReverse.put(i, getHash(kmerR));
		}
		int nextSequenceId = minimizerCountsBySequence.size();
		int end = n - kmerLength - windowLength;
		for(int i=0;i<end;i++) {
			int minimizerI = 0;
			int minJ = -1;
			for(int j=0;j<windowLength;j++) {
				Integer hashForward = hashcodesForward.get(i+j);
				Integer hashReverse = hashcodesReverse.get(i+j);
				if(hashcodesForward==hashcodesReverse) continue;
				if (hashForward!=null && (minJ==-1 || hashForward < minimizerI)) {
					minimizerI = hashForward;
					minJ = j;
				}
				if (hashReverse!=null && (minJ==-1 || hashReverse < minimizerI)) {
					minimizerI = hashReverse;
					minJ = j;
				}
			}
			if(minJ == -1) continue;
			List<MinimizerEntry> minList = minimizersSeq.computeIfAbsent(minimizerI, l -> new ArrayList<MinimizerEntry>());
			for(int j=0;j<windowLength;j++) {
				Integer hashForward = hashcodesForward.get(i+j);
				if(hashForward!=null && hashForward==minimizerI) {
					MinimizerEntry entry = new MinimizerEntry(nextSequenceId, i, false);
					minList.add(entry);
				}
				Integer hashReverse = hashcodesReverse.get(i+j);
				if(hashReverse!=null && hashReverse==minimizerI) {
					MinimizerEntry entry = new MinimizerEntry(nextSequenceId, i, true);
					minList.add(entry);
				}
				
			}
		}
		return minimizersSeq;
	}

	private int getHash(CharSequence kmer) {
		return kmer.toString().hashCode();
	}
	
	public int getMinimizerCountBySequenceId(int sequenceId) {
		int answer = 0;
		Integer count = minimizerCountsBySequence.get(sequenceId);
		if(count!=null) answer = count;
		return answer;
	}

	public Map<Integer,Integer> countMinimizerHits(int queryIdx, CharSequence query) {
		Map<Integer, CharSequence> kmers = KmersExtractor.extractKmersAsMap(query, kmerLength, 1, 0, query.length(), false, true, true);
		Map<Integer, List<MinimizerEntry>> minimizersSeq = computeSequenceMinimizers(kmers, query.length());
		//System.out.println("Computed: "+minimizersSeq.size()+" minimizers for sequence "+queryIdx);
		Map<Integer,Integer> answer = new HashMap<Integer, Integer>();
		for(int minimizer:minimizersSeq.keySet()) {
			List<Integer> entriesMatching = sequencesByMinimizer.get(minimizer);
			if (entriesMatching==null) {
				System.out.println("Warning. Zero matches for minimizer: "+minimizer+" from sequence "+queryIdx+" matches query: "+minimizersSeq.get(minimizer).size());
				continue;
			}
			Set<Integer> targetSeqs = new HashSet<Integer>();
			for(int sequenceId:entriesMatching) {
				//MinimizerEntry matchingEntry = new MinimizerEntry(matchingEntryCode);
				if (sequenceId==queryIdx) continue;
				targetSeqs.add(sequenceId); 
			}
			for(int targetId:targetSeqs) {
				int count = answer.computeIfAbsent(targetId, v->0);
				answer.put(targetId, count+1);
			}
		}
		return answer;
	}

	public Distribution calculateDistributionHits() {
		Distribution dist = new Distribution(1, 50, 1);
		for(List<Integer> list:sequencesByMinimizer.values()) {
			dist.processDatapoint(list.size());
			
		}
		return dist;
	}
}
class MinimizerEntry {
	private int sequenceId;
	private int start;
	private boolean reverse;
	
	
	public MinimizerEntry(int sequenceId, int start, boolean reverse) {
		this.sequenceId = sequenceId;
		this.start = start;
		this.reverse = reverse;
	}

	public MinimizerEntry (long entryCode) {
		reverse = entryCode<0;
		if(reverse) entryCode = -entryCode -1;
		start = (int) (entryCode & Integer.MAX_VALUE);
		sequenceId = (int) (entryCode >> 16);
	}

	public long encode() {
		long code = sequenceId;
		code = code << 16;
		code+=start;
		if(reverse) code = -code -1;
		return code;
	}

	public int getSequenceId() {
		return sequenceId;
	}
	public int getStart() {
		return start;
	}
	public boolean isReverse() {
		return reverse;
	}
	
}

