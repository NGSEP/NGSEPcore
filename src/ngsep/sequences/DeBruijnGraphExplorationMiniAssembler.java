package ngsep.sequences;

import java.util.Comparator;
import java.util.PriorityQueue;

public class DeBruijnGraphExplorationMiniAssembler {
	private KmersMap kmersMap;
	private int minKmerCount = 1;
	
	
	
	public DeBruijnGraphExplorationMiniAssembler(KmersMap kmersMap, int minKmerCount) {
		super();
		this.kmersMap = kmersMap;
		this.minKmerCount = minKmerCount;
	}
	public String assemble(String sourceKmer, String destKmer, int minAssemblyLength, int expectedAssemblyLength, int maxAssemblyLength) {
		int kmerLength = sourceKmer.length();
		if(expectedAssemblyLength<kmerLength) return null;
		
		//Stack<String> agenda = new Stack<>();
		PriorityQueue<String> agenda = new PriorityQueue<String>(new Comparator<String>() {
			@Override
			public int compare(String state1, String state2) {
				int score1 = getScore (state1, destKmer);
				int score2 = getScore (state2, destKmer);
				return score2-score1;
			}	
		});
		
		agenda.add(sourceKmer);
		for (int candidates = 0;agenda.size()>0 && candidates < 5000;candidates++) {
			String nextState = agenda.remove();
			//if(sourceKmerIdx==2528) System.out.println("Next state: "+nextState+" length: "+nextState.length()+" score: "+getScore(nextState, destKmer)+" agenda size: "+agenda.size());
			//Satisfability
			if((destKmer==null && nextState.length()==expectedAssemblyLength)) {
				return nextState;
			}
			if(destKmer!=null && nextState.length()>=minAssemblyLength && nextState.endsWith(destKmer)) {
				return nextState;
			}
			//Viability
			if(nextState.length()>=maxAssemblyLength) continue;
			//Next states
			String lastKmer = nextState.substring(nextState.length()-kmerLength);
			String kMinus1Mer = lastKmer.substring(1);
			String dna = DNASequence.BASES_STRING;
			for(int i=0;i<dna.length();i++) {
				char bp = dna.charAt(i);
				String nextKmer = kMinus1Mer+bp;
				if (kmersMap.getCount(nextKmer)>=minKmerCount) {
					agenda.add(nextState+bp); 
				}
			}
		}
		return null;
	}
	private static int getScore(String state, String destKmer) {
		if(destKmer==null) return 0;
		for(int i=destKmer.length();i>0;i--) {
			if(state.endsWith(destKmer.substring(0,i))) return i;
		}
		return 0;
	}
}
