package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.DNASequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.KmersCounter;
import ngsep.sequences.LimitedSequence;

public class OverlapGraph {
	private List<LimitedSequence> sequences = new ArrayList<>();
	private List<List<ReadOverlap>> overlaps;
	
	public void addSequence (LimitedSequence sequence) {
		sequences.add(sequence);
	}
	
	
	
	
	
	public static void main(String[] args) {
		OverlapGraph instance = new OverlapGraph();
		int numSequences = 40;
		int totalLength = 0;
		long time1 = System.currentTimeMillis();
		Random r = new Random();
		for(int i=0;i<numSequences;i++) {
			StringBuilder randomSequence = new StringBuilder();
			int length = r.nextInt(40000)+10000;
			totalLength+=length;
			for(int j=0;j<length;j++) {
				int bpI = r.nextInt(4);
				randomSequence.append(DNASequence.BASES_STRING.charAt(bpI));
			}
			DNASequence dna = new DNASequence(randomSequence);
			instance.addSequence(dna);
			if((i+1)%1000 == 0) System.out.println("Created "+instance.sequences.size()+" random DNA sequences");
		}
		long time2 = System.currentTimeMillis();
		System.out.println("Created "+instance.sequences.size()+" random DNA sequences. Total length: "+totalLength +" time: "+ (time2 - time1));
		//System.out.println("First sequence: "+instance.sequences.get(0));
		/*time1 = time2;
		Collections.sort(instance.sequences);
		time2 = System.currentTimeMillis();
		System.out.println("Sorting time: "+ (time2 - time1));*/
		time1=time2;
		FMIndex index = new FMIndex();
		index.loadUnnamedSequences("", instance.sequences);
		time2 = System.currentTimeMillis();
		System.out.println("Created FMIndex in time: "+ (time2 - time1));
		time1=time2;
		CharSequence [] kmers = KmersCounter.extractKmers(instance.sequences.get(0), 15, true);
		time2 = System.currentTimeMillis();
		System.out.println("Extracted "+kmers.length+" k-mers from the first sequence in time: "+ (time2 - time1));
		time1=time2;
		for(int i=0;i<kmers.length;i++) {
			String kmer = kmers[i].toString();
			List<ReadAlignment> regions = index.search(kmer);
			//System.out.println("Hits: "+regions.size());
			/*for(GenomicRegion region:regions) {
				int k = Integer.parseInt(region.getSequenceName().substring(1));
				String seqOut = instance.sequences.get(k).subSequence(region.getFirst(), region.getLast()+1).toString();
				if(!seqOut.equals(kmer)) throw new RuntimeException ("Hit error in sequence: "+k+" pos: "+region.getFirst()+" search "+kmer+" found: "+seqOut);
			}*/
		}
		//System.out.println("First k-mer: "+kmers[0]);
		
		time2 = System.currentTimeMillis();
		System.out.println("Search time: "+ (time2 - time1));
		
	}
	

}
class ReadOverlap {
	private int indexSequence2;
	private int first1;
	private int last1;
	private int first2;
	private int last2;
	private boolean negativeStrand;
	public ReadOverlap(int indexSequence2, int first1, int last1, int first2, int last2, boolean negativeStrand) {
		super();
		this.indexSequence2 = indexSequence2;
		this.first1 = first1;
		this.last1 = last1;
		this.first2 = first2;
		this.last2 = last2;
		this.negativeStrand = negativeStrand;
	}
}
