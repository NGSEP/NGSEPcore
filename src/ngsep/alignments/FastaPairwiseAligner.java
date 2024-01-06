package ngsep.alignments;

import java.util.List;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;

public class FastaPairwiseAligner {

	public static void main(String[] args) throws Exception {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> seqs= handler.loadSequences(args[0]);
		PairwiseAligner aligner = new PairwiseAlignerDynamicKmers();
		//PairwiseAligner aligner = new PairwiseAlignerStaticBanded();
		CharSequence seq1 = seqs.get(0).getCharacters();
		CharSequence seq2 = seqs.get(1).getCharacters();
		System.out.println("Length1: "+seq1.length()+" length2: "+seq2.length());
		String [] alignedSeqs = aligner.calculateAlignment(seq1, seq2 );
		System.out.println(alignedSeqs[0]);
		System.out.println(alignedSeqs[1]);
		

	}

}
