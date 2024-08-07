package ngsep.alignments;

import java.util.List;

import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;

public class FastaPairwiseAligner {

	public static void main(String[] args) throws Exception {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> seqs= handler.loadSequences(args[0]);
		//PairwiseAligner aligner = new PairwiseAlignerDynamicKmers();
		PairwiseAligner aligner = new PairwiseAlignerStaticBanded();
		//((PairwiseAlignerStaticBanded)aligner).setForceEnd2(false);
		//PairwiseAligner aligner = new PairwiseAlignerSimpleGap();
		CharSequence seq1 = seqs.get(0).getCharacters();
		CharSequence seq2 = seqs.get(1).getCharacters();
		System.out.println("Length1: "+seq1.length()+" length2: "+seq2.length());
		PairwiseAlignment aln = aligner.calculateAlignment(seq1, seq2 );
		System.out.println(aln.getAlignedSequence1());
		System.out.println(aln.getAlignedSequence2());
		double numMismatches=(new HammingSequenceDistanceMeasure()).calculateDistance(aln.getAlignedSequence1(), aln.getAlignedSequence2());
		System.out.println("Mismatches: "+numMismatches);
		

	}

}
