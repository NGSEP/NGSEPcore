package ngsep.transposons;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.ReadAlignment.Platform;
import ngsep.alignments.ReadsAligner;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class DeNovoTransposableElementsFinderConservedEnds implements DeNovoTransposableElementsFinder {

	private int debugPos = -1;
	private HMMTransposonDomainsFinder domainsFinder = new HMMTransposonDomainsFinder();
	
	public DeNovoTransposableElementsFinderConservedEnds () {
		domainsFinder.loadHMMsFromClasspath();
	}
	
	@Override
	public List<TransposableElementAnnotation> findTransposons(ReferenceGenome genome) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		ReadsAligner aligner = new ReadsAligner(genome,Platform.PACBIO);
		aligner.setMaxAlnsPerRead(100);
		QualifiedSequenceList seqs = genome.getSequencesList();
		for(QualifiedSequence seq:seqs) {
			DNAMaskedSequence dna = (DNAMaskedSequence)seq.getCharacters();
			int n = dna.length();
			int subseqLength = 500;
			int step = 200;
			for(int i=0;i<n-subseqLength;i+=step) {
				DNAMaskedSequence read = (DNAMaskedSequence)dna.subSequence(i,i+subseqLength);
				List<ReadAlignment> alns = aligner.alignRead(new QualifiedSequence("",read));
				List<ReadAlignment> alnsSeq = selectAlnsInSequence(seq.getName(),alns);
				ReadAlignment aln = selectCloseAlignment(i,subseqLength,alnsSeq);
				if(aln==null) continue;
				//Discard possible hits due to tandem repeats
				int countB = countAlnsInBetween(i+subseqLength,aln.getFirst(),alnsSeq);
				if(countB>1) continue;
				TransposableElementAnnotation ltrAnn = inferTEFromEndAlignment(seq,i,aln);
				if(ltrAnn!=null) {
					answer.add(ltrAnn);
					i=ltrAnn.getLast()-step;
				}
			}
		}
		return answer;
	}

	private List<ReadAlignment> selectAlnsInSequence(String name, List<ReadAlignment> alns) {
		List<ReadAlignment> alnsSeq = new ArrayList<ReadAlignment>();
		for(ReadAlignment aln:alns) {
			if(name == aln.getSequenceName()) alnsSeq.add(aln);
		}
		return alnsSeq;
	}

	private ReadAlignment selectCloseAlignment(int start, int length,  List<ReadAlignment> alns) {
		ReadAlignment bestAln = null;
		for(ReadAlignment aln:alns) {
			if(aln.getFirst()<start+length) continue;
			if(aln.getFirst()>start+30000) continue;
			if(improvedAln(start, bestAln, aln)) bestAln = aln;
		}
		return bestAln;
	}

	private boolean improvedAln(int start, ReadAlignment bestAln, ReadAlignment aln) {
		if(bestAln==null) return true;
		int distanceBest = bestAln.getFirst()-start;
		int distanceAln = aln.getFirst()-start;
		int diff5000Best = Math.abs(distanceBest-5000);
		int diff5000Aln = Math.abs(distanceAln-5000);
		return diff5000Aln<diff5000Best;
	}
	
	private int countAlnsInBetween(int start, int end, List<ReadAlignment> alns) {
		int count = 0;
		for(ReadAlignment aln:alns) {
			if(start<=aln.getFirst() && aln.getLast()<=end) count++; 
		}
		return count;
	}

	private TransposableElementAnnotation inferTEFromEndAlignment(QualifiedSequence seq, int start, ReadAlignment aln) {
		int start2 = Math.max(0, start-500);
		int end2 = Math.min(seq.getLength(), aln.getLast()+500);
		QualifiedSequence rawSeq = new QualifiedSequence("",seq.getCharacters().subSequence(start2, end2));
		System.err.println("Checking hit of seq from "+start2+" seqLength: "+rawSeq.getLength()+" last: "+end2+" end aln: "+aln);
		TransposableElement te = new TransposableElement(rawSeq);
		int [] ends = te.alignEnds();
		if(ends == null) return null;
		System.err.println("ends: "+ends[0]+" "+ends[1]+" "+ends[2]+" "+ends[3]+" "+ends[4]);
		TransposableElementAnnotation ann = new TransposableElementAnnotation(seq.getName(), start+ends[0]+1, start+ends[3]);
		
		domainsFinder.assignFamily(ann,(DNAMaskedSequence)rawSeq.getCharacters());
		if(ann.getInferredFamily()==null) {
			domainsFinder.assignFamily(ann,((DNAMaskedSequence)rawSeq.getCharacters()).getReverseComplement());
			if(ann.getInferredFamily()!=null) ann.setNegativeStrand(true);
		}
		ann.setBordersFixed(true);
		return ann;
	}
	public static void main(String[] args) throws Exception {
		ReferenceGenome genome = new ReferenceGenome(args[0]);
		DeNovoTransposableElementsFinderConservedEnds instance = new DeNovoTransposableElementsFinderConservedEnds();
		List<TransposableElementAnnotation> anns = instance.findTransposons(genome);
		try (PrintStream out=new PrintStream(args[1])) {
			for(TransposableElementAnnotation ann:anns) out.println(""+ann.getSequenceName()+"\t"+ann.getFirst()+"\t"+ann.getLast()+"\t"+ann.isPositiveStrand()+"\t"+ann.getInferredFamily());
		}	
	}
}
