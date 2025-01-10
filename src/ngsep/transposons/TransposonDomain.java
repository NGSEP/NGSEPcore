package ngsep.transposons;

import ngsep.hmm.ProfileAlignmentDomain;
import ngsep.sequences.QualifiedSequence;

public class TransposonDomain {
	private QualifiedSequence qseq;
	private int start;
	private int end;
	private boolean reverse = false;
	private ProfileAlignmentDomain alnDomain;
	public TransposonDomain(QualifiedSequence qseq, int start, int end, ProfileAlignmentDomain alnDomain) {
		this.qseq = qseq;
		this.start = start;
		this.end = end;
		this.alnDomain = alnDomain;
	}
	public boolean isReverse() {
		return reverse;
	}
	public void setReverse(boolean reverse) {
		this.reverse = reverse;
	}
	public QualifiedSequence getQseq() {
		return qseq;
	}
	public int getStart() {
		return start;
	}
	public int getEnd() {
		return end;
	}
	public int getLength() {
		return end-start;
	}
	public ProfileAlignmentDomain getAlnDomain() {
		return alnDomain;
	}
	
	
	
}
