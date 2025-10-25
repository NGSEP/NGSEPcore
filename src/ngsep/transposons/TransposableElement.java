package ngsep.transposons;

import java.util.List;

import ngsep.alignments.PairwiseAlignerSimpleGap;
import ngsep.alignments.PairwiseAlignment;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;

public class TransposableElement {
	private String id;
	private CharSequence sequence;
	private TransposableElementFamily family;
	private String taxonomy;
	private List<TransposonDomainAlignment> domainAlignments;
	private int leftEndRepeat = -1;
	private int rightStartRepeat = -1;
	private byte orientation = -1;
	private String tsd = null;
	private static final int ENDS_LENGTH_ALIGNMENT = 1200;
	
	public TransposableElement(String id, CharSequence sequence) {
		super();
		this.id = id;
		this.sequence = sequence;
		decodeSequenceId();
	}
	public TransposableElement(QualifiedSequence seq) {
		this(seq.getName(),seq.getCharacters());
	}
	public TransposableElementFamily getFamily() {
		return family;
	}
	public void setFamily(TransposableElementFamily family) {
		this.family = family;
	}
	
	public String getTaxonomy() {
		return taxonomy;
	}
	public String getId() {
		return id;
	}
	public CharSequence getSequence() {
		return sequence;
	}
	
	public List<TransposonDomainAlignment> getDomainAlignments() {
		return domainAlignments;
	}
	public void setDomainAlignments(List<TransposonDomainAlignment> domainAlignments) {
		this.domainAlignments = domainAlignments;
	}
	public int getLeftEndRepeat() {
		return leftEndRepeat;
	}
	public int getRightStartRepeat() {
		return rightStartRepeat;
	}
	public byte getOrientation() {
		return orientation;
	}
	public void setRepeatLimits(int leftEnd, int rightStart, byte orientation) {
		this.leftEndRepeat = leftEnd;
		this.rightStartRepeat = rightStart;
		this.orientation = orientation;
	}
	
	public String getTsd() {
		return tsd;
	}
	public void setTsd(String tsd) {
		this.tsd = tsd;
	}
	public boolean isValidated() {
		return leftEndRepeat>0 && getFamily()!=null;
	}
	
	private void decodeSequenceId() {
		int i = id.indexOf('#');
		if(i<0) return;
		taxonomy = id.substring(i+1);
		decodeTaxonomy();
	}
	private void decodeTaxonomy() {
		int i;
		i=taxonomy.indexOf('/');
		if(i<0) i=taxonomy.length();
		String orderStr = taxonomy.substring(0,i);
		if (i<taxonomy.length()) {
			String familyInfo2 = taxonomy.substring(i+1);
			i=familyInfo2.indexOf('/');
			if(i<0) i=familyInfo2.length();
			String familyStr = familyInfo2.substring(0,i);
			family = TransposableElementFamily.findFamily(orderStr, familyStr);
		} else {
			family = TransposableElementFamily.findUnknown(orderStr);
		}
	}
	public void modifyIdFromFamily () {
		StringBuilder newTaxonomy = new StringBuilder();
		int i = id.indexOf('#');
		if(i>0) id = id.substring(0,i);
		newTaxonomy.append(family.getOrder());
		newTaxonomy.append('/');
		newTaxonomy.append(family.getId());
		taxonomy = newTaxonomy.toString();
		id = id + '#'+taxonomy;
	}
	
	public int [] alignEnds() {
		PairwiseAlignment alignment = findEndsAlignment(false);
		boolean reverse = (alignment == null);
		if (reverse) alignment = findEndsAlignment(true);
		if (alignment == null) return null;
		return getBordersInSequence(alignment, reverse);
	}
	private PairwiseAlignment findEndsAlignment(boolean reverse) {
		int length = sequence.length();
		
		if(length< 2*ENDS_LENGTH_ALIGNMENT) return null; 
		//if(family == null) return null;
		//if(!family.isLTR()) return null;
		PairwiseAlignerSimpleGap pairAligner = new PairwiseAlignerSimpleGap();
		pairAligner.setLocal(true);
		
		CharSequence leftSeq = sequence.subSequence(0, ENDS_LENGTH_ALIGNMENT);
		CharSequence rightSeq = sequence.subSequence(length - ENDS_LENGTH_ALIGNMENT, length);
		if(reverse) rightSeq = DNAMaskedSequence.getReverseComplement(rightSeq);
		PairwiseAlignment alignment =  pairAligner.calculateAlignment(leftSeq, rightSeq);
		int [] borders = getBordersInSequence(alignment, reverse);
		int length1 = borders[1] - borders[0];
		int length2 = borders[3] - borders[2];
		int minLength = Math.min(length1, length2);
		int maxLength = Math.max(length1, length2);
		System.out.println("End alignment. TE length: "+length+" Lengths: "+minLength+" "+maxLength+". Borders1: "+borders[0]+" "+borders[1]+" Borders2: "+borders[2]+" "+borders[3]+". Mismatches: "+alignment.getMismatches());
		if (minLength > 0.9 * maxLength && minLength > 100 && alignment.getMismatches() < 0.2 * minLength) {
			return alignment;
		}
		return null;
	}
	private int[] getBordersInSequence(PairwiseAlignment alignment, boolean reverse) {
		int length = sequence.length();
		int start1 = alignment.getStart1();
		int end1 = alignment.getEnd1();
		int start2 = alignment.getStart2() + length - ENDS_LENGTH_ALIGNMENT;
		int end2 = alignment.getEnd2() + length - ENDS_LENGTH_ALIGNMENT;
		int revN = (reverse?1:0);
		int [] answer = {start1,end1,start2,end2,revN};
		return answer;
	}
}
