package ngsep.transposons;

import ngsep.sequences.QualifiedSequence;

public class TransposableElement {
	private String id;
	private CharSequence sequence;
	private TransposableElementFamily family;
	private String taxonomy;
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
}
