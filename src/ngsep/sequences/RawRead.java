package ngsep.sequences;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class RawRead extends QualifiedSequence {
	public RawRead(String id, CharSequence sequence, String qualityScores) {
		super(id,sequence,qualityScores);
	}
	public static RawRead load(BufferedReader in) throws IOException {
		String id = in.readLine();
		if(id==null) return null;
		String seq = in.readLine();
		if(seq==null) return null;
		String plus = in.readLine();
		if(plus==null) return null;
		String qs = in.readLine();
		if(qs==null) return null;
		return new RawRead(id, seq, qs);
	}
	public String getSequenceString() {
		return getCharacters().toString();
	}
	public void trimFirstNucleotides (int n) {
		setCharacters(getSequenceString().substring(n));
		setQualityScores(getQualityScores().substring(n));
	}
	public void trimFromSequence(Pattern regexp) {
		String sequence = getSequenceString();
		if(regexp==null) return;
		Matcher m = regexp.matcher(sequence);
		if(m.find()) {
			int start = m.start();
			trimToLength(start);
		}
	}
	public void trimToLength(int newLength) {
		String sequence = getSequenceString();
		setCharacters(sequence.substring(0,newLength));
		setQualityScores(getQualityScores().substring(0,newLength));
	}
	public void save (PrintStream out) {
		out.println(this.getName());
		out.println(this.getCharacters());
		out.println("+");
		out.println(this.getQualityScores());
	}
}
