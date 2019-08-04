package ngsep.assembly;

import java.util.Iterator;
import java.util.AbstractMap.SimpleEntry;
import java.util.Map.Entry;

import ngsep.sequences.DNAMaskedSequence;

public class KmerIterator {

	private AssemblyConfiguration config;

	public KmerIterator(AssemblyConfiguration config) {
		this.config = config;
	}

	public Iterable<Entry<Integer, String>> positiveStrand(CharSequence seq) {
		String str = seq.toString();

		return new Iterable<Entry<Integer, String>>() {
			public Iterator<Entry<Integer, String>> iterator() {
				return new Iterator<Entry<Integer, String>>() {
					int i = 0;

					public boolean hasNext() {
						return i < str.length() - config.overlap().getKmerLength();
					}

					public Entry<Integer, String> next() {
						String str2 = str.substring(i, i + config.overlap().getKmerLength());
						Entry<Integer, String> ans = new SimpleEntry<>(i, str2);
						i += config.overlap().getKmerDistance() + config.overlap().getKmerLength();
						return ans;
					}
				};
			};
		};
	}

	public Iterable<Entry<Integer, String>> negativeStrand(CharSequence seq) {
		String str = DNAMaskedSequence.getReverseComplement(seq).toString();
		return new Iterable<Entry<Integer, String>>() {
			public Iterator<Entry<Integer, String>> iterator() {
				return new Iterator<Entry<Integer, String>>() {
					int i = 0;

					public boolean hasNext() {
						return i < str.length() - config.overlap().getKmerLength();
					}

					public Entry<Integer, String> next() {
						String str2 = str.substring(i, i + config.overlap().getKmerLength());
						Entry<Integer, String> ans = new SimpleEntry<>(i, str2);
						i += config.overlap().getKmerDistance() + config.overlap().getKmerLength();
						return ans;
					}
				};
			};
		};
	}

}
