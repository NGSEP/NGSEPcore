package ngsep.assembly;

import static ngsep.assembly.TimeUtilities.timeIt;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.sequences.io.FastqFileReader;

public class GraphQualityAnalyzer {
	private static final String[] fastq = { ".fastq", ".fastq.gz" };
	private static final String[] fasta = { ".fasta", ".fa" };

	private SimplifiedAssemblyGraph lec;
	private SimplifiedAssemblyGraph ref;
	private List<Sequence> nams;

	public GraphQualityAnalyzer(String pathLects, OverlapConfiguration overlapConfiguration, boolean extrapolate)
			throws IOException {
		this.nams = load(pathLects);
		Collections.sort(nams, (l1, l2) -> l2.sequence.length() - l1.sequence.length());
		int i = 0;
		for (Sequence seq : nams)
			seq.id = i++;

		this.ref = getGraph(nams);

		AssemblyConfiguration ac = new AssemblyConfiguration();
		ac.setOverlap(overlapConfiguration);
		GraphBuilderFMIndex builder = new GraphBuilderFMIndex();
		builder.setConfig(ac);
		this.lec = builder
				.buildSimplifiedAssemblyGraph(this.nams.stream().map(a -> a.sequence).collect(Collectors.toList()));
		
		if (extrapolate)
			timeIt("      Extrapolate Aligns", () -> lec.ExtrapolateAligns());

		System.out.println("-------------PerfectGraph------------------");
		ref.printInfo();
		System.out.println("-------------CurrentGraph------------------");
		lec.printInfo();
		System.out.println("-------------------------------------------");
		emmbededTest();
	}

	public void emmbededTest() {
		Set<Integer> refEmb = new HashSet<>();
		for (Map<Integer, Embedded> a : ref.getEmbbeded().values())
			refEmb.addAll(a.keySet());

		Set<Integer> lectEmb = new HashSet<>();
		for (Entry<Integer, Map<Integer, Embedded>> a : lec.getEmbbeded().entrySet())
			lectEmb.addAll(a.getValue().keySet());

		int trueP = 0, falseP = 0;
		for (Entry<Integer, Map<Integer, Embedded>> a : lec.getEmbbeded().entrySet()) {
			for (Entry<Integer, Embedded> b : a.getValue().entrySet()) {
				int i = b.getKey();
				if (refEmb.contains(i))
					trueP++;
				else {
					int h = nams.get(a.getKey()).pos + nams.get(a.getKey()).len;
					int j = nams.get(i).pos + nams.get(i).len;
					int tope = (int) (ref.getSequences().get(i).length() * 0.01);
					if (Math.abs(h - j) <= tope || Math.abs(nams.get(a.getKey()).pos - nams.get(i).pos) <= tope) {
						trueP++;
					} else {
						falseP++;
						// System.out.println("falso embebido");
						// System.out.println("seq(" + i + ")" + " emb seq(" +
						// a.getKey() + ")");
						// System.out.println(b.getValue().getPos());
						// System.out.println(b.getValue().getRate());
						// System.out.println(b.getValue().isReversed());
						// System.out.println("--------------------------------");
						// System.out.println(nams.get(a.getKey()));
						// System.out.println(nams.get(i));
						// System.out.println("--------------------------------");
						// System.out.println(">ref");
						// System.out.println(ref.getSequences().get(a.getKey()));
						// System.out.println(">lect");
						// System.out.println(ref.getSequences().get(i));
					}
				}
			}

		}

		int falseN = 0;
		for (int i : refEmb) {
			if (!lectEmb.contains(i))
				falseN++;
		}
		int trueN = (ref.getSequences().size() - ref.amuontOfEmbeddedSequences()) - falseP;

		System.out.println("Embedded Sequences");
		System.out.println("false|true");
		System.out.println("neg= " + falseN + "|" + trueN);
		System.out.println("pos= " + falseP + "|" + trueP);
		System.out.println("precision = " + (100 * trueP) / (double) (trueP + falseP) + " %");
		System.out.println("recall = " + (100 * trueP) / (double) (trueP + falseN) + " %");

		Set<String> refEdg = new HashSet<>();
		for (Entry<Integer, Map<Integer, Alignment>> map : ref.getEdges().entrySet()) {
			int id1 = map.getKey();
			if (!lectEmb.contains(id1 >> 1) && !refEmb.contains(id1 >> 1))
				for (int id2 : map.getValue().keySet())
					if (id1 + 1 < id2 && !lectEmb.contains(id2 >> 1) && !refEmb.contains(id2 >> 1))
						refEdg.add(id1 + "-" + id2);
		}

		Set<String> lecEdg = new HashSet<>();
		for (Entry<Integer, Map<Integer, Alignment>> map : lec.getEdges().entrySet()) {
			int id1 = map.getKey();
			if (!lectEmb.contains(id1 >> 1) && !refEmb.contains(id1 >> 1))
				for (int id2 : map.getValue().keySet())
					if (id1 + 1 < id2 && !lectEmb.contains(id2 >> 1) && !refEmb.contains(id2 >> 1))
						lecEdg.add(id1 + "-" + id2);
		}

		System.out.println();
		System.out.println(refEdg.size());
		System.out.println(lecEdg.size());
		trueP = 0;
		falseP = 0;
		for (String i : lecEdg) {
			if (refEdg.contains(i))
				trueP++;
			else
				falseP++;
		}

		falseN = 0;
		for (String i : refEdg) {
			if (!lecEdg.contains(i))
				falseN++;
		}
		trueN = -1;

		System.out.println("edges (without embeddes in both graphs)");
		System.out.println("false|true");
		System.out.println("neg= " + falseN + "|" + trueN);
		System.out.println("pos= " + falseP + "|" + trueP);
		System.out.println("precision = " + (100 * trueP) / (double) (trueP + falseP) + " %");
		System.out.println("recall = " + (100 * trueP) / (double) (trueP + falseN) + " %");

		int s = 0, l = 0;
		for (Entry<Integer, Map<Integer, Alignment>> map : ref.getEdges().entrySet()) {
			int id1 = map.getKey();
			if (lectEmb.contains(id1) && refEmb.contains(id1))
				for (Entry<Integer, Alignment> a : map.getValue().entrySet()) {
					int id2 = a.getKey();
					if (id1 + 1 < id2 && lectEmb.contains(id2) && refEmb.contains(id2))
						if (!lecEdg.contains(id1 + "-" + id2)) {
							s += a.getValue().getOverlap();
							l++;
						}
				}
		}
		System.out.println("mean overlap (edgeFalsePositive) = " + (s / (double) l));

	}

	private SimplifiedAssemblyGraph getGraph(List<Sequence> sequences) throws FileNotFoundException {
		SimplifiedAssemblyGraph sag = new SimplifiedAssemblyGraph(getSequences(sequences));

		Map<String, List<Sequence>> a = new HashMap<>();
		for (Sequence s : sequences)
			a.computeIfAbsent(s.ref, (x) -> new ArrayList<>()).add(s);

		for (List<Sequence> list : a.values())
			Collections.sort(list, (Sequence x, Sequence y) -> {
				int ans = x.pos - y.pos;
				if (ans != 0)
					return ans;
				return y.len - x.len;
			});

		for (List<Sequence> list : a.values()) {
			for (int i = 0; i < list.size() - 1; i++) {
				Sequence lect1 = list.get(i);
				for (int j = i + 1; j < list.size() && list.get(j).pos < lect1.pos + lect1.len; j++) {
					Sequence lect2 = list.get(j);

					int relativePos = lect2.pos - lect1.pos;
					if (relativePos + lect2.len > lect1.len)
						// lect1 -> lect2
						sag.addEdge((lect1.id << 1) + (lect1.rev ? 0 : 1), (lect2.id << 1) + (lect2.rev ? 1 : 0),
								lect1.len - relativePos, 1);
					else {
						// lect2 into lect1
						boolean reversed = lect1.rev ^ lect2.rev;
						relativePos = (lect1.rev) ? lect1.len - lect2.len - relativePos : relativePos;
						sag.addEmbedded(lect1.id, lect2.id, relativePos, reversed, 1);
					}
				}
			}
		}

		sag.removeAllEmbeddedsIntoGraph();
		return sag;
	}

	private List<CharSequence> getSequences(List<Sequence> sequences) {
		return sequences.stream().map((Sequence a) -> a.sequence).collect(Collectors.toList());
	}

	/**
	 * Load the sequences of the file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	public List<Sequence> load(String filename) throws IOException {
		if (Stream.of(fastq)
				.anyMatch((String s) -> filename.endsWith(s.toLowerCase()) || filename.endsWith(s.toUpperCase()))) {
			return loadFastq(filename);
		} else if (Stream.of(fasta)
				.anyMatch((String s) -> filename.endsWith(s.toLowerCase()) || filename.endsWith(s.toUpperCase()))) {
			return loadFasta(filename);
		} else
			throw new IOException("the file not is a fasta or fastq file: " + filename);

	}

	/**
	 * Load the sequences of the Fasta file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	private List<Sequence> loadFasta(String filename) throws IOException {
		List<Sequence> sequences = new ArrayList<>();
		FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList seqsQl = handler.loadSequences(filename);
		int i = 0;
		for (QualifiedSequence seq : seqsQl) {
			DNAMaskedSequence characters = (DNAMaskedSequence) seq.getCharacters();
			String[] args = seq.getName().split("_");
			sequences.add(new Sequence(i++, args[args.length - 3], Integer.valueOf(args[args.length - 2]),
					(Integer.valueOf(args[args.length - 1]) == 1), characters));
		}
		return sequences;
	}

	/**
	 * Load the sequences of the Fastq file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	private List<Sequence> loadFastq(String filename) throws IOException {
		List<Sequence> sequences = new ArrayList<>();
		try (FastqFileReader reader = new FastqFileReader(filename)) {
			reader.setLoadMode(FastqFileReader.LOAD_MODE_FULL);
			reader.setSequenceType(DNAMaskedSequence.class);
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				DNAMaskedSequence characters = (DNAMaskedSequence) read.getCharacters();

				String[] args = read.getName().split("_");
				sequences.add(new Sequence(0, args[args.length - 3], Integer.valueOf(args[args.length - 2]),
						(Integer.valueOf(args[args.length - 1]) == 1), characters));
			}
		}
		return sequences;
	}

	public static void main(String[] args) throws FileNotFoundException, ClassNotFoundException, IOException {
		String pathLects = args[0];
		OverlapConfiguration overlapConfiguration = new OverlapConfiguration(Integer.valueOf(args[1]),
				Integer.valueOf(args[2]), Integer.valueOf(args[3]), Double.valueOf(args[4]));
		new GraphQualityAnalyzer(pathLects, overlapConfiguration, args.length >= 6 && args[5].equals("true"));
		// analizer.emmbededTest();
	}

	static class Sequence {
		int id;
		String ref;
		int pos;
		int len;
		boolean rev;
		CharSequence sequence;

		public Sequence(int id, String ref, int position, boolean isReversed, CharSequence sequence) {
			this.id = id;
			this.ref = ref;
			this.pos = position;
			this.len = sequence.length();
			this.rev = isReversed;
			this.sequence = sequence;
		}

		@Override
		public String toString() {
			return "Sequence [id=" + id + ", ref=" + ref + ", pos=" + pos + ", len=" + len + ", rev=" + rev + "]";
		}
	}
}
