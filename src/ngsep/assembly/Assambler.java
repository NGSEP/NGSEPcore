package ngsep.assembly;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Stream;
import ngsep.alignments.ReadAlignment;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.sequences.io.FastqFileReader;

public class Assambler {
	/**
	 * represent the maximum number of non kmers between two kmers of the same hit
	 */
	private static final int NUMBER_INVALID_KMERS = 5;
	/** represent the maximum bleed + 1 between two kmers of the same hit */
	private static final int MAX_KMER_BLEED = 6;
	/** Possible end of Fastq file */
	private static final String[] FASTQ = { ".fastq", ".fastq.gz" };
	/** Possible end of Fasta file */
	private static final String[] FASTA = { ".fasta" };
	/** format to print the time */
	private static final DecimalFormat FORMAT = new DecimalFormat("0.00000");
	/** manage the hits of the kmers of each sequence */
	private final TreeMap<Integer, int[]>[][] hits;

	private Map<Integer, ReadOverlap> embeddedOverlaps = new Hashtable<>();

	@SuppressWarnings("unchecked")
	public Assambler(final Logger log, String filePath, String path) throws IOException {
		File file = new File(path);
		file.mkdir();

		// load the file
		System.out.println("Loading " + filePath + " ...");
		long ini = System.currentTimeMillis();
		List<DNAMaskedSequence> sequences = load(filePath);
		System.out.println("Time " + FORMAT.format((System.currentTimeMillis() - ini) / (double) 1000) + " s");

		// sort the reads
		System.out.println("Sort the reads " + filePath + " ...");
		ini = System.currentTimeMillis();
		Collections.sort(sequences, new Comparator<DNAMaskedSequence>() {
			public int compare(DNAMaskedSequence arg0, DNAMaskedSequence arg1) {
				return arg0.length() - arg1.length();
			}
		});
		System.out.println("Time " + FORMAT.format((System.currentTimeMillis() - ini) / (double) 1000) + " s");

		// create the hits manager
		hits = new TreeMap[NUMBER_INVALID_KMERS + 1][sequences.size()];
		for (int y = 0; y <= NUMBER_INVALID_KMERS; y++)
			for (int x = 0; x < sequences.size(); x++)
				hits[y][x] = new TreeMap<>();

		// First fmIndex
		System.out.println("          Making the firts FM index");
		ini = System.currentTimeMillis();
		FMIndex fmIndex = fmIndex(sequences);
		System.out.println("Time " + FORMAT.format((System.currentTimeMillis() - ini) / (double) 1000) + " s");

		// Identify the embed sequences
		System.out.println("          Identify the embed sequences positive strand");
		ini = System.currentTimeMillis();
		emmbeddeDetection(sequences, fmIndex);
		System.out.println("Time " + FORMAT.format((System.currentTimeMillis() - ini) / (double) 1000) + " s");

	}

	private void emmbeddeDetection(List<DNAMaskedSequence> sequencesToIndex, FMIndex index) {
		int n;
		int k;
		Iterator<String> iter;
		TreeMap<Integer, int[]> tree;
		boolean[] mark = new boolean[sequencesToIndex.size()];
		Hashtable<Integer, List<int[]>> hitsPoint = new Hashtable<>();
		for (int idSequence = 0; idSequence < sequencesToIndex.size(); idSequence++) {
			System.out.println(idSequence);
			hitsPoint.clear();
			KmerEmmbeddedIterator ki = new KmerEmmbeddedIterator(sequencesToIndex.get(idSequence));
			n = (int) Math.round(ki.getNumber() * KmerEmmbeddedIterator.PORCENTAJE_OF_TOLERANCE);
			if (n <= 0)
				n = 1;

			k = 0;
			iter = ki.firts().iterator();
			for (int i = 0; i < n; i++) {
				for (ReadAlignment aln : index.search(iter.next(), idSequence + 1, sequencesToIndex.size()))
					alignAdd(aln, k);
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}
			for (int i = idSequence + 1; i < mark.length; i++) {
				boolean a = false;
				for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
					if (!hits[j][i].isEmpty()) {
						a = true;
						break;
					}
				mark[i] = a;
			}
			while (iter.hasNext()) {
				for (ReadAlignment aln : index.search(iter.next(), idSequence + 1, sequencesToIndex.size())) {
					if (mark[Integer.parseInt(aln.getSequenceName())])
						align(aln, k);
				}
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}
			for (int i = idSequence + 1; i < mark.length; i++) {
				if (mark[i]) {
					boolean a = false;
					for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
						if (!hits[j][i].isEmpty()) {
							for (int[] array : hits[j][i].values()) {
								if (array[0] + sequencesToIndex.get(idSequence).length() <= sequencesToIndex.get(i)
										.length()) {
									a = true;
									hitsPoint.computeIfAbsent(i, (o) -> new LinkedList<>()).add(array);
								}
							}
							hits[j][i].clear();
						}
					mark[i] = a;
				}
			}
			k = 0;
			iter = ki.firts().iterator();
			for (int i = 0; i < n; i++) {
				for (ReadAlignment aln : index.search(iter.next(), idSequence + 1, sequencesToIndex.size()))
					alignAdd(aln, k);
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}
			for (int i = idSequence + 1; i < mark.length; i++) {
				boolean a = false;
				for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
					if (!hits[j][i].isEmpty()) {
						a = true;
						break;
					}
				mark[i] = a;
			}
			while (iter.hasNext()) {
				for (ReadAlignment aln : index.search(iter.next(), idSequence + 1, sequencesToIndex.size())) {
					if (mark[Integer.parseInt(aln.getSequenceName())])
						align(aln, k);
				}
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}
			for (int i = idSequence + 1; i < mark.length; i++) {
				if (mark[i]) {
					boolean a = false;
					for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
						if (!hits[j][i].isEmpty()) {
							for (int[] array : hits[j][i].values()) {
								if (array[0] + sequencesToIndex.get(idSequence).length() <= sequencesToIndex.get(i)
										.length()) {
									a = true;
									hitsPoint.computeIfAbsent(i, (o) -> new LinkedList<>()).add(array);
								}
							}
							hits[j][i].clear();
						}
					mark[i] = a;
				}
			}

			int falt = ki.getNumber() * KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
			k = sequencesToIndex.get(idSequence).length() - falt;
			iter = ki.lasts().iterator();
			for (int i = 0; i < n; i++) {
				for (ReadAlignment aln : index.search(iter.next(), idSequence + 1, sequencesToIndex.size()))
					if (mark[Integer.parseInt(aln.getSequenceName())]) {
						if (aln.getFirst() + falt <= sequencesToIndex.get(Integer.parseInt(aln.getSequenceName()))
								.length())
							alignAdd(aln, k);
					}
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				falt -= KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}
			for (int i = idSequence + 1; i < mark.length; i++) {
				if (mark[i]) {
					boolean a = false;
					for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
						if (!hits[j][i].isEmpty()) {
							a = true;
							break;
						}
					mark[i] = a;
				}
			}
			while (iter.hasNext()) {
				for (ReadAlignment aln : index.search(iter.next(), idSequence + 1, sequencesToIndex.size())) {
					if (mark[Integer.parseInt(aln.getSequenceName())])
						align(aln, k);
				}
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}

			for (int i = idSequence + 1; i < mark.length; i++) {
				if (mark[i]) {
					for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
						if (!hits[j][i].isEmpty()) {
							for (int[] array : hits[j][i].values()) {
								for (int[] jj : hitsPoint.get(i))
									if (Math.abs(array[0] + (ki.getNumber() * ki.SEARCH_KMER_LENGTH) - jj[0]
											- sequencesToIndex.get(idSequence).length()) < sequencesToIndex
													.get(idSequence).length() / (double) 10) {
										embeddedOverlaps.put(idSequence,
												new ReadOverlap(idSequence, 0,
														sequencesToIndex.get(idSequence).length() - 1, i, jj[0],
														array[0] + (ki.getNumber() * ki.SEARCH_KMER_LENGTH), false));
										break;
									}
							}
							hits[j][i].clear();
						}
				}
			}
		}
	}

	/**
	 * try to align the read whit each possible aligns in a
	 * 
	 * @param aln            the read (one alignment for the kmer)
	 * @param posinRefrenece index of the kmer in the original sequence.
	 */
	private void alignAdd(ReadAlignment aln, int posinRefrenece) {
		int idSequenceAligned = Integer.parseInt(aln.getSequenceName());
		int pos = aln.getFirst();

		int minposkmer = 0;
		int minnode = 0;
		int min = MAX_KMER_BLEED;
		for (int dif = KmerEmmbeddedIterator.SEARCH_KMER_LENGTH, i = 1; i <= NUMBER_INVALID_KMERS; i++, dif += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH) {
			TreeMap<Integer, int[]> treeMap = hits[i][idSequenceAligned];
			if (treeMap.isEmpty())
				continue;

			int t = pos - dif;
			Integer p1 = treeMap.ceilingKey(t);
			Integer p2 = treeMap.floorKey(t);
			if (p1 != null && Math.abs(p1 - t) < min) {
				min = Math.abs(p1 - t);
				minposkmer = i;
				minnode = p1;
			}
			if (p2 != null && Math.abs(p2 - t) < min) {
				min = Math.abs(p2 - t);
				minposkmer = i;
				minnode = p2;
			}
		}

		if (min == MAX_KMER_BLEED) {
			hits[0][idSequenceAligned].put(pos, new int[] { pos, posinRefrenece });
		} else {
			int[] aux = hits[minposkmer][idSequenceAligned].remove(minnode);
			hits[0][idSequenceAligned].put(pos, aux);
		}
	}

	/**
	 * try to align the read whit each possible aligns in a
	 * 
	 * @param aln            the read (one alignment for the kmer)
	 * @param posinRefrenece index of the kmer in the original sequence.
	 */
	private void align(ReadAlignment aln, int posinRefrenece) {
		int idSequenceAligned = Integer.parseInt(aln.getSequenceName());
		int pos = aln.getFirst();

		int minposkmer = 0;
		int minnode = 0;
		int min = MAX_KMER_BLEED;
		for (int dif = KmerEmmbeddedIterator.SEARCH_KMER_LENGTH, i = 1; i <= NUMBER_INVALID_KMERS; i++, dif += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH) {
			TreeMap<Integer, int[]> treeMap = hits[i][idSequenceAligned];
			if (treeMap.isEmpty())
				continue;

			int t = pos - dif;
			Integer p1 = treeMap.ceilingKey(t);
			Integer p2 = treeMap.floorKey(t);
			if (p1 != null && Math.abs(p1 - t) < min) {
				min = Math.abs(p1 - t);
				minposkmer = i;
				minnode = p1;
			}
			if (p2 != null && Math.abs(p2 - t) < min) {
				min = Math.abs(p2 - t);
				minposkmer = i;
				minnode = p2;
			}
		}

		if (min != MAX_KMER_BLEED) {
			int[] aux = hits[minposkmer][idSequenceAligned].remove(minnode);
			hits[0][idSequenceAligned].put(pos, aux);
		}
	}

	/**
	 * rotate the columns of hits
	 */
	private void rotate() {
		TreeMap<Integer, int[]>[] temp = hits[NUMBER_INVALID_KMERS];
		for (TreeMap<Integer, int[]> t : temp)
			t.clear();
		System.arraycopy(hits, 0, hits, 1, NUMBER_INVALID_KMERS);
		hits[0] = temp;
	}

	/**
	 * Create a fmIndex for the sequences
	 * 
	 * @param sequences the sequences
	 * @return the fmIndex
	 */
	private static FMIndex fmIndex(List<DNAMaskedSequence> sequences) {
		FMIndex index = new FMIndex();
		index.loadUnnamedSequences(sequences, 100, 25);
		return index;
	}

	/**
	 * Load the sequences of the file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	private static List<DNAMaskedSequence> load(String filename) throws IOException {
		if (Stream.of(FASTQ)
				.anyMatch((String s) -> filename.endsWith(s.toLowerCase()) || filename.endsWith(s.toUpperCase())))
			return loadFastq(filename);
		else if (Stream.of(FASTA)
				.anyMatch((String s) -> filename.endsWith(s.toLowerCase()) || filename.endsWith(s.toUpperCase())))
			return loadFasta(filename);

		throw new IOException("the file not is a fasta or fastq file: " + filename);
	}

	/**
	 * Load the sequences of the Fasta file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	private static List<DNAMaskedSequence> loadFasta(String filename) throws IOException {
		List<DNAMaskedSequence> sequences = new ArrayList<>();
		FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList seqsQl = handler.loadSequences(filename);
		for (QualifiedSequence seq : seqsQl) {
			DNAMaskedSequence characters = (DNAMaskedSequence) seq.getCharacters();
			sequences.add(characters);
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
	private static List<DNAMaskedSequence> loadFastq(String filename) throws IOException {
		List<DNAMaskedSequence> sequences = new ArrayList<>();
		try (FastqFileReader reader = new FastqFileReader(filename)) {
			reader.setLoadMode(FastqFileReader.LOAD_MODE_MINIMAL);
			reader.setSequenceType(DNAMaskedSequence.class);
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				DNAMaskedSequence characters = (DNAMaskedSequence) read.getCharacters();
				sequences.add(characters);
			}
		}
		return sequences;
	}

	public static void main(String[] args) {
		Logger log = Logger.getLogger(Assambler.class.getSimpleName());
		switch (args.length) {
		case 0:
			log.log(Level.SEVERE, "please enter the file path and the path to answer");
			break;

		case 1:
			log.log(Level.SEVERE, "please enter path to answer");
			break;

		case 2:
			try {
				new Assambler(log, args[0], args[1]);
			} catch (Exception e) {
				log.log(Level.SEVERE, "Error", e);
			}
			break;

		default:
			log.log(Level.SEVERE, "too many arguments");
			break;
		}
	}
}
