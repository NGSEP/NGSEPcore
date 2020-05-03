package ngsep.alignments.test;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import junit.framework.TestCase;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.FMIndexReadAlignmentAlgorithm;
import ngsep.genome.GenomeIndexer;
import ngsep.genome.GenomicRegion;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.ReferenceGenomeFMIndex;

public class ShortSingleReadsAlignerTest extends TestCase {
	private FMIndexReadAlignmentAlgorithm readsAligner;
	public final static String FM_INDEX_PATH= ".\\test\\Saccharomyces_cerevisiae.fmindex";
	public final static String FASTA_PATH= ".\\training\\Saccharomyces_cerevisiae.fa";


	public void setUpReadsAligner() throws IOException {
		ReferenceGenome genome = new ReferenceGenome(FASTA_PATH);
		File f = new File(FM_INDEX_PATH);
		if(!f.exists()) {
			GenomeIndexer genomeIndexer=new GenomeIndexer();
			genomeIndexer.createIndex(FASTA_PATH,FM_INDEX_PATH);	
		}
		ReferenceGenomeFMIndex fmIndex = ReferenceGenomeFMIndex.load(genome, FM_INDEX_PATH);
		readsAligner=new FMIndexReadAlignmentAlgorithm(fmIndex, 15, 3);
	}

	public void afterSetUpReadsAligner() {
		File f = new File(FM_INDEX_PATH);
		f.delete();
	}

	public void testLoadTRF() throws IOException {
		setUpReadsAligner();
		String path = setUpTRF();
		readsAligner.loadSTRsFile(path);
		Map<String, List<GenomicRegion>> map = readsAligner.getKnownSTRs(); 
		assertEquals(false, isOverlappging(map));
		assertEquals(true, map.get("chrI").size()<=3);
		assertEquals(true, map.get("chrII").size()<=3);

		//Case 1
		//Region "chrXI 122235 122272"
		ReadAlignment aln=new ReadAlignment("chrXI", 122176, 122265, 90, 0);
		String read =    "TCGGATCGAAATGAACGATATTCTCCCTTATATTATCAGCCAGTAGCGTATACTCTGGCATTTTTCATTTATTTGACTTATTTTTATTTN";
		//	122176	122235 -> 122176	122234
		GenomicRegion region =readsAligner.findTandemRepeat(aln.getSequenceName(),aln.getFirst(),aln.getLast());
		ReadAlignment newAln=readsAligner.verifyShortTandemRepeats(aln.getSequenceName(),aln.getFirst(),aln.getLast(),read,region);
		assertEquals(122176, newAln.getFirst());
		assertEquals(122234, newAln.getLast());
		assertEquals("59M31S", newAln.getCigarString());

		//Case 2
		//Region "chrII 151275 151285"
		aln=new ReadAlignment("chrII", 151281, 151370, 90, 0);
		read =          "TTTTTATTATGCATTTAAGAGTAGTCTCTACTTATGAACATTTTCTCTGGCCTCTGATCACGTTACTTTATTACCCGGATACTGATCATN";
		//	151281	151370 -> 151286	151370
		region =readsAligner.findTandemRepeat(aln.getSequenceName(),aln.getFirst(),aln.getLast());
		newAln=readsAligner.verifyShortTandemRepeats(aln.getSequenceName(),aln.getFirst(),aln.getLast(),read,region);
		assertEquals(151286, newAln.getFirst());
		assertEquals(151370, newAln.getLast());
		assertEquals("5S85M", newAln.getCigarString());

		//Case 3a
		//Region "chrIX 255901 255918"
		aln=new ReadAlignment("chrIX", 255867, 255956, 90, 0);
		read =         "ATTTTCTTTTATTTTTTTGATAAAACTACTACGCTAAAAATAAAATAAAAATGTATGATTTCCCTCCATTTCCGACCAATTGTATAATTT";
		//Left:  255867 255900
		//Right: 255919 255956
		region =readsAligner.findTandemRepeat(aln.getSequenceName(),aln.getFirst(),aln.getLast());
		newAln=readsAligner.verifyShortTandemRepeats(aln.getSequenceName(),aln.getFirst(),aln.getLast(),read,region);
		assertNotNull(newAln);
		assertEquals(255867, newAln.getFirst());
		assertEquals(255956, newAln.getLast());
		assertEquals("34M18M38M", newAln.getCigarString());

		//Case 3c
		//Region "chrXII 460003 460019"
		aln=new ReadAlignment("chrXII", 459953, 460042, 90, 0);
		read =         "GATATGTACAAACAATATCCTCCTCCGATATTCTCCTCCGATATTCCTACAAAAAAAAAAACACTCCGGTTTTGTTCTCTTCCCTCCATT";
		//Left:  459953 460002
		//Right: 460020 460042
		region =readsAligner.findTandemRepeat(aln.getSequenceName(),aln.getFirst(),aln.getLast());
		newAln=readsAligner.verifyShortTandemRepeats(aln.getSequenceName(),aln.getFirst(),aln.getLast(),read,region);
		assertNotNull(newAln);
		assertEquals(459966, newAln.getFirst());
		assertEquals(460042, newAln.getLast());
		assertEquals("13M3I1M10I23M17M6D14M2I1M1I2M3I", newAln.getCigarString());

		File f = new File(path);
		f.delete();
		afterSetUpReadsAligner();
	}

	private Object isOverlappging(Map<String, List<GenomicRegion>> map) {
		Iterator<String>it=map.keySet().iterator();
		while(it.hasNext()) {
			if(isOverlappging(map.get(it.next()))) {
				return true;
			}
		}
		return false;
	}

	private boolean isOverlappging(List<GenomicRegion>l) {
		for (int i = 0; i < l.size()-1; i++) {
			GenomicRegion current = l.get(i);
			GenomicRegion next = l.get(i+1);
			if(isOverlappging(current, next)) {
				return true;
			}
		}
		return false;
	}

	private boolean isOverlappging(GenomicRegion a,GenomicRegion b) {
		return b.getFirst()>=a.getFirst()&&b.getFirst()<=a.getLast();
	}

	private String setUpTRF() {
		List<String> lines = Arrays.asList(
				"chrXII 460003 460019",
				"chrIX 255901 255918",
				"chrII 151275 151285",
				"chrXI 122235 122272",
				"chrI 2 62 2 33.0 2 75 15 60 40 59 0 0 0.98 CA CACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACA C TCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCT",
				"chrI 2 61 12 5.4 12 84 15 85 40 60 0 0 0.97 CACACCACACCA CACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACAC C ATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTC",
				"chrI 255 277 11 2.1 11 100 0 46 26 43 0 30 1.55 TTACCCTACCA TTACCCTACCATTACCCTACCAT CCACTCACCCACCGTTACCCTCCAATTACCCATATCCAACCCACTGCCAC CCACCATGACCTACTCACCATACTGTTCTTCTACCCACCATATTGAAACG",
				"chrII 1524 1575 12 4.3 12 65 18 36 15 11 42 30 1.82 TGGTGCTAGCAG TGGTGCTAGCAGTGGTAGTGGCATTAGTGCTGGAGTTGGTGCTAGCAGTGGT CTTGGCACTAGCGTTGGTACTTTCAGTGGTAGTGGCATTAGTGCTGGAGT AGTAGCACTAGTGTTGGAGTCGGTACTTTCGGTGGTAGTAGCACTAGTGT",
				"chrII 1883 2014 24 5.5 24 70 13 101 17 14 37 31 1.90 TTGGTGCTGGCAGTGGTAGTAGCA TTGGTGCTGGCAGTGGTAGTAGCATTAGTGCTGGAGTTGGTAGTCGCATTGGTAGTAGCACTAGTCCTGACGTTGGTGCTGGCAGTGGTAGTAGCATTAGTGCTGGAGTTGGTAGTCGCATTGGTACTGGCA CACTAGTCCTGACGTTGGTGCTGGCAGTGGTAGTAGCACTAGTCCTGACG TTAGTGTTGGAGTTGGTACTTTCAGTGGTAGTCGCACTAGTCCTGACGTT",
				"chrII 1955 2016 12 5.2 12 63 7 36 16 12 37 33 1.87 TTGGTACTGGCA TTGGTGCTGGCAGTGGTAGTAGCATTAGTGCTGGAGTTGGTAGTCGCATTGGTACTGGCATT CATTAGTGCTGGAGTTGGTAGTCGCATTGGTAGTAGCACTAGTCCTGACG AGTGTTGGAGTTGGTACTTTCAGTGGTAGTCGCACTAGTCCTGACGTTGA");
		String path = "tandemRepeats.txt";
		Path file = Paths.get(path);
		try {
			Files.write(file, lines, Charset.forName("UTF-8"));
			return path;
		} catch (IOException e) {
			e.printStackTrace();
			return null;	
		}
	}
}
