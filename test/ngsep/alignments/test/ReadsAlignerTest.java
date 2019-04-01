package ngsep.alignments.test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import junit.framework.TestCase;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.ReadsAligner;
import ngsep.genome.GenomicRegion;

public class ReadsAlignerTest extends TestCase {
	private ReadsAligner readsAligner;

	public void setUpReadsAligner() {
		readsAligner= new ReadsAligner();
	}

	public void testCheckPairEnd() {
		setUpReadsAligner();

		//Useful flags tool https://broadinstitute.github.io/picard/explain-flags.html

		//Case 1: ProperPair
		int flags1=ReadAlignment.FLAG_PAIRED;
		ReadAlignment aln1=new ReadAlignment("sequenceName", 5000, 5249, 250, flags1);
		int flags2=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		ReadAlignment aln2=new ReadAlignment("sequenceName", 5500, 5251, 250, flags2);
		aln2.setMateNegativeStrand(true);

		boolean proper =readsAligner.pairAlignMents(aln1, aln2,true);
		assertEquals(true, proper);
		proper =readsAligner.pairAlignMents(aln2, aln1,true);
		assertEquals(true, proper);


		//Case 2: Long deletion
		flags1=ReadAlignment.FLAG_PAIRED;
		aln1=new ReadAlignment("sequenceName", 21000, 21249, 250, flags1);
		flags2=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		aln2=new ReadAlignment("sequenceName", 22000, 21751, 250, flags2);
		aln2.setMateNegativeStrand(true);


		boolean unProper =readsAligner.pairAlignMents(aln1, aln2,true);
		assertEquals(false, unProper);
		unProper =readsAligner.pairAlignMents(aln2, aln1,true);
		assertEquals(false, unProper);

		unProper =readsAligner.pairAlignMents(aln1, aln2,false);
		assertEquals(true, unProper);
		unProper =readsAligner.pairAlignMents(aln2, aln1,false);
		assertEquals(true, unProper);
	}

	public void testLoadTRF() {
		setUpReadsAligner();
		String path = setUpTRF();
		Map<String, List<GenomicRegion>> map = readsAligner.loadTRF(path);
		assertEquals(2, map.size());
		assertEquals(3, map.get("chrI").size());
		assertEquals(3, map.get("chrII").size());
		File f = new File(path);
		f.delete();
	}
	
	public void testCheckReadVsTRF() {
		setUpReadsAligner();
		String path = setUpTRF();
		Map<String, List<GenomicRegion>> map = readsAligner.loadTRF(path);
		
		//map.get(key)
		
				
		
		
		File f = new File(path);
		f.delete();
	}

	private String setUpTRF() {
		List<String> lines = Arrays.asList(
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
