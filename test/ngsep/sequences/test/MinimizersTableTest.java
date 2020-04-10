package ngsep.sequences.test;


import java.util.List;
import java.util.Map;

import junit.framework.TestCase;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.MinimizersTable;
import ngsep.sequences.MinimizersTableEntry;
import ngsep.sequences.UngappedSearchHit;

public class MinimizersTableTest extends TestCase {
	private String sequence = "CTCAACTAGATCGCACAACGTCGGAATGGTTTCATCCACAGATTGAATTTTTGGTTGCTGTATCAGTCCTTGAATGATGTCCATTCTTGATAGGAGGGTTGTTATAGATATTAATCACTCGAAGTCGTGAACAAGAAATTGTCTTCTCTCCAGTATTCAGTCTCTGTGAT";
	public void testEncodeDecode () {
		MinimizersTable table = new MinimizersTable(15, 5);
		Map<Integer, CharSequence> kmers = KmersExtractor.extractKmersAsMap(sequence, 15, 1, 0, sequence.length(), false, true, true);
		List<MinimizersTableEntry> minimizers = table.computeSequenceMinimizers(0, 0, sequence.length(), kmers);
		
		for(MinimizersTableEntry entry:minimizers) {
			long code = entry.encode();
			MinimizersTableEntry copy = new MinimizersTableEntry(entry.getMinimizer(), code);
			assertEquals(entry.getMinimizer(), copy.getMinimizer());
			assertEquals(entry.getSequenceId(), copy.getSequenceId());
			assertEquals(entry.getStart(), copy.getStart());
		}
		
	}
	public void testMassiveEncodeDecode () {
		for (int i=0;i<100000;i++) {
			for(int j=0;j<100000;j++) {
				MinimizersTableEntry entry = new MinimizersTableEntry(0,i,j);
				long code = entry.encode();
				MinimizersTableEntry copy = new MinimizersTableEntry(0, code);
				assertEquals(i, copy.getSequenceId());
				assertEquals(j, copy.getStart());
			}
		}
	}
	public void testSearchSelf () {
		MinimizersTable table = new MinimizersTable(15, 5);
		table.addSequence(232, sequence);
		Map<Integer,List<UngappedSearchHit>> hits = table.match(sequence);
		assertEquals(hits.size(), 1);
		List<UngappedSearchHit> resultsList = hits.get(232);
		//assertEquals(minimizers.size(), resultsList.size());
	}
}
