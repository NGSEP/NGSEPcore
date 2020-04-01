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
		Map<Integer, List<MinimizersTableEntry>> minimizers = table.computeSequenceMinimizers(0, sequence.length(), kmers);
		for (int minimizer:minimizers.keySet()) {
			for(MinimizersTableEntry entry:minimizers.get(minimizer)) {
				long code = entry.encode();
				MinimizersTableEntry copy = new MinimizersTableEntry(minimizer, code);
				assertEquals(entry.getMinimizer(), copy.getMinimizer());
				assertEquals(entry.getSequenceId(), copy.getSequenceId());
				assertEquals(entry.getStart(), copy.getStart());
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
