package ngsep.sequences.test;


import java.util.List;
import java.util.Map;

import junit.framework.TestCase;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.MinimizersTable;
import ngsep.sequences.MinimizersTableEntry;

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
		Map<Integer, CharSequence> kmers = KmersExtractor.extractKmersAsMap(sequence, 15, 1, 0, sequence.length(), false, true, true);
		Map<Integer, List<MinimizersTableEntry>> minimizers = table.computeSequenceMinimizers(345, sequence.length(), kmers);
		Map<Integer, List<MinimizersTableEntry>> results = table.calculateMinimizerHits(345, minimizers);
		assertEquals(results.size(), 1);
		List<MinimizersTableEntry> resultsList = results.get(232);
		//assertEquals(minimizers.size(), resultsList.size());
		for (MinimizersTableEntry resultsEntry:resultsList) {
			int minResult = resultsEntry.getMinimizer();
			List<MinimizersTableEntry> originalEntries = minimizers.get(minResult);
			assertNotNull(originalEntries);
			System.out.println("Minimizer: "+minResult+" original entries: "+originalEntries.size());
			for(MinimizersTableEntry entry:originalEntries) {
				System.out.println("Minimizer: "+minResult+" checking start: "+entry.getStart()+" results entry start: "+resultsEntry.getStart());
				assertEquals(entry.getMinimizer(), resultsEntry.getMinimizer());
				assertEquals(entry.getSequenceId(),345);
				assertEquals(resultsEntry.getSequenceId(),232);
			}
		}
	}
}
