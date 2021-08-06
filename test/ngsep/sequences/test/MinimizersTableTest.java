package ngsep.sequences.test;


import java.util.List;
import java.util.Map;

import junit.framework.TestCase;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.MinimizersTable;
import ngsep.sequences.KmerCodesTableEntry;

public class MinimizersTableTest extends TestCase {
	private String sequence = "CTCAACTAGATCGCACAACGTCGGAATGGTTTCATCCACAGATTGAATTTTTGGTTGCTGTATCAGTCCTTGAATGATGTCCATTCTTGATAGGAGGGTTGTTATAGATATTAATCACTCGAAGTCGTGAACAAGAAATTGTCTTCTCTCCAGTATTCAGTCTCTGTGAT";
	public void testEncodeDecode () {
		MinimizersTable table = new MinimizersTable(15, 5);
		List<KmerCodesTableEntry> minimizers = table.computeSequenceMinimizers(0,sequence,0, sequence.length());
		
		for(KmerCodesTableEntry entry:minimizers) {
			long code = entry.encode();
			KmerCodesTableEntry copy = new KmerCodesTableEntry(entry.getKmerCode(), code);
			assertEquals(entry.getKmerCode(), copy.getKmerCode());
			assertEquals(entry.getSequenceId(), copy.getSequenceId());
			assertEquals(entry.getStart(), copy.getStart());
		}
		
	}
	public void testMassiveEncodeDecode () {
		for (int i=0;i<100000;i++) {
			for(int j=0;j<100;j++) {
				KmerCodesTableEntry entry = new KmerCodesTableEntry(0,i,j);
				long code = entry.encode();
				KmerCodesTableEntry copy = new KmerCodesTableEntry(0, code);
				assertEquals(i, copy.getSequenceId());
				assertEquals(j, copy.getStart());
			}
		}
	}
}
