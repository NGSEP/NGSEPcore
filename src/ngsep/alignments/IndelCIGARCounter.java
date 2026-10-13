package ngsep.alignments;

import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.variants.GenomicVariant;

public class IndelCIGARCounter {

	public static void main(String[] args) throws Exception {
		IndelCIGARCounter instance = new IndelCIGARCounter();
		int [][] indelsDist = instance.processFile(args[0]);
		for(int i=1;i<indelsDist.length;i++) {
			System.out.println(""+i+"\t"+indelsDist[i][0]+"\t"+indelsDist[i][1]);
		}
	}

	public  int[][] processFile(String filename) throws IOException {
		int[][] counts = new int [21][2];
		try (ReadAlignmentFileReader reader = new ReadAlignmentFileReader(filename)) {
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			filterFlags+=ReadAlignment.FLAG_SECONDARY;
			filterFlags+=ReadAlignment.FLAG_SUPPLEMENTARY;
			reader.setFilterFlags(filterFlags);
			
			Iterator<ReadAlignment> it = reader.iterator();
			while (it.hasNext()) {
				ReadAlignment aln = it.next();
				Map<Integer,GenomicVariant> indelCalls = aln.getIndelCalls();
				if(indelCalls==null) continue;
				for(GenomicVariant indel:indelCalls.values()) {
					int length = Math.min(counts.length-1, indel.length());
					if (indel.getLast()==indel.getFirst()+1) {
						(counts[length][0])++;
					} else {
						(counts[length][1])++;
					}
				}
			}
		}
		return counts;
	}

}
