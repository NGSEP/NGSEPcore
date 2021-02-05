
package ngsep.haplotyping;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Copied from SingleIndividualHaplotyper - Efficient heuristic algorithms for the SIH problem
 * @author Jorge Duitama
 */
public interface SIHAlgorithm {
	/**
	 * Builds a haplotype for the given matrix representing a block of overlapping fragments
	 * @param block that represents the matrix of overlapping fragments 
	 */
	public void buildHaplotype (HaplotypeBlock block);
	public Logger getLog();
	public void setLog(Logger log);
	public static List<List<HaplotypeFragment>> buildClusters(HaplotypeBlock block, boolean[] cut) {
		int n = cut.length;
		List<HaplotypeFragment> fragments = block.getFragments();
		List<List<HaplotypeFragment>> answer = new ArrayList<List<HaplotypeFragment>>(2);
		List<HaplotypeFragment> cluster0 = new ArrayList<HaplotypeFragment>();
		answer.add(cluster0);
		List<HaplotypeFragment> cluster1 = new ArrayList<HaplotypeFragment>();
		answer.add(cluster1);
		for(int i=0;i<n;i++) {
			HaplotypeFragment fragment = fragments.get(i);
			if(cut[i]) cluster1.add(fragment);
			else cluster0.add(fragment);
		}
		return answer;
	}
}
