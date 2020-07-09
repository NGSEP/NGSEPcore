/**
 * Copied from SingleIndividualHaplotyper - Efficient heuristic algorithms for the SIH problem
 * Based on the initial implementation of the Refhap algorithm for haplotyping
 * @author Jorge Duitama, Christian Chavarro Espejo, Daniel Bautista
 */

package ngsep.haplotyping;


public class Refhap3SIHAlgorithm  implements SIHAlgorithm 
{
	private boolean [] cut;
	private byte [] haplotype;
	
	/**
	 * Build the haplotypes using a greedy version of the original RefHap Algorithm
	 */
	public void buildHaplotype(HaplotypeBlock block) 
	{
		FragmentsCutBuilder builder = new FragmentsCutBuilder(block);
		builder.calculateMaxCutStrategy3();
		cut = builder.getCut();
		haplotype=CutHaplotypeTranslator.getHaplotype(block, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);
		block.setHaplotype(haplotype);
		
	}

}
