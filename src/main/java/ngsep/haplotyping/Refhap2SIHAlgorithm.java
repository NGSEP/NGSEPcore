package ngsep.haplotyping;

/**
 * Copied from SingleIndividualHaplotyper - Efficient heuristic algorithms for the SIH problem
 * Based on the initial implementation of the Refhap algorithm for haplotyping
 * @author Jorge Duitama, Christian Chavarro Espejo, Daniel Bautista
 */
public class Refhap2SIHAlgorithm implements SIHAlgorithm 
{
	private boolean [] cut;
	private byte [] haplotype;
	
	/**
	 * Build the haplotypes using a greedy version of the original RefHap Algorithm
	 */
	@Override
	public void buildHaplotype(HaplotypeBlock block) 
	{
		FragmentsCutBuilder builder = new FragmentsCutBuilder(block);
		builder.calculateMaxCutStrategy2();
		cut = builder.getCut();
		haplotype=CutHaplotypeTranslator.getHaplotype(block, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);
		block.setHaplotype(haplotype);
	}	
}
