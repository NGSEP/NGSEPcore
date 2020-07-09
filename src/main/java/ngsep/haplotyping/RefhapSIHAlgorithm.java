package ngsep.haplotyping;

/**
 * Copied from SingleIndividualHaplotyper - Efficient heuristic algorithms for the SIH problem
 * @author Jorge Duitama
 */
public class RefhapSIHAlgorithm implements SIHAlgorithm {
	private boolean [] cut;
	private byte [] haplotype;
	
	@Override
	public void buildHaplotype(HaplotypeBlock block) {
		FragmentsCutBuilder builder = new FragmentsCutBuilder(block);
		builder.calculateMaxCut();
		cut = builder.getCut();
		haplotype=CutHaplotypeTranslator.getHaplotype(block, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);
		block.setHaplotype(haplotype);
	}	
}
