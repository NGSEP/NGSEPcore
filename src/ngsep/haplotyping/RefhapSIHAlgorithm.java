package ngsep.haplotyping;

import java.util.logging.Logger;

/**
 * Copied from SingleIndividualHaplotyper - Efficient heuristic algorithms for the SIH problem
 * @author Jorge Duitama
 */
public class RefhapSIHAlgorithm implements SIHAlgorithm {
	private Logger log = Logger.getAnonymousLogger();
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}

	@Override
	public void buildHaplotype(HaplotypeBlock block) {
		FragmentsCutBuilder builder = new FragmentsCutBuilder(block);
		log.info("Built graph");
		builder.calculateMaxCut();
		log.info("Calculated cut");
		boolean [] cut = builder.getCut();
		byte [] haplotype=CutHaplotypeTranslator.getHaplotype(block, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);
		log.info("Calculated haplotypes");
		block.setHaplotype(haplotype);
		block.setFragmentsClusters(SIHAlgorithm.buildClusters(block, cut));
	}	
}
