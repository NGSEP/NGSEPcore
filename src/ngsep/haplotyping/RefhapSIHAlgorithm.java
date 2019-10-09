package ngsep.haplotyping;

import java.util.logging.Logger;

/**
 * Copied from SingleIndividualHaplotyper - Efficient heuristic algorithms for the SIH problem
 * @author Jorge Duitama
 */
public class RefhapSIHAlgorithm implements SIHAlgorithm {
	private static Logger log = Logger.getAnonymousLogger();
	private boolean [] cut;
	private byte [] haplotype;
	
	/**
	 * @return the log
	 */
	public static Logger getLog() {
		return log;
	}

	/**
	 * @param log the log to set
	 */
	public static void setLog(Logger log) {
		RefhapSIHAlgorithm.log = log;
	}



	@Override
	public void buildHaplotype(HaplotypeBlock block) {
		System.err.println("Phasing block");
		FragmentsCutBuilder builder = new FragmentsCutBuilder(block);
		System.err.println("Built graph");
		builder.calculateMaxCut();
		cut = builder.getCut();
		System.err.println("Calculated cut");
		haplotype=CutHaplotypeTranslator.getHaplotype(block, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);
		System.err.println("Assembled haplotypes");
		block.setHaplotype(haplotype);
	}	
}
