package ngsep.transcriptome;

public class VariantAnnotationParameters {
	public static final int DEF_UPSTREAM=1000;
	public static final int DEF_DOWNSTREAM=300;
	public static final int DEF_SPLICE_DONOR=2;
	public static final int DEF_SPLICE_ACCEPTOR=2;
	public static final int DEF_SPLICE_REGION_INTRON=10;
	public static final int DEF_SPLICE_REGION_EXON=2;
	
	private int offsetUpstream = DEF_UPSTREAM;
	private int offsetDownstream = DEF_DOWNSTREAM;
	private int spliceDonorOffset = DEF_SPLICE_DONOR;
	private int spliceAcceptorOffset = DEF_SPLICE_ACCEPTOR;
	private int spliceRegionIntronOffset = DEF_SPLICE_REGION_INTRON;
	private int spliceRegionExonOffset = DEF_SPLICE_REGION_EXON;
	/**
	 * @return the offsetUpstream
	 */
	public int getOffsetUpstream() {
		return offsetUpstream;
	}
	/**
	 * @param offsetUpstream the offsetUpstream to set
	 */
	public void setOffsetUpstream(int offsetUpstream) {
		this.offsetUpstream = offsetUpstream;
	}
	/**
	 * @return the offsetDownstream
	 */
	public int getOffsetDownstream() {
		return offsetDownstream;
	}
	/**
	 * @param offsetDownstream the offsetDownstream to set
	 */
	public void setOffsetDownstream(int offsetDownstream) {
		this.offsetDownstream = offsetDownstream;
	}
	/**
	 * @return the spliceDonorOffset
	 */
	public int getSpliceDonorOffset() {
		return spliceDonorOffset;
	}
	/**
	 * @param spliceDonorOffset the spliceDonorOffset to set
	 */
	public void setSpliceDonorOffset(int spliceDonorOffset) {
		this.spliceDonorOffset = spliceDonorOffset;
	}
	/**
	 * @return the spliceAcceptorOffset
	 */
	public int getSpliceAcceptorOffset() {
		return spliceAcceptorOffset;
	}
	/**
	 * @param spliceAcceptorOffset the spliceAcceptorOffset to set
	 */
	public void setSpliceAcceptorOffset(int spliceAcceptorOffset) {
		this.spliceAcceptorOffset = spliceAcceptorOffset;
	}
	/**
	 * @return the spliceRegionIntronOffset
	 */
	public int getSpliceRegionIntronOffset() {
		return spliceRegionIntronOffset;
	}
	/**
	 * @param spliceRegionIntronOffset the spliceRegionIntronOffset to set
	 */
	public void setSpliceRegionIntronOffset(int spliceRegionIntronOffset) {
		this.spliceRegionIntronOffset = spliceRegionIntronOffset;
	}
	/**
	 * @return the spliceRegionExonOffset
	 */
	public int getSpliceRegionExonOffset() {
		return spliceRegionExonOffset;
	}
	/**
	 * @param spliceRegionExonOffset the spliceRegionExonOffset to set
	 */
	public void setSpliceRegionExonOffset(int spliceRegionExonOffset) {
		this.spliceRegionExonOffset = spliceRegionExonOffset;
	}
}
