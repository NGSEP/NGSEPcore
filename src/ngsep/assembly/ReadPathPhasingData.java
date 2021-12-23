package ngsep.assembly;

public class ReadPathPhasingData {
	private int readId;
	private int pathId;
	private int blockNumber;
	private int phaseWithinBlock = -1;
	private double readDepth=0;
	private boolean inHomozygousRegion = false;
	public ReadPathPhasingData(int readId, int pathId, int blockNumber) {
		super();
		this.readId = readId;
		this.pathId = pathId;
		this.blockNumber = blockNumber;
	}
	public int getPhaseWithinBlock() {
		return phaseWithinBlock;
	}
	public void setPhaseWithinBlock(int phaseWithinBlock) {
		this.phaseWithinBlock = phaseWithinBlock;
	}
	public double getReadDepth() {
		return readDepth;
	}
	public void setReadDepth(double readDepth) {
		this.readDepth = readDepth;
	}
	public int getReadId() {
		return readId;
	}
	public int getPathId() {
		return pathId;
	}
	public int getBlockNumber() {
		return blockNumber;
	}
	public boolean isPhased () {
		return phaseWithinBlock>=0;
	}
	
	public boolean isInHomozygousRegion() {
		return inHomozygousRegion;
	}
	public void setInHomozygousRegion(boolean inHomozygousRegion) {
		this.inHomozygousRegion = inHomozygousRegion;
		if(inHomozygousRegion) phaseWithinBlock = -1;
	}
	
	
	public boolean isOppositePhase(ReadPathPhasingData read2) {
		if(this.pathId != read2.pathId) return false;
		if(this.blockNumber != read2.blockNumber) return false;
		if(this.phaseWithinBlock == read2.phaseWithinBlock) return false;
		return true;
	}
	
}
