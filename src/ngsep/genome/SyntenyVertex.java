package ngsep.genome;


public class SyntenyVertex {
	private LocalHomologyCluster localRegion1;
	private LocalHomologyCluster localRegion2;
	public SyntenyVertex(LocalHomologyCluster localRegion1, LocalHomologyCluster localRegion2) {
		super();
		this.localRegion1 = localRegion1;
		this.localRegion2 = localRegion2;
	}
	public HomologyCluster getHomologyCluster() {
		return localRegion1.getParent();
	}
	public LocalHomologyCluster getLocalRegion1() {
		return localRegion1;
	}
	public LocalHomologyCluster getLocalRegion2() {
		return localRegion2;
	}
	
}
