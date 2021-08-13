package ngsep.genome;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;


/**
 * @author Laura Gonzalez
 */
public class HomologyCluster 
{

	// Constants for default values
	public final static String SCORE = "soft_core";
	public final static String SACCESORY = "soft_accesory";
	public final static String ECORE = "exact_core";
	public final static String EACCESORY = "exact_accesory";
	public final static String PERSISTENT = "persistent";
	public final static String SHELL = "shell";
	public final static String CLOUD = "cloud";
	
	// Attributes
	private List<HomologyUnit> homologyUnitsCluster = new ArrayList<>();
	private Logger log;
	private int clusterId;
	private double frequency;
	private String exactCat;
	private String softCat;
	
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}
	
	public HomologyCluster(int clusterId, List<HomologyUnit> homologyUnitsCluster) {
		this.clusterId = clusterId;
		this.homologyUnitsCluster = homologyUnitsCluster;
	}
	
	public List<HomologyUnit> getHomologyUnitsCluster() {
		return this.homologyUnitsCluster;
	}
	
	public int getClusterId() {
		return this.clusterId;
	}
	
	public double getFrequency() {
		return this.frequency;
	}
	
	public void setFrequency(double frequency) {
		this.frequency=frequency;
	}
	
	public String getExactCategory() {
		return this.exactCat;
	}
	
	public void setExactCategory(String exactCat) {
		this.exactCat = exactCat;
	}
	
	public String getSoftCategory() {
		return this.softCat;
	}
	
	public void setSoftCategory(String softCat) {
		this.softCat = softCat;
	}
}
