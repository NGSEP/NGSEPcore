package ngsep.main;

import java.util.List;

import ngsep.discovery.CountsHelper;

public class CacheRegistry {
	public static void initializeCacheResources() {
		CountsHelper helper= new CountsHelper();
		List<String> alleles = helper.getAllelesList();
		//System.err.println("Initialized counts helper. Alleles: "+alleles);
	}
}
