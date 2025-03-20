package ngsep.transposons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class TransposableElementFamily {
	private String id;
	private String order;
	private List<String> teDomainIds;
	private static final Map<String,List<TransposableElementFamily>> TE_FAMILIES = loadTEFamilies();
	
	public static final String UNKNOWN = "Unknown";
	public static final TransposableElementFamily LTR_UNKNOWN = new TransposableElementFamily("LTR", UNKNOWN, new ArrayList<String>());
	public static final TransposableElementFamily LINE_UNKNOWN = new TransposableElementFamily("LINE", UNKNOWN, new ArrayList<String>());
	public static final TransposableElementFamily SINE_UNKNOWN = new TransposableElementFamily("SINE", UNKNOWN, new ArrayList<String>());
	private TransposableElementFamily(String order, String id, List<String> teDomainIds) {
		super();
		this.order = order;
		this.id = id;
		this.teDomainIds = teDomainIds;
	}
	
	
	public String getOrder() {
		return order;
	}

	public String getId() {
		return id;
	}
	
	private boolean matchStrict (List<String> alignedDomainIds) {
		if(teDomainIds.size()!=alignedDomainIds.size()) return false;
		for(int i=0;i<teDomainIds.size();i++) {
			String domainId = teDomainIds.get(i);
			String alnDId = alignedDomainIds.get(i);
			if(!domainId.equals(alnDId)) return false;
		}
		return true;
	}

	public static TransposableElementFamily matchFamily(List<TransposonDomainAlignment> alns) {
		List<String> alignedDomainIds = new ArrayList<String>();
		int gagPos = -1;
		double gagEvalue = 10;
		int intPos = -1;
		int rtPos = -1;
		int rhPos = -1;
		int envPos = -1;
		int  yrPos = -1;
		boolean endPresent = false;
		boolean hthPresent = false;
		double hthEvalue = 10;
		int pos=0;
		for(TransposonDomainAlignment aln:alns) {
			String domainCode = aln.getDomainCode();
			alignedDomainIds.add(domainCode);
			if("GAG".equals(domainCode)) {
				gagPos = pos;
				gagEvalue = Math.min(gagEvalue,aln.getEvalue());
			}
			if("INT".equals(domainCode)) intPos = pos;
			if("RT".equals(domainCode)) rtPos = pos;
			if("RNASEH".equals(domainCode)) rhPos = pos;
			if("ENV".equals(domainCode)) envPos = pos;
			if("YR".equals(domainCode)) yrPos = pos;
			if("END".equals(domainCode)) endPresent = true;
			if("HTH".equals(domainCode)) {
				hthPresent = true;
				hthEvalue = Math.min(hthEvalue, aln.getEvalue());
			}
			pos++;
		}
		int nAlnDoms = alignedDomainIds.size();
		if(gagPos>=0 && hthPresent) {
			if(gagEvalue< hthEvalue) hthPresent = false;
			else gagPos = -1;
		}
		if(hthPresent) return findFamily(TE_FAMILIES.get("TIR"),"TIR");
		if(gagPos>=0) {
			if(intPos>gagPos) {
				//LTR
				List<TransposableElementFamily> fams = TE_FAMILIES.get("LTR1"); 
				for(TransposableElementFamily family:fams) {
					if(family.matchStrict(alignedDomainIds)) return family;
				}
				//TODO: Return something more specific in these cases
				if(envPos>intPos) return findFamily(fams,"RLR");
				if(envPos>=0 || nAlnDoms==2) return LTR_UNKNOWN;
				if(intPos<nAlnDoms-1) {
					if(rtPos>intPos || rhPos>intPos) return findFamily(fams,"RLC");
					return LTR_UNKNOWN;
				}
				return findFamily(fams,"RLG");
				
			} else if (yrPos>gagPos) {
				//DIRS
				return findFamily(TE_FAMILIES.get("YR"),"RYD");
			} else {
				return LTR_UNKNOWN;
			}
		}
		if(intPos>=0) {
			List<TransposableElementFamily> fams = TE_FAMILIES.get("LTR1");
			//TODO: Return something more specific in these cases
			if(alignedDomainIds.contains("ENV")) return findFamily(fams,"RLR");
			int posInt = alignedDomainIds.indexOf("INT");
			if(posInt<alignedDomainIds.size()-1) return findFamily(fams,"RLC");
			return LTR_UNKNOWN;
		}
		if(yrPos>=0) {
			if(alignedDomainIds.size()==1) return findFamily(TE_FAMILIES.get("YR"),"DYC");
			else return findFamily(TE_FAMILIES.get("YR"),"RYD");
		}
		//LINEs
		if(endPresent) {
			List<TransposableElementFamily> fams = TE_FAMILIES.get("INE"); 
			for(TransposableElementFamily family:fams) {
				if(family.matchStrict(alignedDomainIds)) return family;
			}
			if(endPresent) {
				if(rhPos>=0) return findFamily(fams,"RII");
				return findFamily(fams,"RIT");
			}
			return findFamily(fams,"RIR");
		}
		
		return null;
	}
	
	private static TransposableElementFamily findFamily(List<TransposableElementFamily> fams, String id) {
		for(TransposableElementFamily family:fams) {
			if(id.equals(family.getId())) return family;
		}
		return null;
	}


	private static Map<String,List<TransposableElementFamily>> loadTEFamilies() {
		List<TransposableElementFamily> ltr1Families = new ArrayList<TransposableElementFamily>();
		List<TransposableElementFamily> yrFamilies = new ArrayList<TransposableElementFamily>();
		List<TransposableElementFamily> ineFamilies = new ArrayList<TransposableElementFamily>();
		List<TransposableElementFamily> tirFamilies = new ArrayList<TransposableElementFamily>();
		
		String [] rlrDomains = {"GAG","AP","RT","RNASEH","INT","ENV"};
		ltr1Families.add(new TransposableElementFamily("LTR","RLR", Arrays.asList(rlrDomains)));
		
		String [] rlgDomains = {"GAG","AP","RT","RNASEH","INT"};
		ltr1Families.add(new TransposableElementFamily("LTR","RLG", Arrays.asList(rlgDomains)));
		
		String [] rlcDomains = {"GAG","AP","INT", "RT","RNASEH"};
		ltr1Families.add(new TransposableElementFamily("LTR","RLC", Arrays.asList(rlcDomains)));
		
		String [] rydDomains = {"GAG","AP","RT","RNASEH","YR"};
		yrFamilies.add(new TransposableElementFamily("DIRS","RYD", Arrays.asList(rydDomains)));
		
		String [] cryptonDomains = {"YR"};
		yrFamilies.add(new TransposableElementFamily("Crypton","DYC", Arrays.asList(cryptonDomains)));
		
		String [] rirDomains = {"RT","EN"};
		ineFamilies.add(new TransposableElementFamily("LINE","RIR", Arrays.asList(rirDomains)));
		
		String [] riiDomains = {"APE", "RT", "RH"};
		ineFamilies.add(new TransposableElementFamily("LINE","RII", Arrays.asList(riiDomains)));
		
		String [] ritDomains = {"APE", "RT"};
		ineFamilies.add(new TransposableElementFamily("LINE","RIT", Arrays.asList(ritDomains)));
		
		String [] dtpDomains = {"Tase"};
		tirFamilies.add(new TransposableElementFamily("TIR","TIR", Arrays.asList(dtpDomains)));
		
		Map<String,List<TransposableElementFamily>> answer = new HashMap<String, List<TransposableElementFamily>>();
		answer.put("LTR1",ltr1Families);
		answer.put("YR",yrFamilies);
		answer.put("INE",ineFamilies);
		answer.put("TIR",tirFamilies);
		return answer;
	}


	public static TransposableElementFamily findFamily(String orderStr, String familyStr) {
		for(List<TransposableElementFamily> famList : TE_FAMILIES.values()) {
			for(TransposableElementFamily family:famList) {
				if(family.getOrder().equals(orderStr) && family.getId().equals(familyStr)) return family;
			}
		}
		return null;
	}


	public static TransposableElementFamily findUnknown(String orderStr) {
		if("LTR".equals(orderStr)) return LTR_UNKNOWN;
		if("LINE".equals(orderStr)) return LINE_UNKNOWN;
		if("SINE".equals(orderStr)) return SINE_UNKNOWN;
		return null;
	}
}
