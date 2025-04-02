package ngsep.transposons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.math.PhredScoreHelper;

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
	
	public boolean isLTR() {
		return "LTR".equals(order);
	}
	public String toString() {
		return order+"/"+id;
	}
	
	public int calculateMatchScoreDP(List<TransposonDomainAlignment> alns) {
		int n1 = teDomainIds.size();
		int n2 = alns.size();
		int [][] dp = new int [n1+1][n2+1];
		for(int i=0;i<=n1;i++) {
			for(int j=0;j<=n2;j++) {
				if(i==0 || j==0) {
					dp[i][j] = 0;
					continue;
				}
				String d1 = teDomainIds.get(i-1);
				TransposonDomainAlignment aln = alns.get(j-1);
				int score = PhredScoreHelper.calculatePhredScore(aln.getEvalue());
				if(d1.equals(aln.getDomainCode())) {
					dp[i][j] = dp[i-1][j-1]+score;
				} else {
					dp[i][j] = Math.max(dp[i-1][j]-10, dp[i][j-1]-score);
				}
			}
		}
		return dp[n1][n2];
	}
	
	
	public static TransposableElementFamily matchFamily(List<TransposonDomainAlignment> alns) {
		return matchFamilyDP(alns);
	}
	
	public static TransposableElementFamily matchFamilyDP(List<TransposonDomainAlignment> alns) {
		List<TransposableElementFamily> allFams = new ArrayList<TransposableElementFamily>();
		for(List<TransposableElementFamily> famList : TE_FAMILIES.values()) {
			for(TransposableElementFamily family:famList) {
				allFams.add(family);
			}
		}
		int n = allFams.size();
		int [] scores = new int[n];
		for(int i=0;i<n;i++) {
			TransposableElementFamily family = allFams.get(i);
			scores[i] = family.calculateMatchScoreDP(alns);
			//System.out.println("Score family: "+family+":"+scores[i]);
		}
		double max = 0;
		double secondMax = 0;
		TransposableElementFamily answer = null;
		TransposableElementFamily answer2 = null;
		for(int i=0;i<n;i++) {
			if(scores[i]>=max) {
				secondMax = max;
				max = scores[i];
				answer2 = answer;
				answer = allFams.get(i);
			}
		}
		//System.out.println("Max: "+max+" second: "+secondMax+" order1: "+answer.getOrder()+" order2: "+answer2.getOrder());
		if(max-secondMax<10) {
			if (answer.getOrder().equals(answer2.getOrder())) return findUnknown(answer.getOrder());
			return null;
		}
		return answer;
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
		
		String [] rirDomains = {"RT","END"};
		ineFamilies.add(new TransposableElementFamily("LINE","RIR", Arrays.asList(rirDomains)));
		
		String [] riiDomains = {"END", "RT", "RH"};
		ineFamilies.add(new TransposableElementFamily("LINE","RII", Arrays.asList(riiDomains)));
		
		String [] rilDomains = {"END", "RT"};
		ineFamilies.add(new TransposableElementFamily("LINE","RIL", Arrays.asList(rilDomains)));
		
		String [] ritDomains = {"END", "RT"};
		ineFamilies.add(new TransposableElementFamily("LINE","RIT", Arrays.asList(ritDomains)));
		
		String [] dtpDomains = {"HTH"};
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
