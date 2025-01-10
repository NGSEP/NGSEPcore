package ngsep.transposons;

import java.util.List;

public class TransposableElementFamily {
	private String id;
	private List<String> teDomainIds;
	public TransposableElementFamily(String id, List<String> teDomainIds) {
		super();
		this.id = id;
		this.teDomainIds = teDomainIds;
	}
	
	public String getId() {
		return id;
	}

	public boolean matchFamily(List<TransposonDomainAlignment> alns) {
		if(teDomainIds.size()!=alns.size()) return false;
		for(int i=0;i<alns.size();i++) {
			String domainId = teDomainIds.get(i);
			String alnDId = alns.get(i).getAlnDomain().getHmmID();
			if(!domainId.equals(alnDId)) return false;
		}
		return true;
	}
}
