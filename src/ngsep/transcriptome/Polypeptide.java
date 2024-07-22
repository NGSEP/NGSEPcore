package ngsep.transcriptome;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Polypeptide {
	private String id;
	private Transcript transcript;
	private List<String> products;
	private List<String> ontologyTerms;
	private Map<String,String> pfamTerms= new HashMap<>();
	public Polypeptide(String id, Transcript transcript) {
		super();
		this.id = id;
		this.transcript = transcript;
	}
	public List<String> getProducts() {
		return products;
	}
	public void setProducts(List<String> products) {
		this.products = products;
	}
	public String getId() {
		return id;
	}
	public Transcript getTranscript() {
		return transcript;
	}
	
	public List<String> getOntologyTerms() {
		return ontologyTerms;
	}
	public void setOntologyTerms(List<String> ontologyTerms) {
		this.ontologyTerms = ontologyTerms;
	}
	public Map<String, String> getPfamTerms() {
		return pfamTerms;
	}
	public void addPfamTerm(String id, String description) {
		pfamTerms.put(id, description);
	}
	
	
}
