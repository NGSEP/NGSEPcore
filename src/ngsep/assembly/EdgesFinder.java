package ngsep.assembly;

import java.util.List;
import java.util.Map;

public interface EdgesFinder {
	public Map<Integer, List<Edge>> getEdges();

	public Map<Integer, Edge> getEmbedded();

	public static void printStats(Map<Integer, List<Edge>> edges,
			Map<Integer, Edge> embed) {

	}
}

class Edge {
	private int indexSequence1;
	private int indexSequence2;
	private int first1;
	private int last1;
	private int first2;
	private int last2;
	private boolean negativeStrand;
	private boolean negativeCase;

	
	public Edge(int indexSequence1, int indexSequence2, int first1, int last1,
			int first2, int last2) {
		super();
		this.indexSequence1 = indexSequence1;
		this.indexSequence2 = indexSequence2;
		this.first1 = first1;
		this.last1 = last1;
		this.first2 = first2;
		this.last2 = last2;
	}

	public int getOverlap() {
		return last1 - first1;
	}

	public int getIndexSequence1() {
		return indexSequence1;
	}

	public void setIndexSequence1(int indexSequence1) {
		this.indexSequence1 = indexSequence1;
	}

	public int getIndexSequence2() {
		return indexSequence2;
	}

	public void setIndexSequence2(int indexSequence2) {
		this.indexSequence2 = indexSequence2;
	}

	public int getFirst1() {
		return first1;
	}

	public void setFirst1(int first1) {
		this.first1 = first1;
	}

	public int getLast1() {
		return last1;
	}

	public void setLast1(int last1) {
		this.last1 = last1;
	}

	public int getFirst2() {
		return first2;
	}

	public void setFirst2(int first2) {
		this.first2 = first2;
	}

	public int getLast2() {
		return last2;
	}

	public void setLast2(int last2) {
		this.last2 = last2;
	}

	public boolean isNegativeStrand() {
		return negativeStrand;
	}

	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}

	public boolean isNegativeCase() {
		return negativeCase;
	}

	public void setNegativeCase(boolean negativeCase) {
		this.negativeCase = negativeCase;
	}

	public int length() {
		return last2 - first2;
	}
}
