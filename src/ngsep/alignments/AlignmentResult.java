package ngsep.alignments;

import java.util.LinkedList;

class AlignmentResult {
	private int subjectStartIdx;
	private int subjectLastIdx;
	private int distance;
	private LinkedList<Character> path = new LinkedList<>();

	public void addBacktrack(char decision) {
		path.add(0, decision);
	}

	public String getCigarString () {
		StringBuilder cigar = new StringBuilder();
		int nextCount = 0;
		char next = 0;
		for(char c:path) {
			if(c!=next) {
				if(nextCount>0) {
					cigar.append(nextCount);
					cigar.append(next);
				}
				nextCount = 1;
				next = c;
			} else {
				nextCount++;
			}
		}
		if(nextCount>0) {
			cigar.append(nextCount);
			cigar.append(next);
		}
		return cigar.toString();
	}

	/**
	 * @return the subjectStartIdx
	 */
	public int getSubjectStartIdx() {
		return subjectStartIdx;
	}

	/**
	 * @param subjectStartIdx the subjectStartIdx to set
	 */
	public void setSubjectStartIdx(int subjectStartIdx) {
		this.subjectStartIdx = subjectStartIdx;
	}

	/**
	 * @return the subjectLastIdx
	 */
	public int getSubjectLastIdx() {
		return subjectLastIdx;
	}

	/**
	 * @param subjectLastIdx the subjectLastIdx to set
	 */
	public void setSubjectLastIdx(int subjectLastIdx) {
		this.subjectLastIdx = subjectLastIdx;
	}

	/**
	 * @return the distance
	 */
	public int getDistance() {
		return distance;
	}

	/**
	 * @param distance the distance to set
	 */
	public void setDistance(int distance) {
		this.distance = distance;
	}
}