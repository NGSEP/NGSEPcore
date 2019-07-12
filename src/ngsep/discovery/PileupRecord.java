/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.discovery;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.alignments.ReadAlignment;
import ngsep.math.FisherExactTest;

/**
 * Class to store a pileup from many alignments spanning the same 
 * genomic location
 * @author Jorge Duitama
 */
public class PileupRecord {
	private String sequenceName;
	private int position;
	private List<ReadAlignment> alignmentsList = new ArrayList<>();
	private Map<String,List<ReadAlignment>> alignmentsMap = new HashMap<>();
	private int referenceSpan=1;
	private int numAlignments = 0;
	private int numUniqueAlns = 0;
	private int numNegativeStrandAlns = 0;
	private boolean str = false;
	private boolean newSTR = false;
	private boolean embedded = false;
	
	//DEBUG
	private int posPrint = -1;
	
	
	/**
	 * Creates a pileup record with the given information
	 * @param sequenceName Name of the reference sequence
	 * @param position Position in the reference sequence
	 * @param alleleCalls Allele calls for the given locus in pileup format
	 * @param qualityScores Quality scores of the allele calls
	 */
	public PileupRecord(String sequenceName, int position) {
		super();
		this.sequenceName = sequenceName;
		this.position = position;
	}
	
	public String getSequenceName() {
		return sequenceName;
	}
	/**
	 * Changes the reference sequence name
	 * @param sequenceName New sequence name
	 */
	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}
	
	/**
	 * @return int Position in the reference sequence
	 */
	public int getPosition() {
		return position;
	}
	/**
	 * Changes the position in the reference sequence
	 * @param position New position
	 */
	public void setPosition(int position) {
		this.position = position;
	}
	
	public int getReferenceSpan() {
		return referenceSpan;
	}

	public void setReferenceSpan(int referenceSpan) {
		this.referenceSpan = referenceSpan;
	}
	/**
	 * Calculates the allele calls from the pileup position with the given span
	 * @param referenceSpan Length in reference basepairs of the desired allele calls
	 * @param readGroups to return alignments. If null, allele calls for all alignments of this pileup are returned
	 * @return List<PileupAlleleCall> List of allele calls with the given reference span
	 */
	public List<PileupAlleleCall> getAlleleCalls(int referenceSpan, Set<String> readGroups) {
		List<PileupAlleleCall> alleleCalls = new ArrayList<>();
		if(readGroups==null) return getAlleleCalls(referenceSpan, (String)null);
		for(String readGroup:readGroups) {
			alleleCalls.addAll(getAlleleCalls(referenceSpan, readGroup));
		}
		return alleleCalls;
	}
	/**
	 * Calculates the allele calls from the pileup position with the given span considering all alignments spanning this pileup position
	 * @param referenceSpan Length in reference basepairs of the desired allele calls
	 * @return List<PileupAlleleCall> List of allele calls with the given reference span
	 */
	public List<PileupAlleleCall> getAlleleCalls(int referenceSpan) {
		return getAlleleCalls(referenceSpan, (String)null);
	}
	/**
	 * Calculates the allele calls from the pileup position with the given span
	 * @param referenceSpan Length in reference basepairs of the desired allele calls
	 * @param readGroup to return alignments. If null, allele calls for all alignments of this pileup are returned
	 * @return List<PileupAlleleCall> List of allele calls with the given refernce span
	 */
	public List<PileupAlleleCall> getAlleleCalls(int referenceSpan, String readGroup) {
		List<PileupAlleleCall> alleleCalls = new ArrayList<>();
		List<ReadAlignment> alignments;
		if(readGroup == null) alignments = getAlignments();
		else alignments = alignmentsMap.get(readGroup);
		if(alignments==null) return alleleCalls;
		for(ReadAlignment aln:alignments) { 
			CharSequence alleleCall = aln.getAlleleCall(position);
			if(position==posPrint) System.out.println("getAlleleCalls. Allele call: "+alleleCall+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+". Read name: "+aln.getReadName()+". CIGAR: "+aln.getCigarString()+" refSpan: "+referenceSpan+" negativeStrand: "+aln.isNegativeStrand()+". Ignore start: "+aln.getBasesToIgnoreStart()+" Ignore end: "+aln.getBasesToIgnoreEnd());
			if(alleleCall == null) continue;
			String alnQS = ""+aln.getBaseQualityScore(position);
			if(referenceSpan > 1) {
				int lastBase = position+referenceSpan-1;
				alleleCall = aln.getAlleleCall(position, lastBase);
				//if(position==posPrint) System.out.println("getAlleleCalls. With span: "+referenceSpan+". Allele call: "+alleleCall+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+". Read name: "+aln.getReadName()+". CIGAR: "+aln.getCigarString()+" negativeStrand: "+aln.isNegativeStrand()+". Ignore start: "+aln.getBasesToIgnoreStart()+" Ignore end: "+aln.getBasesToIgnoreEnd());
				if(alleleCall==null) continue;
				alnQS = aln.getBaseQualityScores(position, lastBase);
			} else if (alleleCall.length()>1) continue;
			if(position==posPrint) System.out.println("getAlleleCalls. With span: "+referenceSpan+". Allele call: "+alleleCall+". Quality score: "+alnQS+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+". Read name: "+aln.getReadName()+". CIGAR: "+aln.getCigarString()+" negativeStrand: "+aln.isNegativeStrand()+". Ignore start: "+aln.getBasesToIgnoreStart()+" Ignore end: "+aln.getBasesToIgnoreEnd()+" STR: "+str+" read group: "+aln.getReadGroup());
			PileupAlleleCall call = new PileupAlleleCall(alleleCall.toString(), alnQS);
			call.setReadGroup(aln.getReadGroup());
			call.setNegativeStrand(aln.isNegativeStrand());
			alleleCalls.add(call);
		}
		if(position==posPrint) System.out.println("getAlleleCalls. Final number of allele calls: "+alleleCalls.size());
		return alleleCalls;
	}

	public void addAlignment(ReadAlignment aln) {
		if(aln.getFirst()>position) return;
		if(aln.getLast()<position) return;
		alignmentsList.add(aln);
		List<ReadAlignment> alnsRG = alignmentsMap.get(aln.getReadGroup());
		if(alnsRG==null) {
			alnsRG = new ArrayList<>();
			alignmentsMap.put(aln.getReadGroup(), alnsRG);
		}
		alnsRG.add(aln);
		numAlignments++;
		if(aln.isUnique()) numUniqueAlns++;
		if(aln.isNegativeStrand()) numNegativeStrandAlns++;
	}
	
	public List<ReadAlignment> getAlignments() {
		return alignmentsList;
	}
	public int getNumAlignments() {
		return numAlignments;
	}

	public int getNumUniqueAlns() {
		return numUniqueAlns;
	}

	public boolean isSTR() {
		return str;
	}

	public void setSTR(boolean str) {
		this.str = str;
	}

	public boolean isNewSTR() {
		return newSTR;
	}

	public void setNewSTR(boolean newSTR) {
		this.newSTR = newSTR;
	}
	public boolean isInputSTR() {
		return this.str && !newSTR;
	}

	public boolean isEmbedded() {
		return embedded;
	}

	public void setEmbedded(boolean embedded) {
		this.embedded = embedded;
	}
	
	public double calculateStrandBiasPValue () {
		int a = numNegativeStrandAlns;
		int c = numAlignments-a;
		int b = numAlignments/2;
		int d = numAlignments-b;
		if(a>c) {
			//Set a to be the minimum for efficiency
			int t = a;
			a = c;
			c = t;
		}
		return FisherExactTest.calculatePValue(a, b, c, d);
	}
	
}
