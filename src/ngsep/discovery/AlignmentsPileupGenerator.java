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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.GenomicVariant;


public class AlignmentsPileupGenerator {
	private Logger log = Logger.getLogger(AlignmentsPileupGenerator.class.getName());
	private List<PileupListener> listeners = new ArrayList<PileupListener>();
	private GenomicRegionSortedCollection<GenomicVariant> inputVariants = null;
	
	private List<GenomicVariant> seqInputVariants;
	private int idxNextInputVariant = 0;
	private QualifiedSequenceList sequencesMetadata;
	
	private String querySeq=null;
	private int queryFirst = 0;
	private int queryLast = 1000000000;
	private List<ReadAlignment> pendingAlignments = new ArrayList<ReadAlignment>();
	private List<ReadAlignment> sameStartPrimaryAlignments = new ArrayList<ReadAlignment>();
	private List<ReadAlignment> sameStartSecondaryAlignments = new ArrayList<ReadAlignment>();
	private QualifiedSequence currentReferenceSequence = null;
	private int currentReferencePos = 0;
	private int currentReferenceLast = 0;
	private boolean keepRunning = true;
	private int maxAlnsPerStartPos = 5;
	private boolean processSecondaryAlignments = false;
	private byte basesToIgnore5P = 0;
	private byte basesToIgnore3P = 0;
	private boolean ignoreXSField = false;
	

	public void addListener(PileupListener listener) {
		listeners.add(listener);
	}

	public void setInputVariants(GenomicRegionSortedCollection<GenomicVariant> inputVariants) {
		this.inputVariants = inputVariants;
		seqInputVariants = null;
	}
	
	public String getQuerySeq() {
		return querySeq;
	}


	public void setQuerySeq(String querySeq) {
		this.querySeq = querySeq;
	}


	public int getQueryFirst() {
		return queryFirst;
	}


	public void setQueryFirst(int queryFirst) {
		this.queryFirst = queryFirst;
	}


	public int getQueryLast() {
		return queryLast;
	}


	public void setQueryLast(int queryLast) {
		this.queryLast = queryLast;
	}

	public int getMaxAlnsPerStartPos() {
		return maxAlnsPerStartPos;
	}

	public void setMaxAlnsPerStartPos(int maxAlnsPerStartPos) {
		this.maxAlnsPerStartPos = maxAlnsPerStartPos;
	}

	
	public boolean isProcessSecondaryAlignments() {
		return processSecondaryAlignments;
	}

	public void setProcessSecondaryAlignments(boolean processSecondaryAlignments) {
		this.processSecondaryAlignments = processSecondaryAlignments;
	}
	
	public byte getBasesToIgnore5P() {
		return basesToIgnore5P;
	}

	public void setBasesToIgnore5P(byte basesToIgnore5P) {
		this.basesToIgnore5P = basesToIgnore5P;
	}

	public byte getBasesToIgnore3P() {
		return basesToIgnore3P;
	}

	public void setBasesToIgnore3P(byte basesToIgnore3P) {
		this.basesToIgnore3P = basesToIgnore3P;
	}
	
	public boolean isIgnoreXSField() {
		return ignoreXSField;
	}

	public void setIgnoreXSField(boolean ignoreXSField) {
		this.ignoreXSField = ignoreXSField;
	}

	public boolean isKeepRunning() {
		return keepRunning;
	}

	public void setKeepRunning(boolean keepRunning) {
		this.keepRunning = keepRunning;
	}

	public void processFile(String filename) throws IOException {
		ReadAlignmentFileReader reader = null;
		int processedAlns = 0;
		try {
			reader = new ReadAlignmentFileReader(filename);
			reader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_SEQUENCE);
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			if(!processSecondaryAlignments ) filterFlags+=ReadAlignment.FLAG_SECONDARY;
			reader.setFilterFlags(filterFlags);
			reader.setIgnoreXSField(ignoreXSField);
			sequencesMetadata = reader.getSequences();
			boolean querySeqFound = false;
			Iterator<ReadAlignment> it = reader.iterator();
			//Sequence under processing
			while(it.hasNext() && keepRunning) {
				ReadAlignment aln = it.next();
				if(inputVariants!=null && idxNextInputVariant>=inputVariants.size()) break;
				//System.out.println("Processing alignment at pos: "+alnRecord.getAlignmentStart()+". Seq: "+alnRecord.getReferenceName()+". Read name: "+alnRecord.getReadName());
				if(querySeq!=null) {
					if(querySeq.equals(aln.getSequenceName())) {
						querySeqFound = true;
						//Object reuse to make the equals method O(1) in future queries
						querySeq = aln.getSequenceName();
						if(aln.getFirst()>queryLast) break;
						if(queryFirst > aln.getLast()) continue;
					} else if(querySeqFound) {
						break;
					}
				}
				processAlignment(aln);
				processedAlns++;
				if(processedAlns%1000000 == 0) log.info("Processed "+processedAlns+" alignments");
			}
		} finally {
			if (reader!=null) reader.close();
		}
		
		if(keepRunning) notifyEndOfAlignments();
		else log.warning("Cancelled process");
	}
	
	public void processAlignment(ReadAlignment aln) {
		if(!spanInputVariants(aln)) return;
		
		if(currentReferenceSequence!=null) {
			if(!currentReferenceSequence.getName().equals(aln.getSequenceName())) {
				processPileups(currentReferenceLast+1);
				for(PileupListener listener:listeners) listener.onSequenceEnd(currentReferenceSequence);
				currentReferenceSequence=null;
			} else {
				processPileups(aln.getFirst());
			}	
		} 
		if (currentReferenceSequence==null) {
			startSequence(aln);
			for(PileupListener listener:listeners) listener.onSequenceStart(currentReferenceSequence);
		}
		aln.collapseIndelEvents();
		int alnLast = aln.getLast();
		if(alnLast > currentReferenceLast) currentReferenceLast = alnLast;
		if(aln.isSecondary()) {
			sameStartSecondaryAlignments.add(aln);
		} else {
			sameStartPrimaryAlignments.add(aln);
		}
	}

	private void startSequence(ReadAlignment aln) {
		String seqName = aln.getSequenceName();
		currentReferenceSequence = sequencesMetadata.get(seqName);
		if(currentReferenceSequence == null) {
			currentReferenceSequence = new QualifiedSequence(seqName);
		}
		log.info("Processing sequence "+aln.getSequenceName());
		currentReferencePos = aln.getFirst();
		currentReferenceLast = aln.getLast();
		if(inputVariants!=null) {
			seqInputVariants = inputVariants.getSequenceRegions(seqName).asList();
			idxNextInputVariant = 0;
		}
	}
	
	private boolean spanInputVariants(ReadAlignment aln) {
		if(inputVariants==null) return true;
		while(idxNextInputVariant<seqInputVariants.size())  {
			GenomicVariant variant = seqInputVariants.get(idxNextInputVariant);
			if(GenomicRegionSpanComparator.getInstance().span(variant, aln)) {
				return true;
			} else if (variant.getFirst()>aln.getLast()) break;
			idxNextInputVariant++;
		}
		return false;
	}

	
	public void notifyEndOfAlignments() {
		processPileups(currentReferenceLast+1);
		for(PileupListener listener:listeners) listener.onSequenceEnd(currentReferenceSequence);
	}
	private void processPileups(int alignmentStart) {
		if(alignmentStart==currentReferencePos) return;
		processSameReadPileups(alignmentStart);
		while(currentReferencePos<alignmentStart) {
			if(!processCurrentPosition()) {
				updatePendingAlns();
				if(pendingAlignments.size()==0) currentReferencePos = alignmentStart;
			}
		}
		updatePendingAlns();
	}

	private void updatePendingAlns() {
		List<ReadAlignment> alnsToKeep = new ArrayList<ReadAlignment>();
		for(ReadAlignment aln:pendingAlignments) {
			if(aln.getLast()>=currentReferencePos) alnsToKeep.add(aln);
		}
		pendingAlignments.clear();
		if(alnsToKeep.size()>0) pendingAlignments.addAll(alnsToKeep);
	}

	private void processSameReadPileups(int alignmentStart) {
		int posStart = 0;
		int pAlns = sameStartPrimaryAlignments.size();
		if(pAlns>0) {
			posStart = sameStartPrimaryAlignments.get(0).getFirst();
		} else if(sameStartSecondaryAlignments.size()>0) {
			posStart = sameStartSecondaryAlignments.get(0).getFirst();
		}
		if(posStart == 0 || posStart == alignmentStart) return;
		int remaining = maxAlnsPerStartPos;
		if(remaining<=0) remaining = Integer.MAX_VALUE; 
		for(int i=0;i<pAlns && remaining>0;i++) {
			ReadAlignment aln = sameStartPrimaryAlignments.get(i);
			aln.setBasesToIgnore5P(basesToIgnore5P);
			aln.setBasesToIgnore3P(basesToIgnore3P);
			pendingAlignments.add(aln);
			remaining--;
		}
		for(int i=0;i<sameStartSecondaryAlignments.size() && remaining>0;i++) {
			ReadAlignment aln = sameStartSecondaryAlignments.get(i);
			aln.setBasesToIgnore5P(basesToIgnore5P);
			aln.setBasesToIgnore3P(basesToIgnore3P);
			pendingAlignments.add(aln);
			remaining--;
		}
		sameStartPrimaryAlignments.clear();
		sameStartSecondaryAlignments.clear();
	}

	private boolean processCurrentPosition() {
		int posPrint = -1;
		if(pendingAlignments.size()==0) {
			currentReferencePos++;
			return false;
		}
		if(currentReferencePos==posPrint) System.out.println("Number of pending alignments: "+pendingAlignments.size()+". time: "+System.currentTimeMillis());
		
		PileupRecord pileup = new PileupRecord(currentReferenceSequence.getName(), currentReferencePos);
		for(ReadAlignment aln:pendingAlignments) {
			if(currentReferencePos==posPrint)System.out.println("Next pending: "+aln.getReadName()+" located at "+aln.getSequenceName()+":"+aln.getFirst()+"-"+aln.getLast()+". time: "+System.currentTimeMillis());
			pileup.addAlignment(aln);
		}
		if(currentReferencePos==posPrint)System.out.println("Number of alignments in pileup: "+pileup.getNumAlignments()+". time: "+System.currentTimeMillis());
		processPileup(pileup);
		if(currentReferencePos==posPrint)System.out.println("Processed pileup. time: "+System.currentTimeMillis());
		currentReferencePos++;
		return pileup.getNumAlignments()>0;
	}

	private void processPileup(PileupRecord pileup) {
		for(PileupListener listener:listeners) {
			listener.onPileup(pileup);
		}
	}

	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}
	
	//Gets the index with the maximum from the given values excluding the given index
	public static int getMaxIndex(int[] values,int excludeIndex) {
		int index = 0;
		if(excludeIndex == 0) {
			index=1;
		}
		for (int i=0;i<values.length;i++) {
			if(i!=excludeIndex && values[i]>values[index]) {
				index = i;
			}
		}
		return index;
	}

	
}
