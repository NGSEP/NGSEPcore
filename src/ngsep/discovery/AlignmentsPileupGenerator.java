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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegionComparator;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;


public class AlignmentsPileupGenerator {
	
	public static final int DEF_MAX_ALNS_PER_START_POS = 5;
	
	private Logger log = Logger.getLogger(AlignmentsPileupGenerator.class.getName());
	private List<PileupListener> listeners = new ArrayList<PileupListener>();
	private QualifiedSequenceList sequencesMetadata;
	
	private String querySeq=null;
	private int queryFirst = 0;
	private int queryLast = 1000000000;
	
	private boolean keepRunning = true;
	private int maxAlnsPerStartPos = DEF_MAX_ALNS_PER_START_POS;
	private boolean processNonUniquePrimaryAlignments = false;
	private boolean processSecondaryAlignments = false;
	private byte basesToIgnore5P = 0;
	private byte basesToIgnore3P = 0;
	private int minMQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	// Internal attributes to follow up the pileup process
	private QualifiedSequence currentReferenceSequence = null;
	private int currentReferencePos = 0;
	private int currentReferenceLast = 0;
	private List<ReadAlignment> pendingAlignments = new ArrayList<ReadAlignment>();
	
	private List<ReadAlignment> sameStartPrimaryAlignments = new ArrayList<ReadAlignment>();
	private List<ReadAlignment> sameStartSecondaryAlignments = new ArrayList<ReadAlignment>();
	private int lastReadAlignmentStart = 0;
	

	public void addListener(PileupListener listener) {
		listeners.add(listener);
	}
	
	/**
	 * @return the sequencesMetadata
	 */
	public QualifiedSequenceList getSequencesMetadata() {
		return sequencesMetadata;
	}

	/**
	 * @param sequencesMetadata the sequencesMetadata to set
	 */
	public void setSequencesMetadata(QualifiedSequenceList sequencesMetadata) {
		this.sequencesMetadata = sequencesMetadata;
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

	
	/**
	 * @return the processNonUniquePrimaryAlignments
	 */
	public boolean isProcessNonUniquePrimaryAlignments() {
		return processNonUniquePrimaryAlignments;
	}

	/**
	 * @param processNonUniquePrimaryAlignments the processNonUniquePrimaryAlignments to set
	 */
	public void setProcessNonUniquePrimaryAlignments(boolean processNonUniquePrimaryAlignments) {
		this.processNonUniquePrimaryAlignments = processNonUniquePrimaryAlignments;
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
	
	/**
	 * @return the minMQ
	 */
	public int getMinMQ() {
		return minMQ;
	}
	/**
	 * @param minMQ the minMQ to set
	 */
	public void setMinMQ(int minMQ) {
		this.minMQ = minMQ;
	}

	public boolean isKeepRunning() {
		return keepRunning;
	}

	public void setKeepRunning(boolean keepRunning) {
		this.keepRunning = keepRunning;
	}
	
	/**
	 * Parallel processing of several bam files
	 * PRE: The generator has the sequences metadata information
	 * @param alignmentFiles List of files to process in parallel
	 * @throws IOException If there is an I/O error processing the files
	 */
	public void processFiles(List<String> alignmentFiles) throws IOException {
		processFiles(alignmentFiles, new ArrayList<>());
	}
	/**
	 * Parallel processing of several bam files
	 * PRE: The generator has the sequences metadata information
	 * @param alignmentFiles List of files to process in parallel
	 * @param replaceReadGroups Read groups to use instead of the read groups within each file.
	 * Can not be null. If empty, default read groups are used
	 * @throws IOException If there is an I/O error processing the files
	 */
	public void processFiles(List<String> alignmentFiles, List<String> replaceReadGroups) throws IOException {
		int n = alignmentFiles.size();
		if(n==0) return;
		if(n==1) {
			if(replaceReadGroups.size()==0) processFile(alignmentFiles.get(0));
			else processFile(alignmentFiles.get(0),replaceReadGroups.get(0));
			return;
		}
		
		ReadAlignmentFileReader [] readers = new ReadAlignmentFileReader[n];
		ReadAlignment [] currentAlignments = new ReadAlignment[n];
		List<Iterator<ReadAlignment>> iterators = new ArrayList<>();
		try {
			for(int i=0;i<n;i++) {
				readers[i] = createReader(alignmentFiles.get(i));
				Iterator<ReadAlignment> it = readers[i].iterator();
				iterators.add(it);
				if(it.hasNext()) {
					currentAlignments[i] = it.next();
					if(replaceReadGroups.size()>i) (currentAlignments[i]).setReadGroup(replaceReadGroups.get(i)); 
				} else {
					currentAlignments[i] = null;
				}
			}
			boolean querySeqFound = false;
			GenomicRegionComparator cmp = new GenomicRegionComparator(sequencesMetadata);
			while (keepRunning) {
				ReadAlignment aln = chooseNextAln(iterators,currentAlignments, cmp, replaceReadGroups);
				if(aln==null) break;
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
					} else {
						continue;
					}
				}
				//if(aln.getFirst()==187623) System.out.println("Processing alignment "+aln.getReadName()+" group id: "+aln.getReadGroup()+" End: "+aln.getLast()+" CIGAR: "+aln.getCigarString());
				processAlignment(aln);
			}
			if(keepRunning) notifyEndOfAlignments();
			else log.warning("Cancelled process");
			
		} finally {
			for(int i=0;i<n;i++) {
				if(readers[i]!=null)readers[i].close(); 
			}
		}
	}

	private ReadAlignment chooseNextAln(List<Iterator<ReadAlignment>> iterators, ReadAlignment[] currentAlignments, GenomicRegionComparator cmp, List<String> replaceReadGroups) {
		ReadAlignment answer = null;
		int minPos = -1;
		for(int i=0;i<currentAlignments.length;i++) {
			if(currentAlignments[i]!=null) {
				if(answer == null || cmp.compare(answer, currentAlignments[i])>0) {
					answer = currentAlignments[i];
					minPos = i;
				}
			}
		}
		if(minPos>=0) {
			Iterator<ReadAlignment> it = iterators.get(minPos); 
			if(it.hasNext()) {
				currentAlignments[minPos] = it.next();
				if(replaceReadGroups.size()>minPos) (currentAlignments[minPos]).setReadGroup(replaceReadGroups.get(minPos)); 
			} else {
				currentAlignments[minPos] = null;
			}
		}
		return answer;
	}

	public void processFile(String filename) throws IOException {
		processFile(filename, null);
	}
	/**
	 * Process the given alignments file
	 * @param filename Name of the file to process
	 * @param replaceReadGroup Read group to use for the given file. If null, the default read group will be used
	 * @throws IOException If the file can not be read
	 */
	public void processFile(String filename, String replaceReadGroup) throws IOException {		
		int processedAlns = 0;
		try (ReadAlignmentFileReader reader = createReader(filename)) {
			sequencesMetadata = reader.getSequences();
			boolean querySeqFound = false;
			Iterator<ReadAlignment> it = reader.iterator();
			//Sequence under processing
			while(it.hasNext() && keepRunning) {
				ReadAlignment aln = it.next();
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
					} else {
						continue;
					}
				}
				if(replaceReadGroup!=null) aln.setReadGroup(replaceReadGroup);
				processAlignment(aln);
				processedAlns++;
				if(processedAlns%1000000 == 0) log.info("Processed "+processedAlns+" alignments");
			}
		}
		
		if(keepRunning) notifyEndOfAlignments();
		else log.warning("Cancelled process");
	}
	
	private ReadAlignmentFileReader createReader(String filename) throws IOException {
		ReadAlignmentFileReader reader = new ReadAlignmentFileReader(filename);
		//reader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_SEQUENCE);
		reader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_FULL);
		int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
		if(!processSecondaryAlignments ) {
			filterFlags+=ReadAlignment.FLAG_SECONDARY;
			if(!processNonUniquePrimaryAlignments) filterFlags+=ReadAlignment.FLAG_MULTIPLE_ALN;
		}
		reader.setFilterFlags(filterFlags);
		reader.setMinMQ(minMQ);
		return reader;
	}
	
	public void processAlignment(ReadAlignment aln) {
		if(currentReferenceSequence!=null) {
			boolean sameSequence = currentReferenceSequence.getName().equals(aln.getSequenceName());
			if(!sameSequence || lastReadAlignmentStart !=aln.getFirst()) {
				processSameStartAlns();
				if(!sameSequence) {
					processPileups(currentReferenceLast+1);
					for(PileupListener listener:listeners) listener.onSequenceEnd(currentReferenceSequence);
					currentReferenceSequence=null;
				} else {
					processPileups(aln.getFirst());
				}
			}
		} 
		if (currentReferenceSequence==null) {
			startSequence(aln);
			for(PileupListener listener:listeners) listener.onSequenceStart(currentReferenceSequence);
		}
		int alnLast = aln.getLast();
		if(alnLast > currentReferenceLast) currentReferenceLast = alnLast;
		if(aln.isSecondary()) {
			sameStartSecondaryAlignments.add(aln);
		} else {
			sameStartPrimaryAlignments.add(aln);
		}
		lastReadAlignmentStart = aln.getFirst();
	}
	/**
	 * Selects the alignments starting at the same position that will be processed
	 */
	private void processSameStartAlns() {
		int posStart = 0;
		if(sameStartPrimaryAlignments.size()>0) {
			posStart = sameStartPrimaryAlignments.get(0).getFirst();
		} else if(sameStartSecondaryAlignments.size()>0) {
			posStart = sameStartSecondaryAlignments.get(0).getFirst();
		}
		if(posStart == 0 ) return;
		Map<String,Integer> alnsPerReadGroup = new HashMap<>();
		List<ReadAlignment> allAlnsPos = new ArrayList<>();
		allAlnsPos.addAll(sameStartPrimaryAlignments);
		allAlnsPos.addAll(sameStartSecondaryAlignments);
		sameStartPrimaryAlignments.clear();
		sameStartSecondaryAlignments.clear();
		for(ReadAlignment aln:allAlnsPos) {
			Integer count = alnsPerReadGroup.get(aln.getReadGroup());
			if(count == null) {
				alnsPerReadGroup.put(aln.getReadGroup(), 1);
			} else if (maxAlnsPerStartPos<=0 || count<maxAlnsPerStartPos) {
				alnsPerReadGroup.put(aln.getReadGroup(), count+1);
			} else continue;
			aln.setBasesToIgnore5P(basesToIgnore5P);
			aln.setBasesToIgnore3P(basesToIgnore3P);
			pendingAlignments.add(aln);
		}
		allAlnsPos.clear();
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
	}

	
	public void notifyEndOfAlignments() {
		processSameStartAlns();
		processPileups(currentReferenceLast+1);
		if(currentReferenceSequence!=null) for(PileupListener listener:listeners) listener.onSequenceEnd(currentReferenceSequence);
		currentReferenceSequence=null;
	}
	private void processPileups(int alignmentStart) {
		if(alignmentStart==currentReferencePos) return;
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

	

	private boolean processCurrentPosition() {
		int posPrint = -1;
		if(pendingAlignments.size()==0) {
			currentReferencePos++;
			return false;
		}
		if(currentReferencePos==posPrint) System.out.println("Number of pending alignments: "+pendingAlignments.size());
		
		PileupRecord pileup = new PileupRecord(currentReferenceSequence.getName(), currentReferencePos);
		for(ReadAlignment aln:pendingAlignments) {
			if(currentReferencePos==posPrint)System.out.println("Next pending: "+aln.getReadName()+" located at "+aln.getSequenceName()+":"+aln.getFirst()+"-"+aln.getLast()+". group: "+aln.getReadGroup()+" cigar: "+aln.getCigarString());
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

	
}
