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
package ngsep.transcriptome.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.transcriptome.TranscriptSegment;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.Transcriptome;


/**
 * Implementation of an AbstractTranscriptome that loads the transcripts information
 * from a knownGenes file from UCSC. 
 * @author Jorge Duitama
 *
 */
public class UCSCTranscriptomeHandler {
	private QualifiedSequenceList sequenceNames = new QualifiedSequenceList();
	public QualifiedSequenceList getSequenceNames() {
		return sequenceNames;
	}
	public void setSequenceNames(QualifiedSequenceList sequenceNames) {
		if(sequenceNames==null) this.sequenceNames = new QualifiedSequenceList();
		else this.sequenceNames = sequenceNames;
	}
	/**
	 * Loads the genes, transcripts and exons information from the given map file
	 * @param filename Name of the file with the transcriptome description
	 * @throws IOException If the file can not be read
	 */
	public Transcriptome loadMap(String filename) throws IOException {
		Transcriptome answer = new Transcriptome(sequenceNames);
		try (FileInputStream fis = new FileInputStream(filename);
			 BufferedReader in = new BufferedReader(new InputStreamReader(fis)) ) {
			String line=in.readLine();
			while(line!=null) {
				if(line.charAt(0)!='#') {
					String [] items = line.split("\t");
					String newSeqName = items[1]; 
					QualifiedSequence seq;
					try {
						seq = sequenceNames.addOrLookupName(newSeqName);
					} catch (RuntimeException e) {
						System.err.println("Can not load line: "+line+". Unrecognized sequence name. "+e.getMessage());
						line=in.readLine();
						continue;
					}
					boolean reverse = items[2].charAt(0)=='-';
					//Add 1 to all starts because positions are zero indexed and ends are exclusive
					int start = Integer.parseInt(items[3])+1;
					int end = Integer.parseInt(items[4]);
					int translationStart = Integer.parseInt(items[5])+1;
					int translationEnd = Integer.parseInt(items[6]);
					boolean coding = translationEnd > translationStart;
					String geneId = items[0];
					String geneId2 = null;
					if(items.length>10) {
						geneId2 = items[10];
						if(geneId2.length()>0 && !"n/a".equals(geneId2) ) {
							geneId = geneId2; 
						}
					}
					
					Gene gene = answer.getGene(geneId);
					if(gene == null) {
						gene = new Gene(geneId,geneId,seq.getName(), start,end,reverse);
					} else {
						if(gene.getLast()<end) gene.setLast(end);
						if(gene.getFirst()>start) gene.setFirst(start);
					}
					
					
					
					Transcript transcript = new Transcript(items[0],seq.getName(),start,end,reverse);
					transcript.setGene(gene);
					String [] exonStarts = items[8].split(",");
					String [] exonEnds = items[9].split(",");
					List<TranscriptSegment> tExons = new ArrayList<TranscriptSegment>();
					for(int i=0;i<exonStarts.length;i++) {
						int exonStart = Integer.parseInt(exonStarts[i])+1;
						int exonEnd = Integer.parseInt(exonEnds[i]);
						if(coding) {
							//Split exon if it includes non coding and coding region
							if(translationStart>exonStart && translationStart<exonEnd) {
								TranscriptSegment exon = new TranscriptSegment(transcript,exonStart,translationStart-1);
								if(reverse) {
									exon.setStatus(TranscriptSegment.STATUS_3P_UTR);
								} else {
									exon.setStatus(TranscriptSegment.STATUS_5P_UTR);
								}
								tExons.add(exon);
								exonStart = translationStart;
							}
							if(translationEnd>exonStart && translationEnd<exonEnd) {
								TranscriptSegment exon=new TranscriptSegment(transcript, exonStart, translationEnd);
								exon.setStatus(TranscriptSegment.STATUS_CODING);
								tExons.add(exon);
								exonStart = translationEnd+1;
							}
							TranscriptSegment exon = new TranscriptSegment(transcript,exonStart,exonEnd);
							if(exonEnd < translationStart ) {
								if(reverse) {
									exon.setStatus(TranscriptSegment.STATUS_3P_UTR);
								} else {
									exon.setStatus(TranscriptSegment.STATUS_5P_UTR);
								}
							} else if (exonStart>translationEnd) {
								if(reverse) {
									exon.setStatus(TranscriptSegment.STATUS_5P_UTR);
								} else {
									exon.setStatus(TranscriptSegment.STATUS_3P_UTR);
								}
							} else {
								exon.setStatus(TranscriptSegment.STATUS_CODING);
							}
							tExons.add(exon);
						} else {
							TranscriptSegment exon = new TranscriptSegment(transcript,exonStart,exonEnd);
							exon.setStatus(TranscriptSegment.STATUS_NCRNA);
							tExons.add(exon);
						}
					}
					transcript.setTranscriptSegments(tExons);
					answer.addTranscript(transcript);
				}
				line=in.readLine();
			}
		}
		
		return answer;
	}
	public void loadSequences(Transcriptome transcriptome, String mrnaFilename) throws IOException {
		FileInputStream fis = new FileInputStream(mrnaFilename);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis)); 
		String line=in.readLine();
		while(line!=null) {
			String [] items = line.split("\t");
			Transcript transcript = transcriptome.getTranscript(items[0]);
			if(transcript!=null) {
				transcript.setCDNASequence(new DNAMaskedSequence(items[1].toUpperCase()));
			} else {
				System.err.println("WARN: Transcript not found with id: "+items[0]);
			}
			line=in.readLine();
		}
		fis.close();
	}

}
