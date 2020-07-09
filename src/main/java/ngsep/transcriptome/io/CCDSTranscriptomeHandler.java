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
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.transcriptome.TranscriptSegment;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.Transcriptome;

/**
 * Implementation of an AbstractTranscriptome that loads the transcripts information
 * from a CCDS file. Uses a configurable prefix for sequence names to make them consistent
 * with other reference assemblies as the ones in UCSC 
 * @author Jorge Duitama
 *
 */
public class CCDSTranscriptomeHandler {
	//Prefix for sequenceNames
	private String seqNamePrefix = "chr";
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
	 * @return Transcriptome object with the transcriptome map
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
					boolean reverse = items[6].charAt(0)=='-';
					if(!items[5].startsWith("Withdrawn") && !items[7].equals("-") ) {
						String newSeqName = seqNamePrefix+items[0];
						QualifiedSequence seq;
						try {
							seq = sequenceNames.addOrLookupName(newSeqName);
						} catch (RuntimeException ex) {
							System.err.println("Can not load line: "+line+". Unrecognized sequence name"+ex.getMessage());
							line=in.readLine();
							continue;
						}
						int first = Integer.parseInt(items[7])+1;
						int last = Integer.parseInt(items[8])+1;
						Gene gene =  answer.getGene(items[3]);
						if(gene == null) {
							gene = new Gene(items[3],items[2],seq.getName(),first, last, reverse);
						} else {
							if(gene.getLast()<last) gene.setLast(last);
							if(gene.getFirst()>first) gene.setFirst(first);
						}
						
						//Add 1 to all starts and ends because CCDS positions are zero indexed
						Transcript transcript = new Transcript(items[4],seq.getName(),first,last,reverse);
						String exonsStr = items[9].substring(1,items[9].length()-1);
						String [] items2 = exonsStr.split(", ");
						List<TranscriptSegment> tExons = new ArrayList<TranscriptSegment>();
						for(int i=0;i<items2.length;i++) {
							int index2 = items2[i].indexOf('-');
							int exonStart = Integer.parseInt(items2[i].substring(0,index2))+1;
							int exonEnd = Integer.parseInt(items2[i].substring(index2+1))+1;
							TranscriptSegment exon = new TranscriptSegment(transcript, exonStart,exonEnd);
							tExons.add(exon);
						}
						transcript.setTranscriptSegments(tExons);
						answer.addTranscript(transcript);
					}
				}
				line=in.readLine();
			}
		}
		return answer;
	}
	/**
	 * Loads the transcripts sequences from the given fasta file 
	 * @param filename Name of the fasta file with the transcript sequences
	 * @throws IOException If the file can not be read
	 */
	public void loadSequences(Transcriptome transcriptome, String filename) throws IOException {
		FastaSequencesHandler sequencesHandler = new FastaSequencesHandler();
		sequencesHandler.setSequenceType(DNAMaskedSequence.class);
		QualifiedSequenceList sequences = sequencesHandler.loadSequences(filename);
		
		int n = sequences.size();
		for(int i=0;i<n;i++ ) {
			QualifiedSequence sequence = sequences.get(i);
			Transcript transcript = transcriptome.getTranscript(sequence.getName());
			if(transcript!=null) {
				transcript.setCDNASequence((DNAMaskedSequence)sequence.getCharacters());
			} else {
				System.err.println("WARN: Transcript not found with id: "+sequence.getName());
			}
		}
	}
	/**
	 * Cuts the transcript id from the id used in a fasta file
	 * @param sequenceName Name of the sequence in the fasta file
	 * @return String Transcript id in the map file
	 */
	protected String getTranscriptId(String fastaId) {
		int index = fastaId.indexOf('|');
		return fastaId.substring(0,index);
	}
}
