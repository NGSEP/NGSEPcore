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
package ngsep.transcriptome;


import ngsep.genome.GenomicRegion;

/**
 * Information for an Exon in a transcript. If the same exon is spanned by two different transcripts, two different objects
 * should be created
 * @author Jorge Duitama
 */
public class TranscriptSegment implements GenomicRegion {
	public static final int STATUS_CODING = 0;
	public static final int STATUS_5P_UTR = 1;
	public static final int STATUS_3P_UTR = 2;
	public static final int STATUS_NCRNA = 3;
	private Transcript transcript;
	private int first;
	private int last; //Always greater or equal than start
	private int status = STATUS_CODING;
	/**
	 * Creates a new Exon with the given information
	 * @param gene Gene where the exon is located
	 * @param first position of the exon in the sequence where the gene is located
	 * @param last position of the exon in the sequence where the gene is located
	 */
	public TranscriptSegment(Transcript t, int first, int last) {
		this.transcript = t;
		this.setFirst(first);
		this.setLast(last);
	}
	/**
	 * @return int First position of the exon in the sequence where the gene is located
	 */
	public int getFirst() {
		return first;
	}
	/**
	 * Changes the start of the exon
	 * @param start New Start
	 */
	public void setFirst(int first) {
		this.first = first;
	}
	/**
	 * @return int Last position of the exon in the sequence where the gene is located
	 */
	public int getLast() {
		return last;
	}
	/**
	 * Changes the end of the exon
	 * @param end New end
	 */
	public void setLast(int last) {
		this.last = last;
	}
	/**
	 * @return boolean true if the exon is coding
	 */
	public boolean isCoding() {
		return status == STATUS_CODING;
	}
	/**
	 * @return int length of the exon calculated as last-first+1 
	 */
	public int length() {
		return last-first+1;
	}
	/**
	 * @return String Sequence name where the gene is located
	 */
	public String getSequenceName() {
		return transcript.getSequenceName();
	}
	/**
	 * @return int Exon status
	 */
	public int getStatus() {
		return status;
	}
	/**
	 * Changes the exon status
	 * @param status new status
	 */
	public void setStatus(int status) {
		this.status = status;
	}
	/**
	 * @return Transcript where this exon is included
	 */
	public Transcript getTranscript() {
		return transcript;
	}
	@Override
	public boolean isPositiveStrand() {
		return transcript.isPositiveStrand();
	}
	@Override
	public boolean isNegativeStrand() {
		return transcript.isNegativeStrand();
	}
	
}
