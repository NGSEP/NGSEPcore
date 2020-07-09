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
package ngsep.genome;
/**
 * Interface for regions mapped to a reference genome. Regions can be of many types
 * such as SNPs, Indels, CNVs, Exons, Transcripts, Genes, Contigs or any feature for which 
 * limits within a genome can be defined.
 * This common interface helps to implement common tasks like sorting or overlaps finding 
 * @author Jorge Duitama
 *
 */
public interface GenomicRegion {
	/**
	 * @return String name of the sequence where the region is mapped
	 */
	public String getSequenceName();
	/**
	 * @return One-based first position of the region relative to the sequence
	 */
	public int getFirst();
	/**
	 * @return One-based last position of the region relative to the sequence
	 */
	public int getLast(); //Last position of the region always greater than start
	/**
	 * @return int length of the region. Typically equal to last - first + 1
	 */
	public int length();
	/**
	 * @return boolean true if the feature lies in the positive strand of the sequence.
	 * false if the feature lies in the negative strand or if the strand orientation is
	 * not relevant to the feature
	 */
	public boolean isPositiveStrand();
	/**
	 * @return boolean true if the feature lies in the negative strand of the sequence.
	 * false if the feature lies in the positive strand or if the strand orientation is
	 * not relevant to the feature
	 */
	public boolean isNegativeStrand();
}
