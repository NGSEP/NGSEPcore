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
package ngsep.discovery.rd;

import java.util.List;
import java.util.logging.Logger;

import ngsep.variants.CalledCNV;

public interface SingleSampleReadDepthAlgorithm {
	/**
	 * Changes the object for logging
	 * @param log New logger
	 */
	public void setLog(Logger log);
	/**
	 * Changes the genome size assumed to make calculations
	 * @param genomeSize New genome size
	 */
	public void setGenomeSize(long genomeSize);
	/**
	 * Changes the ploidy in non-CNV regions
	 * @param normalPloidy New normal ploidy
	 */
	public void setNormalPloidy(byte normalPloidy);
	/**
	 * Changes the input read depth distribution
	 * @param distribution New distribution
	 */
	public void setReadDepthDistribution(ReadDepthDistribution distribution);
	/**
	 * Runs the read depth algorithm to find CNVs
	 * @return List<CalledCNV> CNVs identified by the algorithm
	 */
	public List<CalledCNV> callCNVs();
	
}
