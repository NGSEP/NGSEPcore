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
package ngsep.assembly;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import ngsep.sequences.QualifiedSequence;

/**
 * @author David Guevara
 * @author Jorge Duitama
 *
 */
public class ConsensusBuilderBidirectionalSimple implements ConsensusBuilder {
	
	private Logger log = Logger.getLogger(ConsensusBuilderBidirectionalWithPolishing.class.getName());
	private String sequenceNamePrefix = "Contig";
	private int numThreads = 1;
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}
	
	public String getSequenceNamePrefix() {
		return sequenceNamePrefix;
	}

	public void setSequenceNamePrefix(String sequenceNamePrefix) {
		this.sequenceNamePrefix = sequenceNamePrefix;
	}
	
	

	public int getNumThreads() {
		return numThreads;
	}

	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}

	@Override
	public List<QualifiedSequence> makeConsensus(AssemblyGraph graph) 
	{
		AssemblyPathReadsAligner aligner = new AssemblyPathReadsAligner();
		aligner.setLog(log);
		//List of final contigs
		List<QualifiedSequence> consensusList = new ArrayList<QualifiedSequence>();
		List<AssemblyPath> paths = graph.getPaths(); 
		for(int i = 0; i < paths.size(); i++)
		{
			AssemblyPath path = paths.get(i);
			String sequenceName = sequenceNamePrefix+"_"+(i+1);
			path.setPathId(i+1);
			path.setSequenceName(sequenceName);
			log.info("processing path: "+sequenceName+" size: "+path.getPathLength());
			
			aligner.calculateConsensus(path);
			CharSequence consensusPath = path.getConsensus();
			log.info("processed path: "+sequenceName+" consensus length: "+consensusPath.length());
			consensusList.add(new QualifiedSequence(sequenceName,consensusPath));
		}
		
		return consensusList;
	}
}
