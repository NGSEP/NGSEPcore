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

	@Override
	public List<QualifiedSequence> makeConsensus(AssemblyGraph graph) 
	{
		//List of final contigs
		List<QualifiedSequence> consensusList = new ArrayList<QualifiedSequence>();
		List<List<AssemblyEdge>> paths = graph.getPaths(); 
		for(int i = 0; i < paths.size(); i++)
		{
			List<AssemblyEdge> path = paths.get(i);
			String sequenceName = sequenceNamePrefix+"_"+(i+1);
			log.info("processing path: "+sequenceName+" size: "+path.size());
			CharSequence consensusPath = makeConsensus (graph, path, i, sequenceName);
			log.info("processed path: "+sequenceName+" consensus length: "+consensusPath.length());
			consensusList.add(new QualifiedSequence(sequenceName,consensusPath));
		}
		
		return consensusList;
	}
	
	public CharSequence makeConsensus(AssemblyGraph graph, List<AssemblyEdge> path, int sequenceIdx, String sequenceName) 
	{
		AssemblyPathReadsAligner aligner = new AssemblyPathReadsAligner();
		aligner.setLog(log);
		aligner.setOnlyGenerateConsensus(true);
		aligner.alignPathReads(graph, path, sequenceIdx);
		return aligner.getConsensus().toString();
	}
}
