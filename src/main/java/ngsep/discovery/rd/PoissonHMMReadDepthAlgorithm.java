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

import JSci.maths.statistics.PoissonDistribution;
import ngsep.hmm.HMMState;
import ngsep.math.LogMath;

public class PoissonHMMReadDepthAlgorithm extends AbstractHMMReadDepthAlgorithm {
	
	public static final String SOURCE_POISSONHMM = "POISSONHMM";
	@Override
	protected String getSource() {
		return SOURCE_POISSONHMM;
	}
	@Override
	protected HMMState createHMMState(int copies, Double logStart) {
		double avgNormalDepth = this.getReadDepthDistribution().getMeanReadDepth();
		double avgDepthState = avgNormalDepth*copies/getNormalPloidy();
		if(copies==0) avgDepthState = 1;
		HMMState state = new PoissonHMMState(copies, avgDepthState, logStart);
		//System.out.println("Created state "+state.getId()+" with average depth "+avgDepthState+" log start "+logStart+" emission 20 reads: "+state.getEmission(20.0, 0));
		return state; 
	}
}
class PoissonHMMState implements HMMState {

	private int copies;
	private double averageDepth;
	private Double logStart;
	
	
	/**
	 * @param copies
	 * @param averageDepth
	 * @param logStart
	 */
	public PoissonHMMState(int copies, double averageDepth, Double logStart) {
		super();
		this.copies = copies;
		this.averageDepth = averageDepth;
		this.logStart = logStart;
	}

	@Override
	public Double getEmission(Object value, int step) {
		if(value == null || !(value instanceof Double)) return null;
		double depth = (Double)value;
		if(depth<1) depth = 1;
		PoissonDistribution dist = new PoissonDistribution(averageDepth);
		double a = dist.probability(depth);
		//System.out.println("--- depthPoisson ---- " + depth + " ----a--- " + a);
		// double p = dist.cumulative(depth+0.5)-dist.cumulative(depth-0.5);
		// if(copies==0 && p<0.00001) System.out.println("Emission prob "+p+" cumulative 1: "+dist.cumulative(depth-0.05)+"cumulative 2 "+dist.cumulative(depth+0.05)+" depth "+depth);
		return LogMath.log10(a);
	}

	@Override
	public Double getLogStart() {
		return logStart;
	}

	@Override
	public void setLogStart(Double logStart) {
		this.logStart = logStart;
	}

	@Override
	public String getId() {
		return ""+copies;
	}
	
}
