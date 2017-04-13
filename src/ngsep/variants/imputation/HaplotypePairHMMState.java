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
package ngsep.variants.imputation;

import ngsep.hmm.HMMState;
import ngsep.math.LogMath;
import ngsep.variants.CalledSNV;

public class HaplotypePairHMMState implements HMMState {
	private String id = null;
	private int index1;
	private int index2;
	private HaplotypeClusterHMMState state1;
	private HaplotypeClusterHMMState state2;
	private Double logStart=null;
	
	public HaplotypePairHMMState(int index1, HaplotypeClusterHMMState state1, int index2, HaplotypeClusterHMMState state2) {
		super();
		this.index1 = index1;
		this.state1 = state1;
		this.index2 = index2;
		this.state2 = state2;
	}

	@Override
	public Double getEmission(Object value, int step) {
		Byte genotype = getGenotype (value);
		Double answer = null;
		if(genotype!=null) {
			byte a0 = 0;
			byte a1 = 1;
			if(genotype==CalledSNV.GENOTYPE_HOMOREF) answer = LogMath.logProduct(state1.getEmission(a0, step), state2.getEmission(a0, step));
			else if(genotype==CalledSNV.GENOTYPE_HOMOALT) answer = LogMath.logProduct(state1.getEmission(a1, step), state2.getEmission(a1, step));
			else if(genotype==CalledSNV.GENOTYPE_HETERO) {
				double p1 = LogMath.logProduct(state1.getEmission(a0, step), state2.getEmission(a1, step));
				double p2 = LogMath.logProduct(state1.getEmission(a1, step), state2.getEmission(a0, step));
				answer = LogMath.logSum(p1, p2);
			}
		}
		
		if(answer == null) {
			answer = HaplotypeClusterHMMState.LOGPROB_UNEXPECTED;
		} else {
			answer = LogMath.logProduct(answer, HaplotypeClusterHMMState.LOGPROB_EXPECTED);
		}
		return answer;
	}

	private Byte getGenotype(Object value) {
		Byte answer = null;
		if(value == null) return answer;
		if(value instanceof Byte) {
			answer = (Byte)value;
		} else if (value instanceof CalledSNV) {
			answer = ((CalledSNV)value).getGenotype();
		}
		return answer;
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
		return id;
	}
	
	public void setId(String id) {
		this.id = id;
	}

	public int getIndex1() {
		return index1;
	}

	public int getIndex2() {
		return index2;
	}

	public HaplotypeClusterHMMState getState1() {
		return state1;
	}

	public HaplotypeClusterHMMState getState2() {
		return state2;
	}
	
	

}
