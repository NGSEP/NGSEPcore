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

import java.util.List;
import java.util.Map;

import ngsep.hmm.HMMState;
import ngsep.variants.CalledSNV;

public class DiploidGenotypeImputationHMM extends GenotypeImputationHMM {

	public DiploidGenotypeImputationHMM(List<? extends HMMState> states, int numMarkers) {
		super(states, numMarkers);
		// TODO Auto-generated constructor stub
	}

	@Override
	protected void runBaumWelchStep(Map<String, List<CalledSNV>> genotypes,boolean forceHaploid) {
		if(forceHaploid) {
			super.runBaumWelchStep(genotypes,forceHaploid);
			return;
		}
		//TODO: Implement training from diploid data
		
	}

	@Override
	protected void calculatePosteriors(List<CalledSNV> genotypes, double[][] genotypePosteriors) {
		// TODO Auto-generated method stub
		super.calculatePosteriors(genotypes, genotypePosteriors);
	}

	
}
