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
package ngsep.sequences;

import java.util.Arrays;

public class MinimapShortKmerCodesSamplingAlgorithm implements ShortKmerCodesSamplingAlgorithm {
	private int windowLength = ShortKmerCodesSampler.DEF_WINDOW_LENGTH;

	public int getWindowLength() {
		return windowLength;
	}
	public void setWindowLength(int windowLength) {
		this.windowLength = windowLength;
	}

	@Override
	public boolean[] sample(Integer[] hashcodes) {
		boolean[] selected = new boolean[hashcodes.length];
		Arrays.fill(selected, false);
		Integer previousMinimizer = null;
		int previousMinimizerPos = -1;
		for(int i=0;true;i++) {
			int posLastKmerWindow = i+windowLength-1;
			if(posLastKmerWindow>=hashcodes.length) break;
			Integer minimizerI = null;
			int minPos = -1;
			Integer newHash = hashcodes[posLastKmerWindow];
			boolean lastInRange = previousMinimizer!=null && previousMinimizerPos>=i;
			if(lastInRange && (newHash==null || previousMinimizer < newHash)) {
				minimizerI = previousMinimizer;
				minPos = previousMinimizerPos;
			} else if (newHash!=null && (previousMinimizer==null || newHash <= previousMinimizer)) {
				minimizerI = newHash;
				minPos = i+windowLength-1;
			}
			if(minimizerI == null) {
				for(int j=0;j<windowLength;j++) {
					Integer hash = hashcodes[i+j];
					if (hash!=null && (minimizerI==null || hash <= minimizerI)) {
						minimizerI = hash;
						minPos = i+j;
					}
				}
				//if(sequenceId==debugIdx && i>0 && i<3000) System.err.println("Minimizer calculated with cycle. Start: "+i+" New pos: "+minPos+" new minimizer: "+minimizerI+" previous: "+previousMinimizer+" kmer code: "+kmerCodes.get(minPos)+" total: "+minimizersSeq.size());
			}
			if (minimizerI==previousMinimizer) continue;
			if(minimizerI != null) {
				selected[minPos] = true;
			}
			previousMinimizer = minimizerI;
			previousMinimizerPos = minPos;
		}
		return selected;
	}

}
