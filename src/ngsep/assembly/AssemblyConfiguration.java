/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform Copyright 2016 Jorge
 * Duitama
 *
 * This file is part of NGSEP.
 *
 * NGSEP is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * NGSEP is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * NGSEP. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.assembly;

/**
 * this class contains the configuration variables of the assembly
 * 
 * @author _____________________________
 *
 */
public class AssemblyConfiguration {
	private final static double DEFAULT_RATE_OF_CHANGES = 0.04;
	private final static double DEFAULT_RATE_OF_INDELS = 0.03;

	/**
	 * Configuration of Overlap step
	 */
	private OverlapConfiguration overlapConfigurations;
	/**
	 * Configuration of layout step
	 */
	private LayoutConfiguration layoutConfigurations;
	/**
	 * Configuration of consensus step
	 */
	private ConsensusConfiguration consuensusConfigurations;

	AssemblyConfiguration() {
		this(DEFAULT_RATE_OF_CHANGES, DEFAULT_RATE_OF_INDELS);
	}

	public AssemblyConfiguration(double changes, double indels) {
		overlapConfigurations = new OverlapConfiguration(changes, indels);
		layoutConfigurations = new LayoutConfiguration(changes, indels);
		consuensusConfigurations = new ConsensusConfiguration(changes, indels);
	}

	class OverlapConfiguration {
		private static final double LN100000 = 11.51292546; // Natural logarithm of 10000
		private static final double LN2 = 0.693147806;// Natural logarithm of 2
		/**
		 * the length of the kmers
		 */
		private int kmerLength;
		/**
		 * distance between kmers
		 */
		private int KmerDistance;
		/**
		 * maximum difference of relative position between two hits of the same align
		 */
		private int maxKmerDiff;
		/**
		 * the lowest rate of kmers to be considered an align
		 */
		private double minKmerCoverRate;
		/**
		 * the mean of kmers over a letter in the sequence
		 */
		private double rate_of_cover;

		public OverlapConfiguration(double changes, double indels) {
			// Supposing independence
			double rate_of_error = changes + indels - changes * indels;
			kmerLength = (int) (1.25 * (LN2 / rate_of_error));
			maxKmerDiff = (int) (indels * (LN100000 / rate_of_error));
			// empirical
			
			rate_of_cover = Math.max(1, rate_of_error * 20);
			KmerDistance = (int) (kmerLength * ((1 / rate_of_cover) - 1));
			// poison of 0 errors
			minKmerCoverRate = 0.8 * Math.exp(-2 * kmerLength * rate_of_error);

		}

		/**
		 * @return distance between kmers
		 */
		public int getKmerDistance() {
			return KmerDistance;
		}

		/**
		 * @return the length of the kmers
		 */
		public int getKmerLength() {
			return kmerLength;
		}

		/**
		 * @return maximum difference of relative position between two hits of the same
		 *         align
		 */
		public int getMaxKmerDiff() {
			return maxKmerDiff;
		}

		/**
		 * @return the lowest rate of kmers to be considered an align
		 */
		public double getMinKmerCoverRate() {
			return minKmerCoverRate;
		}

		/**
		 * @return the mean of kmers over a letter in the sequence
		 */
		public double getRate_of_cover() {
			return rate_of_cover;
		}
	}

	public OverlapConfiguration overlap() {
		return overlapConfigurations;
	}

	class LayoutConfiguration {
		LayoutConfiguration(double changes, double indels) {

		}

	}

	public LayoutConfiguration layout() {
		return layoutConfigurations;
	}

	class ConsensusConfiguration {
		public ConsensusConfiguration(double changes, double indels) {
		}
	}

	public ConsensusConfiguration consuensus() {
		return consuensusConfigurations;
	}
}
