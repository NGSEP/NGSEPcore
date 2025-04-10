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
package ngsep.math;

/**
 * @author Nicolas Rozo Fajardo
 */

public interface EntropyCalculator {

    public double calculateEntropy(CharSequence sequence);
    public double getMaxEntropy();

    /**
     * Method for normalizing the entropy value
     * 
     * @param entropy Entropy value to be normalized
     * @return Entropy value normalized (using maximum entropy value)
     */
    default double normalizeEntropy(double entropy) {
        double maxEntropy = getMaxEntropy();
        return (maxEntropy == 0) ? 0 : entropy / maxEntropy;
    }

    /**
     * Method for denormalizing the entropy value
     * 
     * @param normalizedEntropy Normalized entropy value (using maximum entropy value)
     * @return Original entropy value, without normalization
     */
    default double denormalizeEntropy(double normalizedEntropy) {
        return normalizedEntropy * getMaxEntropy();
    }
}
