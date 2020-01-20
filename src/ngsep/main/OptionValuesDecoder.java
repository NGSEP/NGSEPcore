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
package ngsep.main;

import java.io.IOException;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;

public class OptionValuesDecoder {
	public static Object decode (String value, Class<?> type) {
		if(Integer.class.equals(type)) {
			return Integer.parseInt(value);
		}
		if(Long.class.equals(type)) {
			return Long.parseLong(value);
		}
		if(Short.class.equals(type)) {
			return Short.parseShort(value);
		}
		if(Byte.class.equals(type)) {
			return Byte.parseByte(value);
		}
		if(Float.class.equals(type)) {
			return Float.parseFloat(value);
		}
		if(Double.class.equals(type)) {
			return Double.parseDouble(value);
		}
		if(String.class.equals(type)) {
			return value;
		}
		throw new IllegalArgumentException("Can not decode value of type: "+type.toString());
	}

	public static ReferenceGenome loadGenome(String genomeFile, Logger log) throws IOException {
		log.info("Loading genome from: "+genomeFile);
		ReferenceGenome genome = new ReferenceGenome(genomeFile);
		log.info("Loaded genome with: "+genome.getNumSequences()+" sequences. Total length: "+genome.getTotalLength()+" from file: "+genomeFile);
		return genome;
	}
}
