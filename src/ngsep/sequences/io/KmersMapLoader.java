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
package ngsep.sequences.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.logging.Logger;

import ngsep.main.io.ConcatGZIPInputStream;
import ngsep.sequences.DefaultKmersMapImpl;
import ngsep.sequences.KmersMap;
import ngsep.sequences.ShortArrayDNAKmersMapImpl;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class KmersMapLoader {
	private Logger log = Logger.getLogger(KmersMapLoader.class.getName());
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}



	public KmersMap loadKmersMap(String kmersMapFile, int kmerLength) throws IOException {
		KmersMap kmersMap;
		log.info("Loading k-mers map from : "+kmersMapFile);
		if(kmerLength<=15) kmersMap = new ShortArrayDNAKmersMapImpl((byte) kmerLength);
		else kmersMap = new DefaultKmersMapImpl();
		try (FileInputStream fis = new FileInputStream(kmersMapFile)) {
			InputStream is=fis;
			if(kmersMapFile.toLowerCase().endsWith(".gz")) {
				is = new ConcatGZIPInputStream(is);
			}
			try (BufferedReader in = new BufferedReader(new InputStreamReader(is))) {
				String line = in.readLine();
				while(line!=null) {
					String [] items = line.split("\t| ");
					String kmer = items[0];
					int count = Integer.parseInt(items[1]);
					kmersMap.setCount(kmer,count);
					line = in.readLine();
				}
			}
		}
		log.info("Extracted "+kmersMap.size()+" k-mers from: " + kmersMapFile);
		return kmersMap;
	}
}
