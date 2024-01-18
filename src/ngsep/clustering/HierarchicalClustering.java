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
package ngsep.clustering;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;

import java.io.PrintStream;
import java.util.Set;
import java.util.logging.Logger;

/**
 * Class that exposes various distance based clustering algorithms
 * @author Sebastian Lemus
 * @author Jorge Duitama
 */
public class HierarchicalClustering {

    // Constants for command options
    private static final int NJ = 0;
    private static final int FASTNJ = 1;
    private static final int UPGMA = 2;
    private static final Set<Integer> ALG_OPTIONS = Set.of(
            NJ,
            FASTNJ,
            UPGMA
    );

    // Logging and progress
    private Logger log = Logger.getLogger(NeighborJoining.class.getName());
    private ProgressNotifier progressNotifier=null;

    //Parameters
    private String inputFile = null;
    private String outputFile = null;
    private int algorithm = NJ;

    // Get and set methods
    public Logger getLog() {
        return log;
    }
    public void setLog(Logger log) {
        this.log = log;
    }

    public ProgressNotifier getProgressNotifier() {
        return progressNotifier;
    }
    public void setProgressNotifier(ProgressNotifier progressNotifier) {
        this.progressNotifier = progressNotifier;
    }

    public String getInputFile() {
        return inputFile;
    }
    public void setInputFile(String inputFile) {
        this.inputFile = inputFile;
    }

    public String getOutputFile() {
        return outputFile;
    }
    public void setOutputFile(String outputFile) {
        this.outputFile = outputFile;
    }

    public int getAlgorithm () {
        return algorithm;
    }
    public void setAlgorithm (int algorithm) {
        this.algorithm = algorithm;
    }
    public void setAlgorithm(String value) {
		setAlgorithm((int)OptionValuesDecoder.decode(value, Integer.class));
	}


    // Constructor
    public HierarchicalClustering () {

    }

    public static void main(String[] args) throws Exception {
        HierarchicalClustering service = new HierarchicalClustering();
        CommandsDescriptor.getInstance().loadOptions(service, args);
        service.run();
    }

    public void run () throws Exception {
        DistanceMatrix matrix;
        if (inputFile != null) {
            log.info("Loading matrix from file "+ inputFile);
            matrix = new DistanceMatrix(inputFile);
        } else {
            log.info("Loading matrix from standard input");
            matrix = new DistanceMatrix(System.in);
        }


        if (!ALG_OPTIONS.contains(algorithm)) {
            throw new Exception("You must provide with a valid option for the algorithm to be used");
        }
        Dendrogram tree = runClustering(matrix);

        if (outputFile == null) {
            tree.printTree(System.out);
        } else {
            tree.printTree(new PrintStream(outputFile));
        }
    }
	public Dendrogram runClustering(DistanceMatrix matrix) {
		DistanceMatrixClustering method = createClusteringService();
        Dendrogram tree = method.buildDendrogram(matrix);
		return tree;
	}
	private DistanceMatrixClustering createClusteringService() {
		DistanceMatrixClustering method = new NeighborJoining();
        if (algorithm == FASTNJ) {
            method = new FastNJ();
        } else if (algorithm == UPGMA) {
        	method = new UPGMA();
        }
		return method;
	}
}
