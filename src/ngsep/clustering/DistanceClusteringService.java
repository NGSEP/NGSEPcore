package ngsep.clustering;

import ngsep.clustering.dendrogram.Dendrogram;
import ngsep.clustering.nj.FastNJ;
import ngsep.clustering.nj.NeighborJoining;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;

import java.io.PrintStream;
import java.util.Set;
import java.util.logging.Logger;

/**
 * Class that exposes various distance based clustering algorithms
 */
public class DistanceClusteringService {

    // Constants for command options
    private static final int NJ = 0;
    private static final int FASTNJ = 1;
    private static final Set<Integer> ALG_OPTIONS = Set.of(
            NJ,
            FASTNJ
    );

    // Logging and progress
    private Logger log = Logger.getLogger(NeighborJoining.class.getName());
    private ProgressNotifier progressNotifier=null;

    //Parameters
    private String inputFile = null;
    private String outputFile = null;
    private int algorithm = 0;

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
    public void setAlgorithm (Integer algorithm) {
        this.algorithm = algorithm;
    }


    // Constructor
    public DistanceClusteringService () {

    }

    public static void main(String[] args) throws Exception {
        DistanceClusteringService service = new DistanceClusteringService();
        CommandsDescriptor.getInstance().loadOptions(service, args);
        service.run();
    }

    private void run () throws Exception {
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
        DistanceMatrixClustering method = new NeighborJoining();
        if (algorithm == FASTNJ) {
            method = new FastNJ();
        }

        Dendrogram tree = method.buildDendrogram(matrix);

        if (outputFile == null) {
            tree.printTree(System.out);
        } else {
            tree.printTree(new PrintStream(outputFile));
        }
    }
}
