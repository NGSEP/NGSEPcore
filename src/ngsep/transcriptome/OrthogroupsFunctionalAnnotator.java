package ngsep.transcriptome;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ngsep.main.CommandsDescriptor;

public class OrthogroupsFunctionalAnnotator {
	  int[] functionalNotationColums = new int[] { 2, 6 };
	  private final int KEGG_COLUMN = 11;
	  private final String USELESS_NOTATION_STRING = ".";
	
	  Map<String, Set<String>> genMapNotations = new HashMap<>();
	  Map<String, Set<String>> clusterMapNotations = new HashMap<>();

	  Map<String, Set<String>> keggNotationsByGen = new HashMap<>();
	  Map<String, Set<String>> clusterMapKeggNotations = new HashMap<>();
	

	  private String clustersFilePath;
	  private String fileOut;
	  private String functionalNotationFilePath;
	  
	  public String getFunctionalNotationFilePath() {
	    return functionalNotationFilePath;
	  }

	  public void setFunctionalNotationFilePath(String functionalNotationFilePath) {
	    this.functionalNotationFilePath = functionalNotationFilePath;
	  }

	  public String getClustersFilePath() {
	    return clustersFilePath;
	  }

	  public void setClustersFilePath(String clustersFilePath) {
	    this.clustersFilePath = clustersFilePath;
	  }

	  public String getFileOut() {
	    return fileOut;
	  }

	  public void setFileOut(String fileOut) {
	    this.fileOut = fileOut;
	  }

	  public void run() {
	    try (BufferedReader functionalNotationFile = new BufferedReader(new FileReader(functionalNotationFilePath));
	        BufferedReader clustersFile = new BufferedReader(new FileReader(clustersFilePath))) {

			processFunctionalAnotationFile(functionalNotationFile);
	      	processClustersFile(clustersFile);
	    } catch (IOException e) {
	      e.printStackTrace();
	    }

	    writeClusterNotationCsv(fileOut);
	  }

	  private void processClustersFile(BufferedReader clustersFile) throws IOException {
	    String line;
	    while ((line = clustersFile.readLine()) != null) {
	      String[] data = line.split("\t");
	      String cluster = data[0];
	      for (int i = 1; i < data.length; i++) {
	        String genName = data[i];

			Set<String> currentGenNotation = genMapNotations.getOrDefault(genName, Collections.emptySet());
	        if (!currentGenNotation.isEmpty()) {
	          clusterMapNotations.computeIfAbsent(cluster, k -> new HashSet<>()).addAll(currentGenNotation);
	        }

			Set<String> currentGenKeggNotation = keggNotationsByGen.getOrDefault(genName, Collections.emptySet());
			if (!currentGenKeggNotation.isEmpty()) {
	          clusterMapKeggNotations.computeIfAbsent(cluster, k -> new HashSet<>()).addAll(currentGenKeggNotation);
	        }
	      }
	    }
	  }
	  
	  private Set<String> createKEGGNotationsSet(String keggNotationLine) {
		Set<String> keggNotations = new HashSet<>();
	    String[] keggNotationArray = keggNotationLine.split("`");
	    for (String keggNotation : keggNotationArray) {
	      keggNotations.add(keggNotation);
	    }
	    
	    return keggNotations;
	  }

	  private void processFunctionalAnotationFile(BufferedReader functionalNotationFile) throws IOException {
	    String line;
	    functionalNotationFile.readLine(); // Skip header
	    while ((line = functionalNotationFile.readLine()) != null) {
	      String[] data = line.split("\t");
	      String gen = data[0];
	      Set<String> functionalNotations = new HashSet<>();

	      for (int columNumber : functionalNotationColums) {
	        String functionalNotation = data[columNumber].trim();

	        if (!functionalNotation.equals(USELESS_NOTATION_STRING)) {
	          functionalNotations.add(cleanBlastResult(functionalNotation));
	        }
	      }

	      if (!functionalNotations.isEmpty()) {
	        genMapNotations.put(gen, functionalNotations);
	      }
	      
	      String KEGGNotation = data[KEGG_COLUMN].trim();
	      if (!KEGGNotation.equals(USELESS_NOTATION_STRING)) {
			  keggNotationsByGen.put(gen, createKEGGNotationsSet(KEGGNotation));
	      }
	    }
	  }

	  private String cleanBlastResult(String blast) {
	    if (blast.equals(".")) {
	      return "";
	    }
	    Matcher matcher = Pattern.compile("RecName: .*?;").matcher(blast);
	    if (matcher.find()) {
	      return matcher.group(0).substring(14);
	    } else {
	      throw new IllegalArgumentException(blast);
	    }
	  }

	  private void writeClusterNotationCsv(String fileOut) {
	    try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(fileOut))) {
	      for (Map.Entry<String, Set<String>> entry : clusterMapNotations.entrySet()) {
	        String cluster = entry.getKey();
	        Set<String> notations = entry.getValue();

	        StringBuilder allNotationsString = new StringBuilder();
	        for (String notation : notations) {
	          allNotationsString.append(notation).append(",");
	        }

			StringBuilder allKeggNotationsString = new StringBuilder();
			Set<String> keggSet =  clusterMapKeggNotations.get(cluster);
			if (keggSet != null) {
				for (String notation : keggSet) {
					allKeggNotationsString.append(notation).append(",");
				}
			}

	        writer.write(cluster + "\t" + allNotationsString + "\t" + allKeggNotationsString);
	        writer.newLine();
	      }
	    } catch (IOException e) {
	      e.printStackTrace();
	    }
	  }

	  public static void main(String[] args) throws Exception{
		OrthogroupsFunctionalAnnotator instance = new OrthogroupsFunctionalAnnotator();
	    CommandsDescriptor.getInstance().loadOptions(instance, args);
	    instance.run();
	  }
}
