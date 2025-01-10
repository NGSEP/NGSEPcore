package ngsep.hmm.io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ngsep.hmm.ProfileAlignmentHMM;
import ngsep.hmm.ProfileAlignmentHMMState;
import ngsep.hmm.ProfileAlignmentNullModel;

public class ProfileAlignmentHMMLoader {
	
	public static ProfileAlignmentHMM loadHMM(String filePath, ProfileAlignmentNullModel nullModel) throws IOException {
		ProfileAlignmentHMM hmm = null;
		ProfileAlignmentHMMState matchState = null;
		ProfileAlignmentHMMState insertionState = null;
		ProfileAlignmentHMMState deletionState = null;
		Double[][][] transitionMatrix = null;
		//TODO: Make dynamic alphabet
		String alphabet = "ACDEFGHIKLMNPQRSTVWY";
		try (FileReader fr = new FileReader(filePath);
			 BufferedReader reader = new BufferedReader(fr)) {
			String line;
			String name = null;
			Double miu=0.0;
			Double lambda=0.0;
			int numSteps = 0;
			List<ProfileAlignmentHMMState> states= new ArrayList<>(3);
			
			while ((line = reader.readLine()) != null) {
				if (line.startsWith("NAME")) {
					name =line.split("\\s+")[1];
                }
            	
            	else if (line.startsWith("NSEQ")) {
            		//nseq = Integer.parseInt(line.split("\\s+")[1]);
            		//System.out.println("nseq: "+nseq);
            	}
            	
            	else if (line.startsWith("STATS LOCAL VITERBI")) {
            		miu = Double.parseDouble(line.split("\\s+")[3]);
            		lambda= Double.parseDouble(line.split("\\s+")[4]);
//            		System.out.println("miu: "+miu);
//            		System.out.println("lambda: "+lambda);
            		hmm=new ProfileAlignmentHMM(name, numSteps, states,nullModel);
            		hmm.setMiu(miu);
            		hmm.setLambda(lambda);
            	}
            	else if (line.startsWith("LENG")) {
                    // Extraer el número de steps
                    numSteps = Integer.parseInt(line.split("\\s+")[1])+1;
                    matchState=new ProfileAlignmentHMMState("match", alphabet, numSteps);
                    insertionState=new ProfileAlignmentHMMState("insetion", alphabet, numSteps);
                    deletionState=new ProfileAlignmentHMMState("deletion", alphabet, numSteps);
                    
                    states.add(matchState);
                    states.add(insertionState);
                    states.add(deletionState);
                    
                    transitionMatrix= new Double [numSteps][3][3];



                }
            	if (line.startsWith("  COMPO")) {
                	line=reader.readLine();
                	int step = 0;
                	String[] insertTokens = line.trim().split("\\s+");
                    Double[] insertEmissions = new Double[alphabet.length()];
                	for (int i = 0; i < insertEmissions.length; i++) {
                        insertEmissions[i] = -Double.parseDouble(insertTokens[i]);
                    }
                    insertionState.setStepEmissions(step, insertEmissions);
                    
                    line=reader.readLine();
                    readTransitionValues(transitionMatrix, line, step);
                }
                else if (line.matches("\\s*\\d+\\s+.*")) {
                    // Aquí leeríamos las probabilidades de emisión y transición por cada step
                    // Implementar el procesamiento de la línea para extraer las probabilidades.
                    String[] tokens = line.trim().split("\\s+");
                    int step = Integer.parseInt(tokens[0]);
                    
                    // Leer las probabilidades de emisión para el estado match
                    Double[] matchEmissions = new Double[20];
                    for (int i = 1; i <= 20; i++) {
                        matchEmissions[i - 1] = -Double.parseDouble(tokens[i]);
                    }
                    //deberia ser de la clase HMMState pero el metodo no esta definido, prefiero no tocar la interfaz
                    matchState.setStepEmissions(step, matchEmissions);
                    // Leer las probabilidades de emisión para el estado de inserción (en la siguiente línea)
                    line = reader.readLine();
                    String[] insertTokens = line.trim().split("\\s+");
                    Double[] insertEmissions = new Double[20];
                    for (int i = 0; i < 20; i++) {
                        insertEmissions[i] = -Double.parseDouble(insertTokens[i]);
                    }
                    insertionState.setStepEmissions(step, insertEmissions);
                    
                    // Read transition probabilities in the next line
                    line = reader.readLine();
                    readTransitionValues(transitionMatrix, line, step);
                }
            }
            hmm.setTransitionMatrix(transitionMatrix);
    	}
        return hmm;
    }

	private static void readTransitionValues(Double[][][] transitionMatrix, String line, int step) {
		String[] transitionTokens = line.trim().split("\\s+");
		// Assign transition values : m->m, m->i, m->d, etc.
		transitionMatrix[step][0][0] = transitionTokens[0].equals("*") ? null : -Double.parseDouble(transitionTokens[0]);  // m -> m
		transitionMatrix[step][0][1] = transitionTokens[1].equals("*") ? null : -Double.parseDouble(transitionTokens[1]);  // m -> i
		transitionMatrix[step][0][2] = transitionTokens[2].equals("*") ? null : -Double.parseDouble(transitionTokens[2]);  // m -> d

		transitionMatrix[step][1][0] =transitionTokens[3].equals("*") ? null : -Double.parseDouble(transitionTokens[3]); // i -> m
		transitionMatrix[step][1][1] = transitionTokens[4].equals("*") ? null : -Double.parseDouble(transitionTokens[4]);// i -> i
		transitionMatrix[step][1][2] = null;  // i -> d not possible

		transitionMatrix[step][2][0] = transitionTokens[5].equals("*") ? null : -Double.parseDouble(transitionTokens[5]);  // d -> m
		transitionMatrix[step][2][1] = null;  // d -> i not possible
		transitionMatrix[step][2][2] = transitionTokens[6].equals("*") ? null : -Double.parseDouble(transitionTokens[6]); // d -> d
	}
}
