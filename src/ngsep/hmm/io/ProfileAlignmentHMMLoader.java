package ngsep.hmm.io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.hmm.ProfileAlignmentHMM;
import ngsep.hmm.ProfileAlignmentHMMState;
import ngsep.hmm.ProfileAlignmentNullModel;

public class ProfileAlignmentHMMLoader {
	private ProfileAlignmentNullModel nullModel;
	private Map<String,String> hmmDomainCodes = new HashMap<String, String>();
	
	public ProfileAlignmentHMMLoader(ProfileAlignmentNullModel nullModel) {
		super();
		this.nullModel = nullModel;
	}
	public void setDomainCode (String id, String code) {
		hmmDomainCodes.put(id, code);
	}
	public void loadDomainCodes (String filePath) throws IOException {
		try (FileReader fr = new FileReader(filePath);
			 BufferedReader reader = new BufferedReader(fr)) {
				String line = reader.readLine();
				while(line!=null) {
					String [] items = line.split("\t");
					if(items.length>=2) setDomainCode(items[0], items[1]);
					line = reader.readLine();
				}
			}
	}
	public List<ProfileAlignmentHMM> loadHMMs (String filePath) throws IOException {
		List<ProfileAlignmentHMM> hmms = new ArrayList<ProfileAlignmentHMM>();
		try (FileReader fr = new FileReader(filePath);
			 BufferedReader reader = new BufferedReader(fr)) {
			String line = reader.readLine();
			while(line!=null) {
				if(line.startsWith("HMMER")) hmms.add(loadHMM(reader));
				line = reader.readLine();
			}
		}
		return hmms;
	}
	private ProfileAlignmentHMM loadHMM(BufferedReader reader) throws IOException {
		ProfileAlignmentHMM hmm = null;
		ProfileAlignmentHMMState matchState = null;
		ProfileAlignmentHMMState insertionState = null;
		ProfileAlignmentHMMState deletionState = null;
		Double[][][] transitionMatrix = null;
		//TODO: Make dynamic alphabet
		String alphabet = "ACDEFGHIKLMNPQRSTVWY";
		String name = null;
		String id = null;
		Double miu=0.0;
		Double lambda=0.0;
		int numSteps = 0;
		List<ProfileAlignmentHMMState> states= new ArrayList<>(3);
		String line = reader.readLine();
		while (line != null && !"//".equals(line) ) {
			String [] items = line.split("\\s+");
			if (line.startsWith("NAME")) {
				name = items[1];
			} else if (line.startsWith("ACC")) {
            	id = items[1];
            } else if (line.startsWith("STATS LOCAL VITERBI")) {
				
				miu = Double.parseDouble(items[3]);
				lambda= Double.parseDouble(items[4]);
				hmm=new ProfileAlignmentHMM(id, numSteps, states,nullModel);
				hmm.setName(name);
				hmm.setDomainCode(hmmDomainCodes.get(id));
				hmm.setMiu(miu);
				hmm.setLambda(lambda);
			}
			else if (line.startsWith("LENG")) {
				numSteps = Integer.parseInt(items[1])+1;
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
			} else if (line.matches("\\s*\\d+\\s+.*")) {
				// read transitions and emmisions for each step
				String[] tokens = line.trim().split("\\s+");
				int step = Integer.parseInt(tokens[0]);
                    
				// Read match emission probabilities
				Double[] matchEmissions = new Double[alphabet.length()];
				for (int i = 1; i <= matchEmissions.length; i++) {
					matchEmissions[i - 1] = -Double.parseDouble(tokens[i]);
				}
				matchState.setStepEmissions(step, matchEmissions);
				// Read emission probabilities for the insertion state (next line)
				line = reader.readLine();
				String[] insertTokens = line.trim().split("\\s+");
				Double[] insertEmissions = new Double[alphabet.length()];
				for (int i = 0; i < insertEmissions.length; i++) {
					insertEmissions[i] = -Double.parseDouble(insertTokens[i]);
				}
				insertionState.setStepEmissions(step, insertEmissions);
                    
				// Read transition probabilities in the next line
				line = reader.readLine();
				readTransitionValues(transitionMatrix, line, step);
			}
            line = reader.readLine();
    	}
		hmm.setTransitionMatrix(transitionMatrix);
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
