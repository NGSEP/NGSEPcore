package ngsep.variants.imputation;

import java.util.Arrays;
import java.util.List;

import ngsep.hmm.RecombinationHMM;
import ngsep.math.LogMath;

public class HaplotypeClustersHMM extends RecombinationHMM {

	
	
	public static final int DEF_STARTS_BAUM_WELCH = 5;
	public static final int DEF_ITER_BAUM_WELCH = 10;
	
	private int iterationsBaumWelch = DEF_ITER_BAUM_WELCH;

	private boolean updateEmissionKnownSites = false;
	
	private List<List<Byte>> referenceHaplotypes = null;
	
	//Local arrays to save reallocation over many runs
	private Double [][] forwardLogs=new Double[0][0];
	private Double [][] backwardLogs=new Double[0][0];
	private Double [] logStarts = new Double [0];
	private Double [][][] logTransitions = new Double [0][0][0];
	private Double [][][] logEmissions = new Double [0][0][0];
	
	
	public HaplotypeClustersHMM(List<? extends HaplotypeClusterHMMState> states, int numMarkers) {
		super(states, numMarkers);
	}
	
	public HaplotypeClustersHMM(List<? extends HaplotypeClusterHMMState> states, int numMarkers, List<Integer> positions) {
		super(states, numMarkers, positions);
	}

	public boolean isUpdateEmissionKnownSites() {
		return updateEmissionKnownSites;
	}

	public void setUpdateEmissionKnownSites(boolean updateEmissionKnownSites) {
		this.updateEmissionKnownSites = updateEmissionKnownSites;
	}
	
	public List<List<Byte>> getReferenceHaplotypes() {
		return referenceHaplotypes;
	}

	public void setReferenceHaplotypes(List<List<Byte>> referenceHaplotypes) {
		this.referenceHaplotypes = referenceHaplotypes;
	}

	public void train() {
		if (referenceHaplotypes==null) throw new IllegalStateException("Training haplotypes must be provided");
		getLog().info("Training model with "+referenceHaplotypes.size()+" haplotypes");
		int n = getNumStates();
		setRandomTransitions();
		//printTransitions(0);
		//printTransitions(2);
		double logUniformStart = Math.log10(1.0/(double)n);
		for(int j=0;j<n;j++) {
			HaplotypeClusterHMMState state = (HaplotypeClusterHMMState) getState(j);
			state.setLogStart(logUniformStart);
			state.setRandomEmissions(updateEmissionKnownSites);
		}
		for(int h = 0; h < iterationsBaumWelch; h++) {
			getLog().info("Running "+h+" Baum-Welch iteration");
			runBaumWelchStep();
		}
		//printTransitions(0);
		//printTransitions(2);
	}
	
	public void printTransitions(int step) {
		int n = getNumStates();
		System.out.println("Transitions step: "+step);
		for(int j=0;j<n;j++) {
			System.out.print(getTransition(j, 0, step));
			for (int k=1;k<n;k++) {
				System.out.print("\t"+getTransition(j, k, step));
			}
			System.out.println();
		}
		
	}

	
	/**
	 * Runs a step of baum-welch training with the attribute reference haplotypes
	 */
	protected void runBaumWelchStep() {
		initArrays(getNumStates(), getSteps());
		Arrays.fill(logStarts, 0.0);
		for(int i=0;i<logTransitions.length;i++) {
			for(int j=0;j<logTransitions[i].length;j++) {
				Arrays.fill(logTransitions[i][j], null);
			}
		}
		for(int j=0;j<logEmissions.length;j++) {
			for(int i=0;i<logEmissions[j].length;i++) {
				logEmissions[j][i][0] = logEmissions[j][i][1] = null;
			}
		}
		for (List<Byte> haplotypesSample:referenceHaplotypes) {
			
			Double logProb = calculateForward(haplotypesSample, forwardLogs);
			calculateBackward(haplotypesSample, backwardLogs);
			//Calculate new starts
			for(int j=0;j<logStarts.length;j++) {
				Object o = haplotypesSample.get(0);
				Double seqProduct = LogMath.logProduct(forwardLogs[0][j], backwardLogs[0][j]);
				seqProduct = LogMath.logProduct(seqProduct, getEmission(j, o, 0));
				seqProduct = LogMath.logProduct(seqProduct, -logProb);
				logStarts[j] = LogMath.logSum(logStarts[j], seqProduct);
			}
			//Calculate new transitions
			if(!isFixedTransitions()) {
				for(int i=0;i<logTransitions.length;i++) {
					for(int j=0;j<logTransitions[i].length;j++) {
						
						for(int k=0;k<logTransitions[i][j].length;k++) {
							Object o1 = haplotypesSample.get(i);
							Object o2 = haplotypesSample.get(i+1);
							
							Double seqProduct = LogMath.logProduct(forwardLogs[i][j], backwardLogs[i+1][k]);
							seqProduct = LogMath.logProduct(seqProduct, getEmission(j, o1, i));
							seqProduct = LogMath.logProduct(seqProduct, getEmission(k, o2, i+1));
							seqProduct = LogMath.logProduct(seqProduct, getTransition(j, k, i));
							seqProduct = LogMath.logProduct(seqProduct, -logProb);
							logTransitions[i][j][k] = LogMath.logSum(logTransitions[i][j][k], seqProduct);
						}
					}
				}
			}
			
			//Calculate new emissions
			for(int j=0;j<logEmissions.length;j++) {
				for(int i=0;i<logEmissions[j].length;i++) {
					byte allele = haplotypesSample.get(i);
					if(allele==0 || allele==1) {
						Double seqProduct = LogMath.logProduct(forwardLogs[i][j], backwardLogs[i][j]);
						seqProduct = LogMath.logProduct(seqProduct, getEmission(j, allele, i));
						seqProduct = LogMath.logProduct(seqProduct, -logProb);
						logEmissions[j][i][allele] = LogMath.logSum(logEmissions[j][i][allele], seqProduct);
					}
				}
			}
			
		}
		//Normalize and update starts
		Double total = null;
		for(int j=0;j<logStarts.length;j++) total = LogMath.logSum(total, logStarts[j]);
		for(int j=0;j<logStarts.length;j++) getState(j).setLogStart(LogMath.logProduct(logStarts[j],-total));
		//Normalize and update transitions
		if(!isFixedTransitions()) {
			for(int i=0;i<logTransitions.length;i++) {
				setTransitions(logTransitions[i], i);
			}
		}
		
		//Normalize and update emissions
		for(int j=0;j<logEmissions.length;j++) {
			//TODO: Replace direct class casting with an interface
			HaplotypeClusterHMMState state = (HaplotypeClusterHMMState)getState(j);
			state.setEmissionLogProbs(logEmissions[j], updateEmissionKnownSites);
		}
	}
	
	public void calculateAllelePosteriors(List<Byte> haplotype, double[][] allelePosteriors) {
		Byte b0 = 0;
		Byte b1 = 1;
		int m = haplotype.size();
		int k = getNumStates();
		initArrays(k, m);
		calculateForward(haplotype, forwardLogs);
		calculateBackward(haplotype, backwardLogs);
		for(int i=0;i<m;i++) {
			Double log0 = null;
			Double log1 = null;
			for(int j=0;j<k;j++) {
				Double fTimesB = LogMath.logProduct(forwardLogs[i][j], backwardLogs[i][j]);
				log0 = LogMath.logSum(log0, LogMath.logProduct(fTimesB, getEmission(j, b0, i)));
				log1 = LogMath.logSum(log1, LogMath.logProduct(fTimesB, getEmission(j, b1, i)));
			}
			//Normalize and raise to calculate final probabilities of genotypes
			Double logSum = LogMath.logSum(log0, log1);
			double prob0 = LogMath.power10(LogMath.logProduct(log0, -logSum));
			double prob1 = LogMath.power10(LogMath.logProduct(log1, -logSum));
			double sum = prob0 + prob1;
			prob0/=sum;
			prob1/=sum;
			allelePosteriors[i][0] = prob0;
			allelePosteriors[i][1] = prob1;
		}
	}
	
	public void calculateStatePosteriors(List<Byte> haplotype, double[][] statePosteriors) {
		int m = haplotype.size();
		int k = getNumStates();
		initArrays(k, m);
		calculateForward(haplotype, forwardLogs);
		calculateBackward(haplotype, backwardLogs);
		Double [] stateLogPosteriors = new Double [k];
		for(int i=0;i<m;i++) {
			byte allele = haplotype.get(i);
			for(int j=0;j<k;j++) {
				Double fTimesB = LogMath.logProduct(forwardLogs[i][j], backwardLogs[i][j]);
				stateLogPosteriors[j] = LogMath.logProduct(fTimesB, getEmission(j, allele, i));
			}
			LogMath.normalizeLogs(stateLogPosteriors);
			for(int j=0;j<k;j++) {
				statePosteriors[i][j] = LogMath.power10(stateLogPosteriors[j]);
			}
		}
	}
	
	private void initArrays(int k, int m) {
		if(forwardLogs.length!=m || forwardLogs[0].length!=k) {
			getLog().info("Creating array for forward probabilities of dimensions "+m+" x "+k);
			forwardLogs = new Double[m][k];
		}
		if(backwardLogs.length!=m || backwardLogs[0].length!=k) {
			getLog().info("Creating array for backward probabilities of dimensions "+m+" x "+k);
			backwardLogs = new Double[m][k];
		}
		if(logTransitions.length!=m-1 || logTransitions[0].length!=k) {
			getLog().info("Creating array for transitions of dimensions "+(m-1)+" x "+k+" x "+k);
			logTransitions = new Double [m-1][k][k];
		}
		if(logEmissions.length!=k || logEmissions[0].length!=m) {
			getLog().info("Creating array for emissions of dimensions "+k+" x "+m+" x 2");
			logEmissions = new Double [k][m][2];
		}
		if(logStarts.length!=k) logStarts = new Double [k];
	}

}
