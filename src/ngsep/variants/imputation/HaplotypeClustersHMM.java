package ngsep.variants.imputation;

import java.util.ArrayList;
import java.util.List;

import ngsep.hmm.RecombinationHMM;
import ngsep.math.LogMath;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;

public class HaplotypeClustersHMM extends RecombinationHMM {

	private boolean updateEmissionKnownSites = false;
	
	//Local arrays to save reallocation over many runs
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

	
	@Override
	protected void randomizeEmissions(int numState) {
		HaplotypeClusterHMMState state = (HaplotypeClusterHMMState) getState(numState);
		state.setRandomEmissions(updateEmissionKnownSites);
	}

	@Override
	protected void initEmissionsBaumWelch() {
		int k = getNumStates();
		int m = getSteps();
		if(logEmissions.length!=k || logEmissions[0].length!=m) {
			getLog().info("Creating array for emissions of dimensions "+k+" x "+m+" x 2");
			logEmissions = new Double [k][m][2];
		}
		for(int j=0;j<logEmissions.length;j++) {
			for(int i=0;i<logEmissions[j].length;i++) {
				logEmissions[j][i][0] = logEmissions[j][i][1] = null;
			}
		}
	}

	@Override
	protected void accumulateEmissionBaumWelch(int step, int stateIndex, Object datum, Double logPosterior) {
		if(datum == null || !(datum instanceof Byte)) return;
		byte allele = (Byte)datum;
		if(allele != CalledGenomicVariant.ALLELE_REFERENCE && allele!=CalledGenomicVariant.ALLELE_ALTERNATIVE) return;
		logEmissions[stateIndex][step][allele] = LogMath.logSum(logEmissions[stateIndex][step][allele], logPosterior);
	}

	@Override
	protected void updateEmissionsBaumWelch(int stateIndex) {
		HaplotypeClusterHMMState state = (HaplotypeClusterHMMState)getState(stateIndex);
		//System.out.println("Emissions state "+stateIndex+": "+logEmissions[stateIndex][0][0]+" "+logEmissions[stateIndex][0][1]);
		state.setEmissionLogProbs(logEmissions[stateIndex], updateEmissionKnownSites);
	}

	public void calculateAllelePosteriors(List<Byte> haplotype, double[][] allelePosteriors) {
		Byte b0 = 0;
		Byte b1 = 1;
		int m = haplotype.size();
		int k = getNumStates();
		Double [][] forwardLogs=calculateForward(haplotype);
		Double [][] backwardLogs=calculateBackward(haplotype);
		//calculateForward(haplotype, forwardLogs);
		//calculateBackward(haplotype, backwardLogs);
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
	
	public static List<Byte> makeHaplotypeWithHomozygous(List<CalledSNV> genotypesSample) {
		List<Byte> hapSample = new ArrayList<>();
		for(CalledSNV call:genotypesSample) {
			byte allele = CalledGenomicVariant.ALLELE_UNDECIDED;
			if(call.isHomozygousReference()) {
				allele = CalledGenomicVariant.ALLELE_REFERENCE;
			} else if (call.isHomozygous()) {
				allele = CalledGenomicVariant.ALLELE_ALTERNATIVE;
			}
			hapSample.add(allele);
		}
		return hapSample;
	}
}
