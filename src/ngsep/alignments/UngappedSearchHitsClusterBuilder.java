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
package ngsep.alignments;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import ngsep.math.CollisionEntropyCalculator;
import ngsep.math.Distribution;
import ngsep.math.EntropyCalculator;
import ngsep.sequences.UngappedSearchHit;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class UngappedSearchHitsClusterBuilder {
	public static final int CLUSTERING_ALGORITHM_KRUSKAL_LIKE = 0;
	public static final int CLUSTERING_ALGORITHM_KMEANS_LIKE = 1;
	private int idxSubjectDebug = -1;
	private int queryLengthDebug = -1;
	private boolean debug = false;
	private int clusteringAlgorithm = CLUSTERING_ALGORITHM_KRUSKAL_LIKE;
	private EntropyCalculator entropyCalculator = new CollisionEntropyCalculator(4);
	
	public int getClusteringAlgorithm() {
		return clusteringAlgorithm;
	}
	public void setClusteringAlgorithm(int clusteringAlgorithm) {
		this.clusteringAlgorithm = clusteringAlgorithm;
	}
	public List<UngappedSearchHitsCluster> clusterRegionKmerAlns(int queryLength, int subjectIdx, int subjectLength, List<UngappedSearchHit> sequenceHits) {
		debug = subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug;
		double minHits = Math.min(20,0.01*queryLength);
		if(sequenceHits.size()<minHits) return new ArrayList<>();
		UngappedSearchHitClusteringAlgorithm alg;
		if(clusteringAlgorithm == CLUSTERING_ALGORITHM_KRUSKAL_LIKE) alg = new UngappedSearchHitClusteringAlgorithmKruskal();
		else alg = new UngappedSearchHitClusteringAlgorithmKmeansLike();
		List<List<UngappedSearchHit>> hitsClusters = alg.clusterLocalSearchHits(sequenceHits);
		
		//Select hits and within clusters and build objects
		Collections.sort(hitsClusters,(l1,l2)->l2.size()-l1.size());
		
		List<UngappedSearchHitsCluster> answer = new ArrayList<>();
		for(List<UngappedSearchHit> hits:hitsClusters) {
			if(debug) System.out.println("Next candidate cluster size: "+hits.size()+" current limit: "+minHits);
			if(hits.size()<minHits) break;
			List<List<UngappedSearchHit>> subclusters = breakByQueryStart(hits);
			if(subclusters.size()>1) System.err.println("WARN. Cluster broken in subclusters by query starts. SubjectIdx: "+subjectIdx+" query length: "+queryLength);
			for(List<UngappedSearchHit> subcluster:subclusters) {
				if(subcluster.size()<minHits) continue;
				//List<UngappedSearchHit> selectedHits = collapseAndSelectSortedHits(queryLength, subjectIdx, subjectLength, subcluster);
				List<UngappedSearchHit> selectedHits = selectHits(queryLength, subjectIdx, subjectLength, subcluster);
				if(selectedHits.size()<minHits) continue;
				UngappedSearchHitsCluster cluster = new UngappedSearchHitsCluster(queryLength, subjectIdx, subjectLength, selectedHits);
				cluster.setRawKmerHits(sequenceHits.size());
				//TODO: Calculate better or remove
				cluster.setRawKmerHitsSubjectStartSD(1);
				if(debug) System.out.println("Next cluster subject predicted coords: "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" subject evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" query evidence: "+cluster.getQueryEvidenceStart()+" "+cluster.getQueryEvidenceEnd()+" number of hits: "+cluster.getCountKmerHitsCluster());
				if(cluster.getCountKmerHitsCluster()>=minHits) answer.add(cluster);
				minHits = Math.max(minHits, cluster.getCountKmerHitsCluster()/2);
			}
		}
		return answer;
	}
	
	
	/**
	 * Breaks a hits cluster based on large differences in query start. This should not be needed if a good initial clustering algorithm is used
	 * @param hits
	 * @return List<List<UngappedSearchHit>> subclusters
	 */
	private List<List<UngappedSearchHit>> breakByQueryStart(List<UngappedSearchHit> hits) {
		List<List<UngappedSearchHit>> answer = new ArrayList<List<UngappedSearchHit>>();
		Collections.sort(hits, (h1,h2)->h1.getQueryStart()-h2.getQueryStart());
		List<UngappedSearchHit> nextSubcluster = new ArrayList<UngappedSearchHit>();
		int lastQueryStart = 0;
		for(UngappedSearchHit hit:hits) {
			//TODO: define better this parameter
			if(hit.getQueryStart()>lastQueryStart+30000) {
				if(nextSubcluster.size()>0) answer.add(nextSubcluster); 
				nextSubcluster = new ArrayList<UngappedSearchHit>();
			}
			nextSubcluster.add(hit);
			lastQueryStart = hit.getQueryStart();
		}
		if(nextSubcluster.size()>0) answer.add(nextSubcluster);
		return answer;
	}
	
	// Metodo para seleccionar los hits a usar en el alineamiento de una lectura. Se seleccionan haciendo uso de un grafo de
	// compatibilidad en el cual se busca el camino mas pesado
	private List<UngappedSearchHit> selectHits(int queryLength, int subjectIdx, int subjectLength, List<UngappedSearchHit> hits) {

		// Se ordenan los hits en funcion de su posicion de inicio en el query
		Collections.sort(hits, (h1, h2) -> Integer.compare(h1.getQueryStart(), h2.getQueryStart()));

		// Estructuras para manejar la informacion del camino más pesado
		double[] maxWeight = new double[hits.size()];
		int[] predecessor = new int[hits.size()];

		// Llenado de las estructuras a los valores por defecto 
		Arrays.fill(maxWeight, Double.NEGATIVE_INFINITY);
		Arrays.fill(predecessor, -1);

		// Se inicializa la minima suma de distancias en infinito positivo
		double  minimumSum = Double.POSITIVE_INFINITY;

		// Se inicializa el peso del primer hit en cero para que el algoritmo tenga un punto de partida
		maxWeight[0] = 0;

		// Se recorre el grafo siguiendo el orden definido por el inicio del hit en el query
		for(int sourceHit = 0; sourceHit < hits.size(); sourceHit++) {
			// Unicamente se recorren los hits que estan a la derecha del hit seleccionado
			for(int destinationHit = sourceHit + 1; destinationHit < hits.size(); destinationHit++) {
				// Existe un eje entre dos hits si el destination hit esta estrictamente a la derecha del source hit en el query y
				// si la distancia no excede el limite planteado por la minima distancia
				int currentSum = calculateDistancesSum(hits.get(sourceHit), hits.get(destinationHit));
				// Si la suma actual es 10 veces mayor (o igual) a la minima suma la ignoro y hago break, cerrando el ciclo interno
				// De paso, tambien se actualiza la minima suma
				if(currentSum < minimumSum) minimumSum = currentSum;
				else if(minimumSum != Double.POSITIVE_INFINITY && currentSum >= 10 * minimumSum) break;
				// Se verifica si se cumplen las condicones para que exista un eje entre los hits
				if(checkIfExistEdge(hits.get(sourceHit), hits.get(destinationHit))) {
					// Se calcula el nuevo peso teniendo en cuenta el puntaje escogido
					double newWeight = maxWeight[sourceHit] + calculateScoreWithoutPenalty(hits.get(sourceHit), hits.get(destinationHit));
					// Se lleva a cabo el proceso de relajacion de los pesos
					if(newWeight > maxWeight[destinationHit]){
						maxWeight[destinationHit] = newWeight;
						predecessor[destinationHit] = sourceHit;
					}
				}
			}
		}

		// Se obtiene el hit en el cual termina el mejor camino
		double finalMaxWeight = Double.NEGATIVE_INFINITY;
		int finalHitIndex = -1;
		// Para encontra el hit en el cual termina el mejor camino se itera
		// sobre todos los hits y se actualiza las variables de informacion
		for (int i = 0; i < maxWeight.length; i++) {
    		if (maxWeight[i] > finalMaxWeight) {
        		finalMaxWeight = maxWeight[i];
        		finalHitIndex = i;
    		}
		}

		// Estructuras para almacenar el camino que logra el mayor peso
		List<UngappedSearchHit> path = new ArrayList<>();
    	int actualHitIndex = finalHitIndex;

		// Se reconstruye el camino mas pesado usando bracktraking
    	while (actualHitIndex != -1) {
        	path.add(hits.get(actualHitIndex));
        	actualHitIndex = predecessor[actualHitIndex];
    	}

		// Se invierte el orden del camino para que sea correcto
    	Collections.reverse(path);
    	return path;
	}

	// Funcion para verificar si existe un eje entre dos hits, el destination hit debe estar estrictamente a la derecha
	// del subjet hit en el query 
	private boolean checkIfExistEdge(UngappedSearchHit sourceHit, UngappedSearchHit destinationHit) {
		return (destinationHit.getQueryStart() - sourceHit.getQueryStart()) > 0;
	}

	private int calculateDistancesSum(UngappedSearchHit sourceHit, UngappedSearchHit destinationHit) {
		// Para la distancia en el query no se utiliza el valor absoluto ya que por precondicion se garantiza que esta distancia
		// siempre es positiva. Para la distancia en el subject si se utiliza el valor absoluto ya que pueden existir inversiones
		int queryDistance = destinationHit.getQueryStart() - sourceHit.getQueryStart();		
		int subjectDistance = Math.abs(destinationHit.getSubjectStart() - sourceHit.getSubjectStart());

		return queryDistance + subjectDistance;
	}

	// Funcion para calcular el peso de un eje, teniendo en cuenta la entropia de los hits y la distancia entre los mismos 
	// sin penalidad
	private double calculateScoreWithoutPenalty(UngappedSearchHit sourceHit, UngappedSearchHit destinationHit) {
		// Extraccion de las entropias de los hits
		double sourceHitEntropy = entropyCalculator.denormalizeEntropy(sourceHit.getWeight());
		double destinationHitEntropy = entropyCalculator.denormalizeEntropy(destinationHit.getWeight());

		// Se calcula la suma de las distancias en el subject y en el query
		int distancesSum = calculateDistancesSum(sourceHit, destinationHit);

		// Se  calcula el puntaje de los hits, se multiplica el numerador por 100 para evitar problemas de precision asociado a los 
		// decimales, a mayor puntaje mejor calidad de hit
		return (100 * sourceHitEntropy * destinationHitEntropy) / distancesSum;
	}

	// Funcion para calcular el peso de un eje, teniendo en cuenta la entropia de los hits y la distancia entre los mismos.
	// Se aplica un puntaje de penalizacion usando la consistencia de las distancias
	private double calculateScoreWithPenalty(UngappedSearchHit sourceHit, UngappedSearchHit destinationHit) {
		// Extraccion de las entropias de los hits
		double sourceHitEntropy = entropyCalculator.denormalizeEntropy(sourceHit.getWeight());
		double destinationHitEntropy = entropyCalculator.denormalizeEntropy(destinationHit.getWeight());

		// Para la distancia en el query no se utiliza el valor absoluto ya que por precondicion se garantiza que esta distancia
		// siempre es positiva. Para la distancia en el subject si se utiliza el valor absoluto ya que pueden existir inversiones
		int queryDistance = destinationHit.getQueryStart() - sourceHit.getQueryStart();
		int subjectDistance = Math.abs(destinationHit.getSubjectStart() - sourceHit.getSubjectStart());

		// Se calcula la penalizacion asociada a la consistencia de las distancias
		double distancePenalty = 1.0 + Math.pow(Math.abs(queryDistance - subjectDistance), 2); 

		// Se  calcula el puntaje de los hits, se multiplica el numerador por 100 para evitar problemas de precision asociado a los 
		// decimales, a mayor puntaje mejor calidad de hit
		return (100 * sourceHitEntropy * destinationHitEntropy * distancePenalty) / (queryDistance + subjectDistance);
	}

	// Funcion para calcular el peso de un eje, teniendo en cuenta la distancia entre los hits sin penalidad
	private double calculateScoreWithoutEntropy(UngappedSearchHit sourceHit, UngappedSearchHit destinationHit) {
		// Se calcula la suma de las distancias en el subject y en el query
		int distancesSum = calculateDistancesSum(sourceHit, destinationHit);
		
		// Se  calcula el puntaje de los hits, se multiplica el numerador por 100 para evitar problemas de precision asociado a los 
		// decimales, a mayor puntaje mejor calidad de hit
		return (1d / distancesSum);
	}
	
	/**
	 * 
	 * @param queryLength
	 * @param subjectLength
	 * @param hits sorted by queryidx
	 * @return
	 */
	private List<UngappedSearchHit> collapseAndSelectSortedHits(int queryLength, int subjectIdx, int subjectLength, List<UngappedSearchHit> inputHits) {
		
		if(debug) System.out.println("KmerHitsCluster. Clustering "+inputHits.size()+" hits. Subject idx: "+subjectIdx);
		
		
		double [] stats = calculateStatsEstimatedSubjectStart(inputHits);
		
		int median = (int) Math.round(stats[1]);
		//System.out.println("Sum: "+sum+" sum2: "+sum2);
		double variance = stats[2];
		double rawKmerHitsSubjectStartSD = Math.max(1, stats[3]);
		Distribution dist = new Distribution(-500, 500, 50);
		Distribution distAbs = new Distribution(0, 500, 1);
		for(UngappedSearchHit hit:inputHits) {
			int start = hit.estimateSubjectStart();
			int distance = start-median;
			int distanceAbs = Math.abs(distance);
			if (distanceAbs < 2*rawKmerHitsSubjectStartSD) {
				dist.processDatapoint(distance);
				distAbs.processDatapoint(distanceAbs);
			}
		}
		int modeDist = (int) Math.round(dist.getLocalMode(-400, 400));
		if(Math.abs(modeDist)>=100) median+=(modeDist/2);
		
		//System.out.println(subjectStarts);
		if(debug) dist.printDistributionInt(System.out);
		int maxDistance = (int) Math.max(distAbs.getAverage(), Math.sqrt(distAbs.getVariance()+1));
		maxDistance *=5;
		maxDistance=Math.max(queryLength/40,maxDistance);
		maxDistance=Math.min(300,maxDistance);
		maxDistance=Math.max(10,maxDistance);
		//if(maxDistance<0.01*query.length()) maxDistance*=2;
		if(debug) System.out.println("KmerHitsCluster. Num hits: "+inputHits.size()+" median: "+median+" average: "+stats[0]+" variance: "+variance+" stdev: "+rawKmerHitsSubjectStartSD+" abs distance avg: "+distAbs.getAverage()+" stdev "+Math.sqrt(distAbs.getVariance())+" max distance: "+maxDistance);
		//Index hits by query kmer start
		Map<Integer, List<UngappedSearchHit>> hitsMultiMap = indexByQueryStart(inputHits);
		if(debug) System.out.println("KmerHitsCluster. Num different kmers: "+hitsMultiMap.size());
		List<UngappedSearchHit> selectedHits = selectHitsByDistanceWithMedian(hitsMultiMap, median, maxDistance);
		if(selectedHits.size()<1) { 
			if(debug) System.err.println("WARN. Empty list of selected hits for subject: "+subjectIdx);
			return selectedHits;
		}
		
		List<UngappedSearchHit> filteredHits = new ArrayList<UngappedSearchHit>();
		/*Set<Integer> outliersToRemove =*/ replaceHitsByLocalAgreement(selectedHits, hitsMultiMap, median, maxDistance, queryLength);
		//if(outliersToRemove.size()>0 && outliersToRemove.size()<0.2*selectedHits.size()) {
			//for(int i=0;i<selectedHits.size();i++) {
				//if(!outliersToRemove.contains(i)) filteredHits.add(selectedHits.get(i));
				//else if (debug) System.out.println("Removing outlier at qpos "+selectedHits.get(i).getQueryStart()+" hit: "+selectedHits.get(i).getSubjectStart());
				
			//}
		//} else filteredHits.addAll(selectedHits);
		filteredHits.addAll(selectedHits);
		filteredHits = removeDisorganized (filteredHits, median);
		if(filteredHits.size()<1) {
			if(debug) System.err.println("WARN. Empty list of sorted hits for subject: "+subjectIdx+" selected hits: "+selectedHits.size()+" query length: "+queryLength);
			
		}
		return filteredHits;
	}
	public static double [] calculateStatsEstimatedSubjectStart(List<UngappedSearchHit> clusterHits) {
		double sum = 0;
		int n=0;
		List<Integer> starts = new ArrayList<>();
		for(UngappedSearchHit hit:clusterHits) {
			int start = hit.estimateSubjectStart();
			starts.add(start);
			sum+=start;
			n++;
		}
		Collections.sort(starts);
		double [] answer = {0,0,1,1,0};
		if(n>0) {
			answer[0] = sum/n;
			answer[1] = starts.get(n/2);
			if(n%2==0) {
				answer[1] += starts.get(n/2-1);
				answer[1] /=2;
			}
		}
		if(n>1) {
			//Better than sum squares because numbers are big
			double sumd2 = 0;
			Distribution dist = new Distribution(-500, 500, 50);
			for(int start:starts) {
				double d = start-answer[0];
				dist.processDatapoint(d);
				sumd2+=d*d;
			}
			answer[2] = sumd2/(n-1);
			//if(answer[2]<0) answer[2]=0;
			answer[3] = Math.sqrt(answer[2]);
			int modeDist = (int) Math.round(dist.getLocalMode(-500, 500));
			answer[4] = modeDist;
		}
		return answer;
	}
	
	private Map<Integer, List<UngappedSearchHit>> indexByQueryStart(List<UngappedSearchHit> inputHits) {
		Map<Integer,List<UngappedSearchHit>> hitsMultiMap = new TreeMap<Integer, List<UngappedSearchHit>>();
		for(UngappedSearchHit hit:inputHits) {
			//if (subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" kmer: "+hit.getQuery()+ " estq: "+estimateQueryStart(hit)+" - "+estimateQueryEnd(hit)+" estS: "+estimateSubjectStart(hit)+" - "+estimateSubjectEnd(hit));
			List<UngappedSearchHit> list = hitsMultiMap.computeIfAbsent(hit.getQueryStart(), l -> new ArrayList<UngappedSearchHit>());
			list.add(hit);
		}
		return hitsMultiMap;
	}
	private List<UngappedSearchHit> selectHitsByDistanceWithMedian(Map<Integer, List<UngappedSearchHit>> hitsMultiMap, int median, int maxDistance) {
		List<UngappedSearchHit> selectedHits = new ArrayList<UngappedSearchHit>(hitsMultiMap.size());
		for(List<UngappedSearchHit> hits:hitsMultiMap.values()) {
			UngappedSearchHit hit = selectHit(hits, median, maxDistance);
			if(hit!=null) {
				if (debug) System.out.println("Selected hits. Next qpos "+hit.getQueryStart()+" hit: "+hit.getSubjectStart()+" estq: "+hit.estimateQueryStart()+" estS: "+hit.estimateSubjectStart()+" all starts: "+calculateHitStarts(hits));
				selectedHits.add(hit);
			}
		}
		return selectedHits;
	}

	private List<Integer> calculateHitStarts(List<UngappedSearchHit> hits) {
		List<Integer> hitStarts = new ArrayList<Integer>();
		for(UngappedSearchHit hit2:hits) hitStarts.add(hit2.estimateSubjectStart());
		return hitStarts;
	}
	
	private UngappedSearchHit selectHit(List<UngappedSearchHit> hits, int median, int maxDistance) {
		UngappedSearchHit answer = null;
		int minCost = 0;
		
		for(UngappedSearchHit hit:hits) {
			int estStart = hit.estimateSubjectStart();
			int distance = Math.abs(estStart-median);
			int cost = distance;
			//if(debug && hit.getQueryIdx()==5314) System.out.println("Selecting hit. Next hit: "+hit.getStart()+" distance: "+distance+" max "+maxDistance+" cost: "+cost);
			if(distance <= maxDistance) {
				if(answer == null || minCost>cost) {
					answer = hit;
					minCost = cost;
				}
			}
		}
		return answer;
	}
	
	private Set<Integer> replaceHitsByLocalAgreement(List<UngappedSearchHit> selectedHits, Map<Integer, List<UngappedSearchHit>> hitsMultiMap, int median, int maxDistance, int queryLength) {
		//Find trustable site
		int minHitPos = -1;
		int minCost = -1;
		int midPoint = queryLength/2;
		for(int i=0;i<selectedHits.size();i++) {
			UngappedSearchHit hit = selectedHits.get(i);
			int estStart = hit.estimateSubjectStart();
			int distance = Math.abs(estStart-median);
			int d2 = Math.abs(hit.getQueryStart()-midPoint);
			int cost = distance + d2/10;
			if(minHitPos==-1 || cost<minCost) {
				minHitPos = i;
				minCost = cost;
			}
		}
		if(minHitPos==-1) return new HashSet<Integer>();
		UngappedSearchHit minHit = selectedHits.get(minHitPos);
		if (debug) System.out.println("Hit closest to median. qpos "+minHit.getQueryStart()+" hit: "+minHit.getSubjectStart()+" cost: "+minCost+" estq: "+minHit.estimateQueryStart()+" estS: "+minHit.estimateSubjectStart());
		int vicinityEstStart = minHit.estimateSubjectStart();
		LinkedList<Integer> localValues = new LinkedList<Integer>();
		localValues.add(vicinityEstStart);
		Set<Integer> posHitsToRemove = new HashSet<Integer>();
		for(int i=minHitPos-1;i>=0;i--) {
			UngappedSearchHit nextHit = selectedHits.get(i);
			int estStart = nextHit.estimateSubjectStart();
			if(Math.abs(estStart-vicinityEstStart)> 5) {
				List<UngappedSearchHit> hits = hitsMultiMap.get(nextHit.getQueryStart());
				vicinityEstStart = getMedian(localValues);
				UngappedSearchHit updatedHit = selectHit(hits, vicinityEstStart, 10);
				if (debug) System.out.println("Trying to replace outlier at qpos "+nextHit.getQueryStart()+" hit: "+nextHit.getSubjectStart()+" ests: "+estStart+" local median: "+vicinityEstStart+" all starts: "+calculateHitStarts(hits));
				if(updatedHit!=null && updatedHit!=nextHit) {
					if (debug) System.out.println("Replacing hit. qpos "+updatedHit.getQueryStart()+" hit: "+updatedHit.getSubjectStart()+" estq: "+updatedHit.estimateQueryStart()+" estS: "+updatedHit.estimateSubjectStart());
					selectedHits.set(i, updatedHit);
					estStart = updatedHit.estimateSubjectStart();
				} else {
					posHitsToRemove.add(i);
				}
			} else {
				vicinityEstStart = estStart;
			}
			localValues.add(estStart);
			if(localValues.size()>10) localValues.removeFirst();
		}
		vicinityEstStart = minHit.estimateSubjectStart();
		localValues.clear();
		localValues.add(vicinityEstStart);
		for(int i=minHitPos+1;i<selectedHits.size();i++) {
			UngappedSearchHit nextHit = selectedHits.get(i);
			int estStart = nextHit.estimateSubjectStart();
			if(Math.abs(estStart-vicinityEstStart)> 5) {
				List<UngappedSearchHit> hits = hitsMultiMap.get(nextHit.getQueryStart());
				vicinityEstStart = getMedian(localValues);
				UngappedSearchHit updatedHit = selectHit(hits, vicinityEstStart, 10);
				if (debug) System.out.println("Trying to replace outlier at qpos "+nextHit.getQueryStart()+" hit: "+nextHit.getSubjectStart()+" ests: "+estStart+" local median: "+vicinityEstStart+" all starts: "+calculateHitStarts(hits));
				if(updatedHit!=null && updatedHit!=nextHit) {
					if (debug) System.out.println("Replacing hit. qpos "+updatedHit.getQueryStart()+" hit: "+updatedHit.getSubjectStart()+" estq: "+updatedHit.estimateQueryStart()+" estS: "+updatedHit.estimateSubjectStart());
					selectedHits.set(i, updatedHit);
					estStart = updatedHit.estimateSubjectStart();
				} else {
					posHitsToRemove.add(i);
				}
			} else {
				vicinityEstStart = estStart;
			}
			localValues.add(estStart);
			if(localValues.size()>10) localValues.removeFirst();
		}
		return posHitsToRemove;
	}

	private int getMedian(LinkedList<Integer> localValues) {
		ArrayList<Integer> copy = new ArrayList<Integer>(localValues);
		Collections.sort(copy);
		return copy.get(copy.size()/2);
	}

	private List<UngappedSearchHit> removeDisorganized(List<UngappedSearchHit> selectedHits, int median) {
		int n = selectedHits.size();
		//if(n<=30) return selectSorted(selectedHits, 0, n-1, median);
		boolean [] vicinityConsistent = new boolean [n];
		Arrays.fill(vicinityConsistent, false);
		for(int i=10;i<n-10;i++) {
			UngappedSearchHit hit = selectedHits.get(i);
			int start = hit.getSubjectStart();
			boolean consistent = true;
			for(int j=i-1;j>=i-10 && consistent;j--) {
				UngappedSearchHit hit2 = selectedHits.get(j);
				if(hit2.getSubjectStart()>=start) consistent = false; 
			}
			for(int j=i+1;j<=i+10 && consistent;j++) {
				UngappedSearchHit hit2 = selectedHits.get(j);
				if(hit2.getSubjectStart()<=start) consistent = false; 
			}
			vicinityConsistent[i] = consistent;
			if(!consistent) {
				for(int j=i-1;j>=i-10;j--) vicinityConsistent[j] = false;
				for(int j=i+1;j<=i+10 ;j++) vicinityConsistent[j] = false;
				i+=10;
			}
		}
		List<UngappedSearchHit> dpSortedHits = new ArrayList<UngappedSearchHit>();
		double sum = 0;
		double sum2 = 0;
		double nS = 0;
		int i=0;
		while(i<n) {
			UngappedSearchHit hit = selectedHits.get(i);
			if(vicinityConsistent[i]) {
				//if(debug) System.out.println("Adding consistent hit at "+hit.getQueryIdx());
				dpSortedHits.add(hit);
				double s = hit.estimateSubjectStart();
				sum+=s;
				sum2+=s*s;
				nS++;
				i++;
			} else {
				int j;
				for(j=i+1;j<n && !vicinityConsistent[j];j++);
				dpSortedHits.addAll(selectSorted(selectedHits,i,j-1,median));
				i=j;
			}
		}
		//Remove final outliers
		int n2= dpSortedHits.size();
		if(debug) System.out.println("Initial hits: "+n+" consistent: "+nS+" dpSorted: "+n2);
		if(n2<30 || nS < 20) return dpSortedHits;
		int averageStart = (int) (sum/nS);
		double variance = Math.max(1, (sum2-sum*sum/nS)/(nS-1));
		double stdev = Math.sqrt(variance);
		double maxDistance = Math.max(100, 3*stdev);
		int countOutliers =0;
		for(UngappedSearchHit hit:dpSortedHits) {
			double distance = Math.abs(hit.estimateSubjectStart()-averageStart);
			if(distance > maxDistance) countOutliers++;
		}
		
		if(debug) System.out.println("Average start: "+averageStart+" stdev: "+stdev+" max distance: "+maxDistance+" outliers: "+countOutliers+" n "+n2+" nconsistent: "+nS);
		List<UngappedSearchHit> answer = new ArrayList<UngappedSearchHit>(n2);
		if(countOutliers<0.02*n2) {
			for(UngappedSearchHit hit:dpSortedHits) {
				double distance = Math.abs(hit.estimateSubjectStart()-averageStart);
				if(distance <= maxDistance) answer.add(hit);
				else if(debug) System.out.println("Removing outlier: "+hit.getQueryStart()+" subject: "+hit.getSubjectStart()+" estimated start: "+hit.estimateSubjectStart());
			}
			return answer;
		}
		Collections.sort(dpSortedHits,(h1,h2)->Math.abs(h2.estimateSubjectStart()-averageStart)-Math.abs(h1.estimateSubjectStart()-averageStart));
		//if(debug) System.out.println("Sorted by distance");
		for (int k=0;k<dpSortedHits.size();k++) {
			UngappedSearchHit hit = dpSortedHits.get(k);
			double distance = Math.abs(hit.estimateSubjectStart()-averageStart);
			if(distance <= maxDistance || k>=0.02*n2) answer.add(hit);
			else if(debug) System.out.println("Removing outlier after sorting: "+hit.getQueryStart()+" subject: "+hit.getSubjectStart()+" estimated start: "+hit.estimateSubjectStart());
		}
		Collections.sort(answer,(h1,h2)->h1.getQueryStart()-h2.getQueryStart());
		return answer;
	}


	private List<UngappedSearchHit> selectSorted(List<UngappedSearchHit> selectedHits, int first, int last, int median) {
		if(debug) System.out.println("DP Sorting hits from "+first+" to "+last+" median start: "+median);
		//Check sorted first
		List<UngappedSearchHit> regionHits = new ArrayList<UngappedSearchHit>(last-first+1);
		int lastStartSubject = -1;
		int lastStartQuery = -1;
		boolean sortedSubject = true;
		boolean sortedQuery = true;
		for(int i=first;i<=last;i++) {
			UngappedSearchHit hit = selectedHits.get(i);
			regionHits.add(hit);
			int startS = hit.getSubjectStart();
			if(lastStartSubject>=startS) sortedSubject = false;
			lastStartSubject = startS;
			int startQ = hit.getQueryStart();
			if(lastStartQuery>=startQ) sortedQuery = false;
			lastStartQuery = startQ;
		}
		if(!sortedQuery) {
			System.err.println("WARN. Disorganized query hits list");
			Collections.sort(regionHits,(h1,h2)->h1.getQueryStart()-h2.getQueryStart());
			
		} else if (sortedSubject) return regionHits;
		List<UngappedSearchHit> regionHitsBySubject = new ArrayList<UngappedSearchHit>();
		regionHitsBySubject.addAll(regionHits);
		Collections.sort(regionHitsBySubject,(h1,h2)->h1.getSubjectStart()-h2.getSubjectStart());
		int n = regionHits.size();
		//Dynamic programming to select the best organized subset
		int [][] scores = new int [n+1][n+1];
		for(int i=0;i<=n;i++) {
			for(int j=0;j<=n;j++) {
				if(i==0 || j==0) scores[i][j] = 0;
				else {
					UngappedSearchHit hit1 = regionHits.get(i-1);
					UngappedSearchHit hit2 = regionHitsBySubject.get(j-1);
					if(hit1==hit2) scores[i][j] = calculateScore(hit1,median) + scores[i-1][j-1];
					else scores[i][j] = Math.max(scores[i-1][j], scores[i][j-1]);
				}
			}
			 
		}
		//if(debug) printMatrix(scores);
		List<UngappedSearchHit> answer = new ArrayList<UngappedSearchHit>();
		int i= n;
		int j = n;
		while(i>0 && j>0) {
			UngappedSearchHit hit1 = regionHits.get(i-1);
			UngappedSearchHit hit2 = regionHitsBySubject.get(j-1);
			if(hit1==hit2) {
				//if(debug) System.out.println("Added DP sorted hit at "+hit1.getQueryIdx());
				answer.add(hit1);
				i--;
				j--;
			} else if (scores[i][j] == scores[i-1][j]) i--;
			else j--;
		}
		Collections.reverse(answer);
		return answer;
	}
	private int calculateScore(UngappedSearchHit hit, int median) {
		/*int start = estimateSubjectStart(hit);
		int distance = Math.abs(start-median);
		int score = (int) (1000*Math.max(0, 1-((double)distance)/1000.0));
		if(debug) System.out.println("Score for hit: "+hit.getQueryIdx()+" "+score+" start. "+start+" distance: "+distance);
		return score;*/
		return 1;
	}
	public void printMatrix(int[][] scores) {
		for(int i=0;i<scores.length;i++) {
			for(int j=0;j<scores[i].length;j++) {
				if(j>0) System.out.print("\t");
				System.out.print(""+scores[i][j]);
			}
			System.out.println();
		}
		
	}
}
