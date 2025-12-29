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
package ngsep.transposons;

import java.util.ArrayList;
import java.util.List;

import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNASequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.KmersMap;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class DeNovoTransposableElementsFinderKmerAbundance implements DeNovoTransposableElementsFinder {

	
	
	/**
	 * Find deNovo transposons given a genome
	 * @param genome to which the transposons will be found 
	 * @return List<TransposableElementAnnotation> list of transposons found with chromosome, starting and ending position
	 * @throws InterruptedException if this exception is caught, it is because for the kmer table the genome could not be loaded
	 */
	@Override
	public List<TransposableElementAnnotation> findTransposons(ReferenceGenome genome) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		
		QualifiedSequenceList list = genome.getSequencesList();
		KmersExtractor tableKmer = new KmersExtractor();
		tableKmer.processQualifiedSequences(list);
		KmersMap kmap= tableKmer.getKmersMap();
		int totalTranspoSize= 0;
	    //traverses the genome, chromosome by chromosome 
		for(QualifiedSequence chromosome:list)
		{
			//gives the chromosome sequence
			System.err.println("Processing sequence: "+chromosome.getName());
			CharSequence chromosomeSequence = chromosome.getCharacters();
			int n = chromosome.getLength();
			int regionStart=0;
			int regionEnd=1;
			int regionFrequency=0;			
			//extract kmer by kmer of the entire genome
			for(int i=0;i<n-15;i++) 
			{
				int end = i+15;
				String kmer = chromosomeSequence.subSequence(i, end).toString();
				//check if kmer is DNA
				if(DNASequence.isDNA(kmer)) 
				{
					int kmerFrequency = kmap.getCount(kmer); 
					//filter kmers with given frequency
					if (kmerFrequency >=10)
					{ 
						int dif= Math.abs(regionEnd-i);
						// report regions
						// does the kmer belongs to the region?
						if (dif<=50)
						{
							regionEnd=end;
							
							// the frequency of the region will be the lowest frequency of the kmers
							if(regionFrequency==0)
							{
								regionFrequency= kmerFrequency;
							}
							if (kmerFrequency<regionFrequency)
							{
								regionFrequency= kmerFrequency;
							}
							// validate if the kmer is the last int the chromosome 
							if(regionEnd==n-1)
							{
								int sizeRegion= regionEnd-regionStart;
								totalTranspoSize+= sizeRegion;
								TransposableElementAnnotation transposon = new TransposableElementAnnotation(chromosome.getName(), regionStart, regionEnd); 
								answer.add(transposon);
								//System.out.println("" + chromosome.getName()+ "   " + regionStart + " "  + regionEnd + "   " +  regionFrequency + "  " + sizeRegion) ;
								// update 
								regionStart=0;
								regionEnd= 1;
								regionFrequency=0;	
							}
						}
						// find a new region	
						else 
						{
							int sizeRegion= regionEnd-regionStart;
							totalTranspoSize+= sizeRegion;
							TransposableElementAnnotation transposon = new TransposableElementAnnotation(chromosome.getName(), regionStart, regionEnd); 
							answer.add(transposon);
							//System.out.println("" + chromosome.getName()+ "   " + regionStart + " "  + regionEnd + "   " +  regionFrequency + "  " + sizeRegion) ;	
							// update
							regionStart=i;
							regionEnd= end;
							regionFrequency=kmerFrequency;	
						}
					
					}
				}
			
			}
 	
		}
		//System.out.println(totalTranspoSize); 	
		return answer;
	}
	public static void main(String[] args) throws Exception {
		ReferenceGenome genome = new ReferenceGenome(args[0]);
		DeNovoTransposableElementsFinderKmerAbundance instance = new DeNovoTransposableElementsFinderKmerAbundance();
		List<TransposableElementAnnotation> anns = instance.findTransposons(genome);
		for(TransposableElementAnnotation ann:anns) {
			System.out.println(""+ann.getSequenceName()+"\t"+ann.getFirst()+"\t"+ann.getLast()+"\t"+ann.getCount());
		}
	}
}
