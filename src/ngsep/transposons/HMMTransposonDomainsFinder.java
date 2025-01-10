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

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import ngsep.hmm.ProfileAlignmentDomain;
import ngsep.hmm.ProfileAlignmentHMM;
import ngsep.hmm.ProfileAlignmentNullModel;
import ngsep.hmm.io.ProfileAlignmentHMMLoader;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.transcriptome.ProteinTranslator;

/**
 * @author Leonidas Villamil
 * @author Jorge Duitama
 */
public class HMMTransposonDomainsFinder {
	private List<ProfileAlignmentHMM> hmms = new ArrayList<ProfileAlignmentHMM>();
	public List<TransposonDomain> findDomains(QualifiedSequence qdnaSequence) {
		List<TransposonDomain> answer = new ArrayList<TransposonDomain>();
		DNAMaskedSequence seqForward = (DNAMaskedSequence) qdnaSequence.getCharacters();
		DNAMaskedSequence seqReverse = seqForward.getReverseComplement();
		Map<Integer,String> orfsForward = calculateORFs(seqForward);
		Map<Integer,String> orfsReverse = calculateORFs(seqReverse);
		int tlf = calculateTotalLength(orfsForward);
		int tlr = calculateTotalLength(orfsReverse);
		int mlf = calculateMaxLength(orfsForward);
		int mlr = calculateMaxLength(orfsReverse);
		System.out.println("Forward orfs: "+orfsForward.size()+" length: " +tlf+" max "+mlf+" Reverse orfs: "+orfsReverse.size()+" length: " +tlr+" max "+mlr);
		//Map<Integer,String> orfs = orfsForward;
		//boolean reverse = tlr>tlf;
		//boolean reverse = true;
		//if(reverse) orfs = orfsReverse;
		//answer.addAll(findDomainsFromORFs(qdnaSequence, orfs, reverse));
		answer.addAll(findDomainsFromORFs(qdnaSequence, orfsForward, false));
		answer.addAll(findDomainsFromORFs(qdnaSequence, orfsReverse, true));
		Collections.sort(answer,(d1,d2)-> d1.getStart()-d2.getStart());
		return answer;
	}
	private List<TransposonDomain> findDomainsFromORFs(QualifiedSequence qdnaSequence, Map<Integer, String> orfs, boolean reverse) {
		List<TransposonDomain> answer = new ArrayList<TransposonDomain>();
		for(Map.Entry<Integer, String> orf:orfs.entrySet()) {
			int startDNA = orf.getKey();
			String aaSeq = orf.getValue();
			//System.out.println("Testing orf at start: "+startDNA+" Seq len: "+aaSeq.length()+" seq: "+aaSeq+ " reverse: "+reverse);
			for (ProfileAlignmentHMM hmm :hmms) {
				//System.out.println("Testing orf at start: "+startDNA+" Next HMM: "+hmm.getId());
				ProfileAlignmentDomain aaDomain = hmm.findDomain(aaSeq);
				if(aaDomain!=null) {
					TransposonDomain domain = buildTEDomain(qdnaSequence, startDNA, aaDomain, reverse);
					answer.add(domain);
				}
			}
		}
		return answer;
	}
	private Map<Integer, String> calculateORFs(CharSequence seq) {
		Map<Integer, String> answer = new LinkedHashMap<Integer, String>();
		ProteinTranslator translator = ProteinTranslator.getInstance();
		int n = seq.length();
		int i=0;
		while(i<n) {
			String orf = translator.getProteinSequence(seq, i);
			if(orf.length()>50) {
				answer.put(i, orf);
				i+=3*orf.length();
			} else {
				i++;
			}
		}
		return answer;
	}
	private int calculateTotalLength(Map<Integer, String> orfs) {
		int total = 0;
		for(String orf:orfs.values()) total+=orf.length();
		return total;
	}
	private int calculateMaxLength(Map<Integer, String> orfs) {
		int max = 0;
		for(String orf:orfs.values()) max=Math.max(max, orf.length());
		return max;
	}
	private TransposonDomain buildTEDomain(QualifiedSequence qdnaSequence, int startDNA, ProfileAlignmentDomain aaDomain, boolean reverseStrand) {
		int start = startDNA+3*aaDomain.getStart();
		int end = startDNA+3*aaDomain.getEnd();
		if(reverseStrand) {
			end = qdnaSequence.getLength() - 1 - startDNA - 3*aaDomain.getStart();
			start = qdnaSequence.getLength() - 1 - startDNA - 3*aaDomain.getEnd(); 
		}
		TransposonDomain domain = new TransposonDomain(qdnaSequence, start, end, aaDomain);
		domain.setReverse(reverseStrand);
		return domain;
	}
	public static void main(String[] args) throws Exception {
		String filename = args[0];
		String domainsDir = args[1];
		String outPrefix = args[2];
		String domainsOutput=outPrefix+"_domainsOutputnew.txt";
		String familiesOutput=outPrefix+"_familiesOutputnew.txt";
		HMMTransposonDomainsFinder instance = new HMMTransposonDomainsFinder();
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> sequences = handler.loadSequences(filename);
		if (sequences.size() == 0) throw new Exception("No sequences found in file: " + filename);
		ProteinNullModel nullModel = new ProteinNullModel();
        
		// Load available HMMs
		String [] hmmFilenames = {"GAGPF03732.hmm","APPF13650.hmm","ENVPF00517.hmm","INTPF00665.hmm","RHPF00075.hmm","RTPF00078.hmm"};
		for(String hmmFilename:hmmFilenames) {
			ProfileAlignmentHMM hmm = ProfileAlignmentHMMLoader.loadHMM(domainsDir+File.separator+hmmFilename,nullModel);
			instance.hmms.add(hmm);
		}
		
		long start = System.currentTimeMillis();
		try (PrintStream outDomains = new PrintStream(domainsOutput);
			 PrintStream outFamilies = new PrintStream(familiesOutput)) {
			// Escribir el encabezado
			outDomains.println("id\tStart\tLength\tClass\tE-value");
			outFamilies.println("id\tStart\tLength\tFamily");
			int i=0;
			for (QualifiedSequence seq : sequences) {
				System.out.println(seq.getName()+" "+i);
				List<TransposonDomain> domains = instance.findDomains(seq);
				instance.printDomains(domains,outDomains);
				i++;
			}
		}
		long end = System.currentTimeMillis();
        long time = end - start;
        System.out.println("Process finished. Time (s): " + (time/1000.0));

    }
	private void printDomains(List<TransposonDomain> domains, PrintStream out) {
		for(TransposonDomain domain:domains) {
			out.print(domain.getQseq().getName());
			out.print("\t"+domain.getQseq().getLength());
			out.print("\t"+domain.getStart());
			out.print("\t"+domain.getLength());
			out.print("\t"+(domain.isReverse()?"-":"+"));
			out.print("\t"+domain.getAlnDomain().getHmmID());
			out.println("\t"+domain.getAlnDomain().getEvalue());
		}
		
		
	}
}
class ProteinNullModel implements ProfileAlignmentNullModel {
	private final Map<Character, Double> aminoAcidProbabilities;
    //private final double log2 = Math.log(2);

    public ProteinNullModel() {
        aminoAcidProbabilities = new HashMap<>();
        aminoAcidProbabilities.put('A', 0.0777);
        aminoAcidProbabilities.put('C', 0.0157);
        aminoAcidProbabilities.put('D', 0.053);
        aminoAcidProbabilities.put('E', 0.0656);
        aminoAcidProbabilities.put('F', 0.0405);
        aminoAcidProbabilities.put('G', 0.0704);
        aminoAcidProbabilities.put('H', 0.0231);
        aminoAcidProbabilities.put('I', 0.0484);
        aminoAcidProbabilities.put('K', 0.0692);
        aminoAcidProbabilities.put('L', 0.096);
        aminoAcidProbabilities.put('M', 0.0238);
        aminoAcidProbabilities.put('N', 0.0427);
        aminoAcidProbabilities.put('P', 0.0469);
        aminoAcidProbabilities.put('Q', 0.0393);
        aminoAcidProbabilities.put('R', 0.0526);
        aminoAcidProbabilities.put('S', 0.0694);
        aminoAcidProbabilities.put('T', 0.055);
        aminoAcidProbabilities.put('V', 0.0667);
        aminoAcidProbabilities.put('W', 0.0118);
        aminoAcidProbabilities.put('Y', 0.0311);
    }

    public Double calculateScore(String sequence) {
        double score = 0.0;
        for (char aa : sequence.toCharArray()) {
        	//System.out.print(aa);
            Double prob = aminoAcidProbabilities.get(aa);
            if (prob != null) {
                score += Math.log(prob);
            }
        }
        return score;
    }

    //public double calculateBitScore(double hmmScore, String sequence) {
    //    double nullScore = calculateNullScore(sequence);
    //    return (hmmScore - nullScore) / log2;
    //}
}
