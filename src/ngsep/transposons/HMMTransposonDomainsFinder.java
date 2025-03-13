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
import ngsep.math.PhredScoreHelper;
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
	public List<TransposonDomainAlignment> findDomains(QualifiedSequence qdnaSequence) {
		List<TransposonDomainAlignment> answer = new ArrayList<TransposonDomainAlignment>();
		DNAMaskedSequence seqForward = (DNAMaskedSequence) qdnaSequence.getCharacters();
		DNAMaskedSequence seqReverse = seqForward.getReverseComplement();
		//Map<Integer,String> orfsForward = calculateORFs(seqForward);
		//Map<Integer,String> orfsReverse = calculateORFs(seqReverse);
		Map<Integer,String> orfsForward = calculateORFConcat(seqForward);
		Map<Integer,String> orfsReverse = calculateORFConcat(seqReverse);
		int tlf = calculateTotalLength(orfsForward);
		int tlr = calculateTotalLength(orfsReverse);
		int mlf = calculateMaxLength(orfsForward);
		int mlr = calculateMaxLength(orfsReverse);
		//System.out.println("Forward orfs: "+orfsForward.size()+" length: " +tlf+" max "+mlf+" Reverse orfs: "+orfsReverse.size()+" length: " +tlr+" max "+mlr);
		//Map<Integer,String> orfs = orfsForward;
		//boolean reverse = tlr>tlf;
		//boolean reverse = true;
		//if(reverse) orfs = orfsReverse;
		//answer.addAll(findDomainsFromORFs(qdnaSequence, orfs, reverse));
		List<TransposonDomainAlignment> domainsForward = findDomainsFromORFs(qdnaSequence, orfsForward, false);
		List<TransposonDomainAlignment> domainsReverse = findDomainsFromORFs(qdnaSequence, orfsReverse, true);
		answer = selectStrand(domainsForward,domainsReverse);
		return sortAndFilterDomains(answer);
	}
	
	private List<TransposonDomainAlignment> findDomainsFromORFs(QualifiedSequence qdnaSequence, Map<Integer, String> orfs, boolean reverse) {
		List<TransposonDomainAlignment> answer = new ArrayList<TransposonDomainAlignment>();
		for(Map.Entry<Integer, String> orf:orfs.entrySet()) {
			int startDNA = orf.getKey();
			String aaSeq = orf.getValue();
			//System.out.println("Testing of at start: "+startDNA+" Seq len: "+aaSeq.length()+" seq: "+aaSeq+ " reverse: "+reverse);
			for (ProfileAlignmentHMM hmm :hmms) {
				//System.out.println("Testing orf at start: "+startDNA+" Next HMM: "+hmm.getId()+" code "+hmm.getDomainCode()+ " length "+ hmm.getSteps());
				ProfileAlignmentDomain aaDomain = hmm.findDomain(aaSeq);
				if(aaDomain!=null) {
					TransposonDomainAlignment domain = buildTEDomain(qdnaSequence, startDNA, aaDomain, reverse);
					//System.out.println("Found domain at: "+startDNA+" HMM: "+hmm.getId()+" code "+hmm.getDomainCode()+" protein limits: "+aaDomain.getStart()+" "+aaDomain.getEnd()+" pos domain: "+domain.getStart());
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
	private Map<Integer, String> calculateORFConcat(CharSequence seq) {
		Map<Integer, String> answer = new LinkedHashMap<Integer, String>();
		ProteinTranslator translator = ProteinTranslator.getInstance();
		int n = seq.length();
		int i = 0;
		int iA = 0;
		StringBuilder orfConcatenated = new StringBuilder();
		Map<Integer, String> independentORFs = calculateORFs(seq);
		for(Map.Entry<Integer, String> entry:independentORFs.entrySet()) {
			int start = entry.getKey();
			String orf = entry.getValue();
			if(start  - i > 20) {
				String orfInternal = calculateLongestORF(seq.subSequence(i, start));
				orfConcatenated.append(orfInternal);
			}
			if(i==0) iA =i-orfConcatenated.length(); 
			i=start;
			orfConcatenated.append(orf);
			i+=3*orf.length();
		}
		answer.put(0,orfConcatenated.toString());
		return answer;
	}
	private String calculateLongestORF(CharSequence seq) {
		Map<Integer, String> orfs = calculateORFs(seq);
		String answer = "";
		int max = 0;
		for (String orf:orfs.values()) {
			if(orf.length()>max) {
				max = orf.length();
				answer = orf;
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
	private TransposonDomainAlignment buildTEDomain(QualifiedSequence qdnaSequence, int startDNA, ProfileAlignmentDomain aaDomain, boolean reverseStrand) {
		int start = startDNA+3*aaDomain.getStart();
		int end = startDNA+3*aaDomain.getEnd();
		if(reverseStrand) {
			end = qdnaSequence.getLength() - 1 - startDNA - 3*aaDomain.getStart();
			start = qdnaSequence.getLength() - 1 - startDNA - 3*aaDomain.getEnd(); 
		}
		TransposonDomainAlignment domain = new TransposonDomainAlignment(qdnaSequence, start, end, aaDomain);
		domain.setReverse(reverseStrand);
		return domain;
	}
	private List<TransposonDomainAlignment> selectStrand(List<TransposonDomainAlignment> domainsForward, List<TransposonDomainAlignment> domainsReverse) {
		double scoreF = calculateScore(domainsForward);
		double scoreR = calculateScore(domainsReverse);
		if(scoreF>scoreR) return domainsForward;
		return domainsReverse;
	}
	private double calculateScore(List<TransposonDomainAlignment> domains) {
		double score = 0;
		for(TransposonDomainAlignment aln:domains) {
			score+= PhredScoreHelper.calculatePhredScore(aln.getEvalue());
		}
		return score;
	}
	private List<TransposonDomainAlignment> sortAndFilterDomains(List<TransposonDomainAlignment> domains) {
		List<TransposonDomainAlignment> answer = new ArrayList<TransposonDomainAlignment>(domains.size());
		Collections.sort(domains,(d1,d2)-> d1.getStart()-d2.getStart());
		TransposonDomainAlignment last = null;
		for(TransposonDomainAlignment aln:domains) {
			//System.out.println("Next aln: "+aln.getDomainCode()+" "+aln.getHmmID()+" limits: "+aln.getStart()+" "+aln.getEnd()+" evalue: "+aln.getEvalue());
			if(!overlap(last,aln)) {
				answer.add(aln);
				last = aln;
			} else if (last.getEvalue()>aln.getEvalue()) {
				answer.set(answer.size()-1, aln);
				last = aln;
			}
		}
		return answer;
	}
	private boolean overlap(TransposonDomainAlignment last, TransposonDomainAlignment aln) {
		if(last==null) return false;
		int overlap = last.getEnd()-aln.getStart();
		return overlap > 0.3*aln.getLength();
	}

	public static void main(String[] args) throws Exception {
		String filename = args[0];
		String domainsName = args[1];
		String domainsMetadata = args[2];
		String outPrefix = args[3];
		String domainsOutput=outPrefix+"_domainsOutputnew.txt";
		String familiesOutput=outPrefix+"_familiesOutputnew.txt";
		HMMTransposonDomainsFinder instance = new HMMTransposonDomainsFinder();
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> sequences = handler.loadSequences(filename);
		if (sequences.size() == 0) throw new Exception("No sequences found in file: " + filename);
		//ProteinNullModel nullModel = new ProteinNullModel();
		NaiveProteinNullModel nullModel = new NaiveProteinNullModel();
        
		// Load available HMMs
		File domainsFile = new File(domainsName);
		if(!domainsFile.exists()) throw new Exception("The file " + domainsFile+" was not found");
		String [] hmmFilenames;
		if(domainsFile.isDirectory()) {
			hmmFilenames = domainsFile.list((d,n)->n.endsWith(".hmm"));
			for(int i=0;i<hmmFilenames.length;i++) hmmFilenames[i] = domainsName+File.separator+hmmFilenames[i];
		} else {
			hmmFilenames = new String [1];
			hmmFilenames[0] = domainsName;
		}
		ProfileAlignmentHMMLoader hmmLoader = new ProfileAlignmentHMMLoader(nullModel);
		hmmLoader.loadDomainCodes(domainsMetadata);
		for(String hmmFilename:hmmFilenames) {
			List<ProfileAlignmentHMM> hmmsFile = hmmLoader.loadHMMs(hmmFilename);
			for(ProfileAlignmentHMM hmm:hmmsFile) System.out.println("Loaded hmm: "+hmm.getId()+" name: "+hmm.getName()+" domainCode: "+hmm.getDomainCode());
			instance.hmms.addAll(hmmsFile);
		}
		long start = System.currentTimeMillis();
		try (PrintStream outDomains = new PrintStream(domainsOutput);
			 PrintStream outFamilies = new PrintStream(familiesOutput)) {
			// Escribir el encabezado
			outDomains.println("id\tStart\tLength\tClass\tProfileID\tE-value");
			outFamilies.println("id\tOrder\tFamily");
			int i=0;
			for (QualifiedSequence seq : sequences) {
				System.out.println(seq.getName()+" "+i);
				List<TransposonDomainAlignment> domains = instance.findDomains(seq);
				instance.printDomains(domains,outDomains);
				if(domains.size()>0) {
					//TODO: Better calculation of strand
					if(domains.get(0).isReverse()) {
						Collections.reverse(domains);
					}
					TransposableElementFamily family = TransposableElementFamily.matchFamily(domains);
					if(family!=null) {
						outFamilies.println(seq.getName()+"\t"+family.getOrder()+"\t"+family.getId());
					} else {
						outFamilies.println(seq.getName()+"\tUnknown\tUnknown");
					}
				}
				i++;
			}
		}
		long end = System.currentTimeMillis();
        long time = end - start;
        System.out.println("Process finished. Time (s): " + (time/1000.0));

    }
	private void printDomains(List<TransposonDomainAlignment> domains, PrintStream out) {
		for(TransposonDomainAlignment domain:domains) {
			out.print(domain.getQseq().getName());
			out.print("\t"+domain.getQseq().getLength());
			out.print("\t"+domain.getStart());
			out.print("\t"+domain.getLength());
			out.print("\t"+(domain.isReverse()?"-":"+"));
			out.print("\t"+domain.getDomainCode());
			out.print("\t"+domain.getHmmID());
			out.println("\t"+domain.getEvalue());
		}
	}
}
class ProteinNullModel implements ProfileAlignmentNullModel {
	private final Map<Character, Double> aminoAcidLogProbs;
    //private final double log2 = Math.log(2);

    public ProteinNullModel() {
        aminoAcidLogProbs = new HashMap<>();
        aminoAcidLogProbs.put('A', Math.log10(0.0777));
        aminoAcidLogProbs.put('C', Math.log10(0.0157));
        aminoAcidLogProbs.put('D', Math.log10(0.053));
        aminoAcidLogProbs.put('E', Math.log10(0.0656));
        aminoAcidLogProbs.put('F', Math.log10(0.0405));
        aminoAcidLogProbs.put('G', Math.log10(0.0704));
        aminoAcidLogProbs.put('H', Math.log10(0.0231));
        aminoAcidLogProbs.put('I', Math.log10(0.0484));
        aminoAcidLogProbs.put('K', Math.log10(0.0692));
        aminoAcidLogProbs.put('L', Math.log10(0.096));
        aminoAcidLogProbs.put('M', Math.log10(0.0238));
        aminoAcidLogProbs.put('N', Math.log10(0.0427));
        aminoAcidLogProbs.put('P', Math.log10(0.0469));
        aminoAcidLogProbs.put('Q', Math.log10(0.0393));
        aminoAcidLogProbs.put('R', Math.log10(0.0526));
        aminoAcidLogProbs.put('S', Math.log10(0.0694));
        aminoAcidLogProbs.put('T', Math.log10(0.055));
        aminoAcidLogProbs.put('V', Math.log10(0.0667));
        aminoAcidLogProbs.put('W', Math.log10(0.0118));
        aminoAcidLogProbs.put('Y', Math.log10(0.0311));
    }

    public Double calculateScore(String sequence) {
        double score = 0.0;
        for (char aa : sequence.toCharArray()) {
        	//System.out.print(aa);
            Double logProb = aminoAcidLogProbs.get(aa);
            if (logProb != null) {
                score += logProb;
            }
        }
        return score;
    }

    //public double calculateBitScore(double hmmScore, String sequence) {
    //    double nullScore = calculateNullScore(sequence);
    //    return (hmmScore - nullScore) / log2;
    //}
}
class NaiveProteinNullModel implements ProfileAlignmentNullModel {
	private final Map<Character, Double> aminoAcidLogProbs;
    //private final double log2 = Math.log(2);

    public NaiveProteinNullModel() {
        aminoAcidLogProbs = new HashMap<>();
        String alphabet = "ACDEFGHIKLMNPQRSTVWY";
        for(int i=0;i<alphabet.length();i++) {
        	aminoAcidLogProbs.put(alphabet.charAt(i), Math.log10(1.0/alphabet.length()));
        }
    }

    public Double calculateScore(String sequence) {
        double score = 0.0;
        for (char aa : sequence.toCharArray()) {
        	//System.out.print(aa);
            Double logProb = aminoAcidLogProbs.get(aa);
            if (logProb != null) {
                score += logProb;
            }
        }
        return score;
    }

    //public double calculateBitScore(double hmmScore, String sequence) {
    //    double nullScore = calculateNullScore(sequence);
    //    return (hmmScore - nullScore) / log2;
    //}
}
