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
package ngsep.transcriptome;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;

/**
 * Implements a transcriptome that handles sequence names and coordinates in a genomic fashion.
 * This calls is useful to transform coordinates relative to a transcript library to coordinates
 * relative to the corresponding genomic assembly. 
 * @author Jorge Duitama
 */
public class Transcriptome {
	public static final String ANNOTATION_CODING="Coding";
	public static final String ANNOTATION_INTRON="Intron";
	public static final String ANNOTATION_INTERGENIC="Intergenic";
	public static final String ANNOTATION_5P_UTR="FivePrimeUTR";
	public static final String ANNOTATION_3P_UTR="ThreePrimeUTR";
	public static final String ANNOTATION_UPSTREAM="Upstream";
	public static final String ANNOTATION_DOWNSTREAM="Downstream";
	public static final String ANNOTATION_NONCODINGRNA="NCRNA";
	public static final String ANNOTATION_SYNONYMOUS="Synonymous";
	public static final String ANNOTATION_MISSENSE="Missense";
	public static final String ANNOTATION_NONSENSE="Nonsense";
	public static final String ANNOTATION_FRAMESHIFT="Frameshift";
	public static final String ANNOTATION_JUNCTION="ExonJunction";
	
	
	//Genes indexed by id
	Map<String, Gene> genesMap = new TreeMap<String, Gene>();
	//Transcripts indexed by id
	private Map<String, Transcript> transcriptsMap = new TreeMap<String, Transcript>();
	//Transcripts sorted by absolute position and indexed by chromosome
	private GenomicRegionSortedCollection<Transcript> sortedTranscripts;
	//Transcripts sorted by absolute position and indexed by chromosome
	private GenomicRegionSortedCollection<Exon> sortedExons;
	//Transcripts by gene
	private Map<String, List<Transcript>> transcriptsByGene = new TreeMap<String, List<Transcript>>();
	//RNA to protein translator
	private ProteinTranslator proteinTranslator = new ProteinTranslator();
	
	public Transcriptome() {
		sortedTranscripts = new GenomicRegionSortedCollection<Transcript>();
		sortedExons = new GenomicRegionSortedCollection<Exon>();
	}
	public Transcriptome (QualifiedSequenceList sequenceNames) {
		sortedTranscripts = new GenomicRegionSortedCollection<Transcript>(sequenceNames);
		sortedExons = new GenomicRegionSortedCollection<Exon>(sequenceNames);
	}
	public Gene getGene (String id) {
		return genesMap.get(id);
	}
	/**
	 * Adds the given transcript to the transcriptome
	 * @param t New transcript
	 */
	public void addTranscript (Transcript t) {
		if(transcriptsMap.get(t.getId()) == null) {
			Gene g = t.getGene();
			List<Transcript> transcriptsG = transcriptsByGene.get(g.getId());
			if(g!=null && genesMap.get(g.getId())==null) {
				genesMap.put(g.getId(), g);
				transcriptsG = new ArrayList<Transcript>();
				transcriptsByGene.put(g.getId(), transcriptsG);
			}
			transcriptsG.add(t);
			transcriptsMap.put(t.getId(), t);
			sortedTranscripts.add(t);
			sortedExons.addAll(t.getExons());
		}
	}
	/**
	 * Assigns the given sequence to the transcript with the given id
	 * @param transcriptId Id of the transcript to assign
	 * @param sequence Sequence to associate with the transcript
	 */
	public void fillTranscript(String transcriptId, DNAMaskedSequence sequence) {
		Transcript transcript = transcriptsMap.get(transcriptId);
		if(transcript!=null) {
			transcript.setCDNASequence(sequence);
		} else {
			System.err.println("WARN: Transcript not found with id: "+transcriptId);
		}
	}
	/**
	 * Returns the transcript with the given id
	 * @param id of the transcript to retireve
	 * @return Transcript with the given id or null if it is not found
	 */
	public Transcript getTranscript(String id) {
		return transcriptsMap.get(id);
	}
	/**
	 * Retrieves all transcripts in the given sequenceName
	 * @param sequenceName Sequence name of the region to look for
	 * @return GenomicRegionSortedCollection<Transcript> All transcripts in the given genomic sequence
	 */
	public GenomicRegionSortedCollection<Transcript> getTranscripts (String sequenceName) {
		return sortedTranscripts.getSequenceRegions(sequenceName);
	}
	/**
	 * Retrieves the transcripts spanning the given position
	 * @param sequenceName Sequence name of the region to look for
	 * @param position Genomic location to look for 
	 * @return GenomicRegionSortedCollection<Transcript> Transcripts spanning the given coordinate
	 */
	public GenomicRegionSortedCollection<Transcript> getTranscripts (String sequenceName, int position) {
		return sortedTranscripts.findSpanningRegions(sequenceName, position);
	}
	/**
	 * Retrieves the transcripts spanning the region delimited by the given coordinates
	 * @param sequenceName Sequence name of the region to look for
	 * @param first First genomic coordinate of the region to look for 
	 * @param last Last genomic  coordinate of the region to look for
	 * @return GenomicRegionSortedCollection<Transcript> Transcripts spanning the given region
	 */
	public GenomicRegionSortedCollection<Transcript> getTranscripts (String sequenceName, int first, int last) {
		return sortedTranscripts.findSpanningRegions(sequenceName, first, last);
	}
	/**
	 * Retrieves the transcripts spanning the given genomic region
	 * @param region Region to look for
	 * @return GenomicRegionSortedCollection<Transcript> Transcripts spanning the given region
	 */
	public GenomicRegionSortedCollection<Transcript> getTranscripts (GenomicRegion region) {
		return sortedTranscripts.findSpanningRegions(region);
	}
	/**
	 * Return the transcripts for a gene with the given id
	 * @param geneId Id of the gene
	 * @return List<Transcript> transcripts of the gene with the given id
	 */
	public List<Transcript> getTranscriptsByGene(String geneId) {
		return transcriptsByGene.get(geneId);
	}
	/**
	 * Retrieves all exons in the given sequence
	 * @param sequenceName Name of the sequence to look for
	 * @return GenomicRegionSortedCollection<Exon> All exons in the given genomic sequence
	 */
	public GenomicRegionSortedCollection<Exon> getExons(String sequenceName) {
		return sortedExons.getSequenceRegions(sequenceName);
	}
	/**
	 * Retrieves all exons spanning a given genomic location
	 * @param sequenceName Name of the sequence to look for 
	 * @param position Genomic location to look for 
	 * @return GenomicRegionSortedCollection<Exon> Exons spanning the given genomic coordinate
	 */
	public GenomicRegionSortedCollection<Exon> getExons(String sequenceName, int position) {
		return sortedExons.findSpanningRegions(sequenceName, position);
	}
	/**
	 * Retrieves all exons spanning the region delimited by the given criteria
	 * @param sequenceName Name of the sequence to look for
	 * @param first First genomic coordinate of the region to look for 
	 * @param last Last genomic  coordinate of the region to look for
	 * @return GenomicRegionSortedCollection<Exon> Exons spanning the genomic region defined by the input
	 * parameters
	 */
	public GenomicRegionSortedCollection<Exon> getExons (String sequenceName, int first, int last) {
		return sortedExons.findSpanningRegions(sequenceName, first, last);
	}
	/**
	 * Retrieves all exons spanning the region delimited by the given criteria
	 * @param genomicRegion Region to look for
	 * @return GenomicRegionSortedCollection<Exon> Exons spanning the given genomic region
	 * parameters
	 */
	public GenomicRegionSortedCollection<Exon> getExons (GenomicRegion region) {
		return sortedExons.findSpanningRegions(region);
	}
	
	/**
	 * @return QualifiedSequenceList Names of the sequences in the transcriptome
	 */
	public QualifiedSequenceList getSequenceNames() {
		return sortedTranscripts.getSequenceNames();
	}


	public char getReferenceBase (String seqName, int absolutePosition) {
		GenomicRegionSortedCollection<Transcript> transcripts = getTranscripts(seqName, absolutePosition); 
		for(Transcript t:transcripts) {
			char base = t.getReferenceBase(absolutePosition);
			if(base!=0) {
				return base;
			}
		}
		return 0;
	}
	public void setReferenceBase(String seqName, int absolutePosition, char base) {
		GenomicRegionSortedCollection<Transcript> transcripts = getTranscripts(seqName, absolutePosition);
		for(Transcript t:transcripts) {
				t.setReferenceBase(absolutePosition, base);
		}
	}
	/**
	 * Retrieves the reference base at the region enclosed by the given parameters
	 * @param sequenceName Name of the sequence to look for
	 * @param first First position to include
	 * @param last Last position to include. It must be greater than first
	 * @return String Sequence of the requested region
	 */
	public String getReference(String sequenceName, int first, int last) {
		GenomicRegionSortedCollection<Transcript> transcripts = getTranscripts(sequenceName, first); 
		for(Transcript t:transcripts) {
			char base = t.getReferenceBase(last);
			if(base!=0) {
				return t.getReference(first, last);
			}
		}
		return null;
	}
	public List<Transcript> getAllTranscripts () {
		return sortedTranscripts.asList();
	}
	/**
	 * Calculates an annotation for the given variant based on its two first alleles
	 * @param variant Genomic variant to annotate
	 * @param offsetUpstream Offset to call the variant upstream of a gene
	 * @param offsetDownstream Offset to call the variant downstream of a gene
	 * @return List<GenomicVariantAnnotation> Information for functional annotations related with the variant 
	 */
	public List<GenomicVariantAnnotation> calculateAnnotation(GenomicVariant variant, int offsetUpstream, int offsetDownstream) {
		String [] alleles = variant.getAlleles();
		String reference = alleles[0];
		String alternative = alleles[1];
		String alternativeR = DNAMaskedSequence.getReverseComplement(alternative);
		
		int differenceBases = alternative.length() - reference.length();
		int expectexProteinIncrease = differenceBases/3;
		int maxOffset = Math.max(offsetUpstream, offsetDownstream);
		List<GenomicVariantAnnotation> annotationCodingExon = new ArrayList<GenomicVariantAnnotation>();
		List<GenomicVariantAnnotation> annotationUTR = new ArrayList<GenomicVariantAnnotation>();
		List<GenomicVariantAnnotation> annotationNCRNA = new ArrayList<GenomicVariantAnnotation>();
		List<GenomicVariantAnnotation> annotationClose = new ArrayList<GenomicVariantAnnotation>();
		List<GenomicVariantAnnotation> annotationIntron = new ArrayList<GenomicVariantAnnotation>();
		for(Transcript t:getTranscripts(variant.getSequenceName(), variant.getFirst()-maxOffset, variant.getLast()+maxOffset)) {
			//System.out.println("Transcript: "+id+". Coding: "+coding+". Reverse: "+reverse);
			Exon exStart = t.getExonByAbsolutePosition(variant.getFirst());
			Exon exEnd = t.getExonByAbsolutePosition(variant.getLast());
			if(exStart!=exEnd) {
				annotationCodingExon = makeAnnotation(variant, ANNOTATION_JUNCTION, t, null);
			} else if (exStart == null) {
				//Cases outside the transcript or in introns
				if(t.getFirst()<=variant.getFirst() && t.getLast()>=variant.getLast()) {
					annotationIntron = makeAnnotation(variant, ANNOTATION_INTRON, t, null);
				} else if (variant.getLast()<t.getFirst() ) {
					//Variant before the transcript in genomic location. Check upstream for positive strand and downstream for negative strand
					if(t.isPositiveStrand() && variant.getLast()>t.getFirst()-offsetUpstream) {
						annotationClose = makeAnnotation(variant, ANNOTATION_UPSTREAM, t, null);
					}
					else if(annotationClose.isEmpty() && t.isNegativeStrand() && variant.getLast()>t.getFirst()-offsetDownstream) {
						annotationClose = makeAnnotation(variant, ANNOTATION_DOWNSTREAM, t, null);
					}
				} else if (variant.getFirst()> t.getLast()) {
					//Variant before the transcript in genomic location. Check upstream for negative strand and downstream for positive strand
					if(annotationClose.isEmpty() && t.isPositiveStrand() && variant.getFirst()<t.getLast()+offsetDownstream) {
						annotationClose = makeAnnotation(variant, ANNOTATION_DOWNSTREAM, t, null);
					}
					else if(t.isNegativeStrand() && variant.getFirst()<t.getLast()+offsetUpstream) {
						annotationClose = makeAnnotation(variant, ANNOTATION_UPSTREAM, t, null);
					}
				}
			} else if(exStart.isCoding()) {
				int transcriptionStart = t.getCodingRelativeStart();
				int absoluteFirst = variant.getFirst();
				if(t.isNegativeStrand()) {
					absoluteFirst = variant.getLast();
				}
				int varTranscriptStart = t.getRelativeTranscriptPosition(absoluteFirst);
				int varTranscriptEnd = varTranscriptStart+(variant.getLast()-variant.getFirst())+1;
				int varCodingStart = varTranscriptStart-transcriptionStart;
				int codon = varCodingStart/3+1;
				int module = varCodingStart%3;
				//System.out.println("Transcript: "+t.getId()+". Relative start: "+varTranscriptStart+". Coding start: "+transcriptionStart+". Codon: "+codon);
				String codonStr = ""+codon+"."+(module+1);
				if(differenceBases%3!=0){
					//Frameshift mutation
					annotationCodingExon = makeAnnotation(variant, ANNOTATION_FRAMESHIFT, t, codonStr);
				} else {
					DNAMaskedSequence cdnaSequence = t.getCDNASequence();
					if(cdnaSequence!=null && varTranscriptStart<cdnaSequence.length()) {
						int startTest = varTranscriptStart - module;
						int endTest = Math.min(cdnaSequence.length(),varTranscriptEnd+3);
						String testReference = cdnaSequence.subSequence(startTest,endTest).toString();
						String testA = alternative;
						if(t.isNegativeStrand()) testA = alternativeR;
						String testVariant = cdnaSequence.subSequence(startTest,varTranscriptStart)+testA;
						if(endTest > varTranscriptEnd) testVariant+=cdnaSequence.subSequence(varTranscriptEnd,endTest);
						//System.out.println("Ref: "+testReference+". Var: "+testVariant+". Transcript start:"+varTranscriptStart);
						String refProt = proteinTranslator.getProteinSequence(testReference);
						String varProt = proteinTranslator.getProteinSequence(testVariant);
						//System.out.println("Ref prot: "+refProt+". Varprot: "+varProt);
						if(refProt.equals(varProt)) {
							if(annotationCodingExon.isEmpty()) {
								annotationCodingExon = makeAnnotation(variant, ANNOTATION_SYNONYMOUS, t, codonStr);
							}
						} else if (refProt.length()+expectexProteinIncrease==varProt.length()) {
							annotationCodingExon = makeAnnotation(variant, ANNOTATION_MISSENSE, t, codonStr);
						} else {
							annotationCodingExon = makeAnnotation(variant, ANNOTATION_NONSENSE, t, codonStr);
						}
					} else if(annotationCodingExon.isEmpty()) {
						annotationCodingExon = makeAnnotation(variant, ANNOTATION_CODING, t, codonStr);
					}
				}
			} else if(exStart.getStatus()==Exon.STATUS_5P_UTR) {
				annotationUTR = makeAnnotation(variant, ANNOTATION_5P_UTR, t, null);
			} else if (annotationUTR.isEmpty() && exStart.getStatus()==Exon.STATUS_3P_UTR) {
				annotationUTR = makeAnnotation(variant, ANNOTATION_3P_UTR, t, null);
			} else if(exStart.getStatus()==Exon.STATUS_NCRNA)  {
				annotationNCRNA = makeAnnotation(variant, ANNOTATION_NONCODINGRNA, t, null);
			}
		}
		//Priority management
		if(!annotationCodingExon.isEmpty()) return annotationCodingExon;
		if (!annotationUTR.isEmpty()) return annotationUTR;
		if (!annotationNCRNA.isEmpty()) return annotationNCRNA;
		if(!annotationClose.isEmpty()) return annotationClose;
		if(!annotationIntron.isEmpty()) return annotationIntron; 
		return makeAnnotation(variant, ANNOTATION_INTERGENIC, null, null);
	}
	private List<GenomicVariantAnnotation> makeAnnotation(GenomicVariant var, String type, Transcript t, String codon) {
		//Remove previous annotations
		List<GenomicVariantAnnotation> answer = new ArrayList<GenomicVariantAnnotation>();
		answer.add(new GenomicVariantAnnotation(var, GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ANNOTATION, type));
		if(t!=null) {
			answer.add(new GenomicVariantAnnotation(var, GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ID, t.getId()));
			answer.add(new GenomicVariantAnnotation(var, GenomicVariantAnnotation.ATTRIBUTE_GENE_NAME, t.getGeneName()));
			if(codon!=null) answer.add(new GenomicVariantAnnotation(var, GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_CODON, codon));
		}
		return answer;
	}
	public void fillSequenceTranscripts(ReferenceGenome genome) {
		for(Transcript t:sortedTranscripts) {
			StringBuilder transcriptSeq = new StringBuilder();
			boolean exonNotFound = false;
			List<Exon> exons = new ArrayList<Exon>(t.getExons());
			if(t.isNegativeStrand()) Collections.reverse(exons);
			for (Exon e: exons) {
				CharSequence genomicSeq = genome.getReference(t.getSequenceName(),e.getFirst(),e.getLast());
				if(genomicSeq == null) {
					System.err.println("WARN: Exon in location "+t.getSequenceName()+":"+e.getFirst()+"-"+ e.getLast()+" not found for transcript: "+t.getId());
					exonNotFound = true;
					break;
				}
				if(t.isNegativeStrand()) {
					genomicSeq = DNAMaskedSequence.getReverseComplement(genomicSeq);
				}
				transcriptSeq.append(genomicSeq);
			}
			if(!exonNotFound) t.setCDNASequence(new DNAMaskedSequence(transcriptSeq.toString().toUpperCase()));
		}
	}
}
