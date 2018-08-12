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

/**
 * Implements a transcriptome that handles sequence names and coordinates in a genomic fashion.
 * This calls is useful to transform coordinates relative to a transcript library to coordinates
 * relative to the corresponding genomic assembly. 
 * @author Jorge Duitama
 */
public class Transcriptome {
	
	private QualifiedSequenceList sequenceNames;
	//Genes indexed by id
	private Map<String, Gene> genesMap = new TreeMap<String, Gene>();
	//Transcripts indexed by id
	private Map<String, Transcript> transcriptsMap = new TreeMap<String, Transcript>();
	//Transcripts sorted by absolute position and indexed by chromosome
	private GenomicRegionSortedCollection<Transcript> sortedTranscripts;
	
	//Transcripts by gene
	private Map<String, List<Transcript>> transcriptsByGene = new TreeMap<String, List<Transcript>>();
	//RNA to protein translator
	private ProteinTranslator proteinTranslator = new ProteinTranslator();
	
	public Transcriptome (QualifiedSequenceList sequenceNames) {
		this.sequenceNames = sequenceNames;
		sortedTranscripts = new GenomicRegionSortedCollection<Transcript>(sequenceNames);
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
	 * @return QualifiedSequenceList Names of the sequences in the transcriptome
	 */
	public QualifiedSequenceList getSequenceNames() {
		return sequenceNames;
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
	 * Calculates the annotations for the given variant based on their alternative alleles
	 * @param variant Genomic variant to annotate
	 * @param parameters Object with the parameters to perform the annotation
	 * @return List<GenomicVariantAnnotation> Functional annotations of the effect of the alternative alleles 
	 */
	public List<VariantFunctionalAnnotation> calculateAnnotations(GenomicVariant variant, VariantAnnotationParameters parameters) {
		List<VariantFunctionalAnnotation> annotations = new ArrayList<>();
		int offsetUpstream = parameters.getOffsetUpstream();
		int offsetDownstream = parameters.getOffsetDownstream();
		int maxOffset = Math.max(offsetUpstream, offsetDownstream);
		for(Transcript t:getTranscripts(variant.getSequenceName(), variant.getFirst()-maxOffset, variant.getLast()+maxOffset)) {
			//if(variant.getFirst()==1096) System.err.println("Transcript: "+t.getId()+". Coding: "+t.isCoding()+". Reverse: "+t.isNegativeStrand()+" at "+t.getSequenceName()+": "+t.getFirst()+"-"+t.getLast());
			TranscriptSegment segmentStart = t.getTranscriptSegmentByAbsolutePosition(variant.getFirst());
			TranscriptSegment segmentEnd = t.getTranscriptSegmentByAbsolutePosition(variant.getLast());
			if(segmentStart!=segmentEnd) {
				VariantFunctionalAnnotation annotationSplice = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_SPLICE_REGION);
				annotationSplice.setTranscript(t);
				annotations.add(annotationSplice);
			} else if (segmentStart == null) {
				//Cases outside the transcript or in introns
				if(t.getFirst()<=variant.getFirst() && t.getLast()>=variant.getLast()) {
					VariantFunctionalAnnotation annotationIntron = makeIntronAnnotation(variant, t, parameters);
					annotations.add(annotationIntron);
				} else {
					VariantFunctionalAnnotation annotationClose = makeAnnotationClose(variant, t, offsetUpstream, offsetDownstream);
					if(annotationClose!=null) annotations.add(annotationClose);
				}
					
			} else {
				if(segmentStart.isCoding()) {
					annotations.addAll(makeCodingAnnotations(variant, segmentStart, parameters));
				} else {
					VariantFunctionalAnnotation annotation =null;
					if(segmentStart.getStatus()==TranscriptSegment.STATUS_5P_UTR) {
						annotation = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_5P_UTR);
					} else if (segmentStart.getStatus()==TranscriptSegment.STATUS_3P_UTR) {
						annotation = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_3P_UTR);
					} else if(segmentStart.getStatus()==TranscriptSegment.STATUS_NCRNA)  {
						annotation = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_NONCODINGRNA);
					}
					if(annotation!=null) {
						annotation.setTranscript(t);
						annotations.add(annotation);
					}
				}
				//Check Exon splice variant
				VariantFunctionalAnnotation annotationExonSplice = makeAnnotationExonSplice(variant, segmentStart, parameters.getSpliceRegionExonOffset());
				if(annotationExonSplice!=null) annotations.add(annotationExonSplice);
			}
		}
		if(annotations.size()==0) {
			annotations.add(new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_INTERGENIC));
		}
		return annotations;
	}
	private VariantFunctionalAnnotation makeAnnotationExonSplice(GenomicVariant variant, TranscriptSegment segment, int offset) {
		int diffFirst = variant.getFirst() - segment.getFirst() + 1;
		int diffLast = segment.getLast() - variant.getLast() + 1;
		if(diffFirst>offset && diffLast>offset) return null;
		if((diffFirst<=offset && segment.hasIntronLeft()) || (diffLast<=offset && segment.hasIntronRight()) ) {
			VariantFunctionalAnnotation answer = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_EXONIC_SPLICE_REGION);
			answer.setTranscript(segment.getTranscript());
			return answer;
		}
		return null;
	}
	private List<VariantFunctionalAnnotation> makeCodingAnnotations(GenomicVariant variant, TranscriptSegment segment, VariantAnnotationParameters parameters) {
		List<VariantFunctionalAnnotation> annotationsCoding = new ArrayList<>();
		String [] alleles = variant.getAlleles();
		String reference = alleles[0];
		Transcript t = segment.getTranscript();
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
		DNAMaskedSequence cdnaSequence = t.getCDNASequence();
		if(cdnaSequence==null || varTranscriptStart>=cdnaSequence.length()) {
			VariantFunctionalAnnotation annotationCoding = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_CODING);
			annotationCoding.setTranscript(t);
			annotationCoding.setCodonNumber(codon);
			annotationCoding.setCodonPosition((byte) (module+1));
			annotationsCoding.add(annotationCoding);
			return annotationsCoding;
		}
		int startTest = varTranscriptStart - module;
		int endTest = Math.min(cdnaSequence.length(),varTranscriptEnd+3);
		String testReference = cdnaSequence.subSequence(startTest,endTest).toString();
		for(int i=1;i<alleles.length;i++) {
			VariantFunctionalAnnotation annotationCoding;
			byte altAlleleIdx = (byte) i;
			String alternative = alleles[altAlleleIdx];
			String alternativeR = DNAMaskedSequence.getReverseComplement(alternative);
			int differenceBases = alternative.length() - reference.length();
			int expectedProteinIncrease = differenceBases/3;
			String testA = alternative;
			if(t.isNegativeStrand()) testA = alternativeR;
			String testVariant = cdnaSequence.subSequence(startTest,varTranscriptStart)+testA;
			if(endTest > varTranscriptEnd) testVariant+=cdnaSequence.subSequence(varTranscriptEnd,endTest);
			//System.out.println("Ref: "+testReference+". Var: "+testVariant+". Transcript start:"+varTranscriptStart);
			String refProt = proteinTranslator.getProteinSequence(testReference);
			String varProt = proteinTranslator.getProteinSequence(testVariant);
			String change = refProt;
			if(change.length()==0) change="*";
			change+=codon;
			if(varProt.length()>0) change+=varProt;
			else change+="*";
			//System.out.println("Ref prot: "+refProt+". Varprot: "+varProt+" Change: "+change );
			if (differenceBases!=0) {
				//Indel events
				if(Math.abs(differenceBases)%3!=0){
					//Frameshift mutation
					annotationCoding = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT);
				} else {
					annotationCoding = new VariantFunctionalAnnotation(variant, differenceBases>0?VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS:VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL);
				}
			} else if(refProt.equals(varProt)) {
				annotationCoding = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS);
			} else if (refProt.length()+expectedProteinIncrease==varProt.length()) {
				//TODO: Ask for a start codon and not just for an M
				if(startTest==0 && refProt.charAt(0)=='M' && (varProt.length()==0 || varProt.charAt(0)!='M')) {
					annotationCoding = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_START_LOST);
				} else {
					annotationCoding = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_MISSENSE);
				}
			} else if (refProt.length()==0 && varProt.length()>0) {
				annotationCoding = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST);
			} else {
				annotationCoding = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_NONSENSE);
			}
			annotationCoding.setAltAlleleIdx(altAlleleIdx);
			annotationCoding.setTranscript(t);
			annotationCoding.setCodonNumber(codon);
			annotationCoding.setCodonPosition((byte) (module+1));
			annotationCoding.setAminoacidChange(change);
			annotationsCoding.add(annotationCoding);
		}
		return annotationsCoding;
	}
	private VariantFunctionalAnnotation makeAnnotationClose(GenomicVariant variant, Transcript t, int offsetUpstream, int offsetDownstream) {
		VariantFunctionalAnnotation annotationClose = null;
		if (variant.getLast()<t.getFirst() ) {
			//Variant before the transcript in genomic location. Check upstream for positive strand and downstream for negative strand
			if(t.isPositiveStrand() && variant.getLast()>=t.getFirst()-offsetUpstream) {
				annotationClose = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_UPSTREAM);
			} else if(t.isNegativeStrand() && variant.getLast()>=t.getFirst()-offsetDownstream) {
				annotationClose = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_DOWNSTREAM);
			}
		} else if (variant.getFirst()>t.getLast()) {
			//Variant after the transcript in genomic location. Check upstream for negative strand and downstream for positive strand
			if(t.isPositiveStrand() && variant.getFirst()<=t.getLast()+offsetDownstream) {
				annotationClose = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_DOWNSTREAM);
			}
			else if(t.isNegativeStrand() && variant.getFirst()<=t.getLast()+offsetUpstream) {
				annotationClose = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_UPSTREAM);
			}
		}
		if(annotationClose!=null) annotationClose.setTranscript(t);
		return annotationClose;
	}
	private VariantFunctionalAnnotation makeIntronAnnotation(GenomicVariant variant, Transcript t, VariantAnnotationParameters parameters) {
		VariantFunctionalAnnotation annotationIntron;
		int spliceRegionIntronOffSet = parameters.getSpliceRegionIntronOffset();
		TranscriptSegment closeSegmentLeft = t.getTranscriptSegmentByAbsolutePosition(variant.getFirst()-spliceRegionIntronOffSet);
		TranscriptSegment closeSegmentRight = t.getTranscriptSegmentByAbsolutePosition(variant.getLast()+spliceRegionIntronOffSet);
		if(closeSegmentLeft!=null) {
			int distance = variant.getFirst()-closeSegmentLeft.getLast(); 
			if(t.isNegativeStrand() && distance<=parameters.getSpliceAcceptorOffset()) annotationIntron = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_SPLICE_ACCEPTOR);
			else if (!t.isNegativeStrand() && distance<=parameters.getSpliceDonorOffset()) annotationIntron = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_SPLICE_DONOR);
			else annotationIntron = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_SPLICE_REGION);
		} else if (closeSegmentRight!=null) {
			int distance = closeSegmentRight.getFirst()-variant.getLast(); 
			if (t.isNegativeStrand() && distance<=parameters.getSpliceDonorOffset()) annotationIntron = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_SPLICE_DONOR);
			else if (!t.isNegativeStrand() && distance<=parameters.getSpliceAcceptorOffset()) annotationIntron = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_SPLICE_ACCEPTOR);
			else annotationIntron = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_SPLICE_REGION);
		} else {
			annotationIntron = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_INTRON);
		}
		annotationIntron.setTranscript(t);
		return annotationIntron;
	}
	public void fillSequenceTranscripts(ReferenceGenome genome) {
		if(sortedTranscripts==null)System.err.println("Null sorted transcripts");
		for(Transcript t:sortedTranscripts) {
			StringBuilder transcriptSeq = new StringBuilder();
			boolean segmentNotFound = false;
			List<TranscriptSegment> segments = new ArrayList<TranscriptSegment>(t.getTranscriptSegments());
			if(t.isNegativeStrand()) Collections.reverse(segments);
			for (TranscriptSegment e: segments) {
				CharSequence genomicSeq = genome.getReference(t.getSequenceName(),e.getFirst(),e.getLast());
				if(genomicSeq == null) {
					System.err.println("WARN: Transcript segment at genomic location "+t.getSequenceName()+":"+e.getFirst()+"-"+ e.getLast()+" not found for transcript: "+t.getId());
					segmentNotFound = true;
					break;
				}
				if(t.isNegativeStrand()) {
					genomicSeq = DNAMaskedSequence.getReverseComplement(genomicSeq);
				}
				transcriptSeq.append(genomicSeq);
			}
			if(!segmentNotFound) t.setCDNASequence(new DNAMaskedSequence(transcriptSeq.toString().toUpperCase()));
		}
	}
}
