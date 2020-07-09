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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
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
	
	public List<Gene> getAllGenes () {
		return new ArrayList<>(genesMap.values());
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
				//System.err.println("Segment: "+segmentStart.getSequenceName()+":"+segmentStart.getFirst()+"-"+segmentStart.getLast()+" coding: "+segmentStart.isCoding()+" status "+segmentStart.getStatus());
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
		int transcriptionEnd = t.getCodingRelativeEnd();
		int absoluteFirst = variant.getFirst();
		if(t.isNegativeStrand()) {
			absoluteFirst = variant.getLast();
		}
		int varTranscriptStart = t.getRelativeTranscriptPosition(absoluteFirst);
		int varTranscriptEnd = varTranscriptStart+(variant.getLast()-variant.getFirst())+1;
		int varCodingStart = varTranscriptStart-transcriptionStart;
		if(varCodingStart<0) {
			//Weird case of variant just before transcription start
			VariantFunctionalAnnotation annotation = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_5P_UTR);
			annotation.setTranscript(t);
			annotationsCoding.add(annotation);
			return annotationsCoding;
		}
		if(varTranscriptStart>transcriptionEnd) {
			VariantFunctionalAnnotation annotation = new VariantFunctionalAnnotation(variant, VariantFunctionalAnnotationType.ANNOTATION_3P_UTR);
			annotation.setTranscript(t);
			annotationsCoding.add(annotation);
			return annotationsCoding;
		}
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
			String alternativeR = DNAMaskedSequence.getReverseComplement(alternative).toString();
			int differenceBases = alternative.length() - reference.length();
			int expectedProteinIncrease = differenceBases/3;
			String testA = alternative;
			if(t.isNegativeStrand()) testA = alternativeR;
			//System.err.println("Transcript: "+t.getId()+" Ref: "+testReference+". Start test: "+startTest+" module: "+module+" "+". Transcript start:"+varTranscriptStart+" coding start: "+varCodingStart);
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
	public void fillSequenceTranscripts(ReferenceGenome genome, Logger log) {
		if(sortedTranscripts==null) throw new RuntimeException("Null sorted transcripts");
		for(Transcript t:sortedTranscripts) {
			try {
				t.fillCDNASequence(genome);
			} catch (Exception e) {
				if(log!=null) log.warning(e.getMessage());
				else e.printStackTrace();
			}
		}
	}
	public void removeTranscript(Transcript t) {
		transcriptsMap.remove(t.getId());
		sortedTranscripts.remove(t);
		String geneId = t.getGeneId();
		List<Transcript> transcriptsGene = transcriptsByGene.get(geneId);
		if(transcriptsGene!=null) {
			transcriptsGene.remove(t);
			if(transcriptsGene.size()==0) { 
				transcriptsByGene.remove(geneId);
				genesMap.remove(geneId);
			}
		}
	}
	public void loadBulkTranscripts(Map<String,List<Transcript>> transcriptsBySeqName, String newGenesPrefix, Logger log ) {
		for(QualifiedSequence seq:sequenceNames) {
			List<Transcript> transcripts = transcriptsBySeqName.get(seq.getName());
			Collections.sort(transcripts,GenomicRegionPositionComparator.getInstance());
			List<Transcript> nextOverlappingCluster = new ArrayList<>();
			int lastCluster = 0;
			for(Transcript t:transcripts) {
				int cdsAbsStart = t.getCodingAbsoluteStart();
				int cdsAbsEnd = t.getCodingAbsoluteEnd();
				int cdsMin = Math.min(cdsAbsStart, cdsAbsEnd);
				int cdsMax = Math.max(cdsAbsStart, cdsAbsEnd);
				if(cdsAbsStart<t.getFirst() || cdsAbsStart>t.getLast() || cdsAbsEnd<t.getFirst() || cdsAbsEnd>t.getLast() ) {
					if(log!=null) log.warning("Improper CDS limits "+cdsAbsStart+"-"+cdsAbsEnd+" for transcript "+t.getId()+" at "+t.getSequenceName()+":"+t.getFirst()+"-"+t.getLast()+". Ignoring");
				}
				if(lastCluster <= cdsMin) {
					createOverlappingTranscripts(nextOverlappingCluster, newGenesPrefix, log);
					nextOverlappingCluster.clear();
				}
				nextOverlappingCluster.add(t);
				lastCluster = Math.max(lastCluster, cdsMax);
			}
			createOverlappingTranscripts(nextOverlappingCluster, newGenesPrefix, log);
		}
	}
	
	private void createOverlappingTranscripts(List<Transcript> transcripts, String newGenesPrefix, Logger log) {
		if(transcripts == null || transcripts.size()==0) {
			return;
		}
		Transcript first = transcripts.get(0);
		if(transcripts.size()==1) {
			addTranscriptWithNewGene(first, newGenesPrefix);
			if(log!=null) log.info("Added single transcript with gene id "+first.getGeneId());
			return;
		}
		//Find out direction and geneIds
		List<Transcript> positiveTranscripts = new ArrayList<>();
		List<Transcript> negativeTranscripts = new ArrayList<>();
		Set<String> positiveGeneIds = new HashSet<>();
		Set<String> negativeGeneIds = new HashSet<>();
		int positiveFirst = -1;
		int positiveLast = 0;
		int negativeFirst = -1;
		int negativeLast = 0;
		Transcript longestCDS = null;
		int maxCDSLen = 0;
		Set<String> cdsSequences = new HashSet<>();
		
		for(Transcript t: transcripts) {
			CharSequence cdsSeq = t.getCDSSequence();
			if(cdsSeq == null) continue;
			int cdsLength = cdsSeq.length();
			if(cdsLength == 0) continue;
			if(cdsSequences.contains(cdsSeq.toString())) continue;
			cdsSequences.add(cdsSeq.toString());
			if(maxCDSLen<cdsLength) {
				maxCDSLen = cdsLength;
				longestCDS = t;
			}
			if(t.isPositiveStrand()) {
				positiveTranscripts.add(t);
				positiveGeneIds.add(t.getGeneId());
				if(positiveFirst==-1 || positiveFirst>t.getFirst()) positiveFirst = t.getFirst();
				if(positiveLast<t.getLast()) positiveLast = t.getLast();
			} else {
				negativeTranscripts.add(t);
				negativeGeneIds.add(t.getGeneId());
				if(negativeFirst==-1 || negativeFirst>t.getFirst()) negativeFirst = t.getFirst();
				if(negativeLast<t.getLast()) negativeLast = t.getLast();
			}
		}
		if (negativeTranscripts.size()==0) {
			//Treat as single gene
			String geneId = null;
			if(positiveGeneIds.size()==1) geneId = positiveGeneIds.iterator().next();
			if(geneId==null || genesMap.containsKey(geneId)) {
				// Create new gene id
				geneId = createGeneId(newGenesPrefix);
			}
			Gene nextGene = new Gene(geneId,geneId,first.getSequenceName(),positiveFirst,positiveLast,false);
			for(Transcript t:positiveTranscripts) {
				t.setGene(nextGene);
				addTranscript(t);
			}
			if(log!=null) log.info("Added gene with "+positiveTranscripts.size()+" positive transcripts. Gene id "+geneId);
		} else if(positiveTranscripts.size()==0) {
			//Treat as single gene
			String geneId = null;
			if(negativeGeneIds.size()==1) geneId = negativeGeneIds.iterator().next();
			
			if(geneId==null || genesMap.containsKey(geneId)) {
				// Create new gene id
				geneId = createGeneId(newGenesPrefix);
				
			}
			Gene nextGene = new Gene(geneId,geneId,first.getSequenceName(),negativeFirst,negativeLast,true);
			for(Transcript t:negativeTranscripts) {
				t.setGene(nextGene);
				addTranscript(t);	
			}
			if(log!=null) log.info("Added gene with "+negativeTranscripts.size()+" negative transcripts. Gene id "+geneId);
		} else {
			addTranscriptWithNewGene(longestCDS, newGenesPrefix);
			if(log!=null) log.info("Overlapping positive and negative transcripts in "+longestCDS.getSequenceName()+":"+Math.min(positiveFirst, negativeFirst)+"-"+Math.max(positiveLast, negativeLast)+" choosing best transcript "+longestCDS.getId()+" from "+transcripts.size()+" transcripts. Added with gene id "+longestCDS.getGeneId());
		}
	}
	private void addTranscriptWithNewGene(Transcript transcript, String genesPrefix) {
		String geneId = transcript.getGeneId();
		if(genesMap.containsKey(geneId)) {
			// Create new gene id
			geneId = createGeneId(genesPrefix);
		}
		Gene gene = new Gene(geneId,geneId,transcript.getSequenceName(),transcript.getFirst(),transcript.getLast(),transcript.isNegativeStrand());
		transcript.setGene(gene);
		addTranscript(transcript);
	}
	private int newGeneId = 1;
	private String createGeneId(String genesPrefix) {
		String id = genesPrefix+newGeneId;
		while (genesMap.containsKey(id)) {
			newGeneId++;
			id = genesPrefix+newGeneId;
		}
		return id;
	}
}
