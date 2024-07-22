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
package ngsep.transcriptome.io;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.genome.io.GFF3GenomicFeature;
import ngsep.genome.io.GFF3GenomicFeatureLine;
import ngsep.genome.io.GFF3Loader;
import ngsep.main.io.ConcatGZIPInputStream;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.transcriptome.TranscriptSegment;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.Polypeptide;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.Transcriptome;

public class GFF3TranscriptomeHandler {
	private QualifiedSequenceList sequenceNames;
	
	private boolean loadTextAnnotations = false;
	
	private Logger log = Logger.getLogger(GFF3TranscriptomeHandler.class.getName());
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	
	public boolean isLoadTextAnnotations() {
		return loadTextAnnotations;
	}
	public void setLoadTextAnnotations(boolean loadTextAnnotations) {
		this.loadTextAnnotations = loadTextAnnotations;
	}
	public GFF3TranscriptomeHandler () {
		sequenceNames = new QualifiedSequenceList();
	}
	
	public GFF3TranscriptomeHandler (QualifiedSequenceList sequenceNames) {
		setSequenceNames(sequenceNames);
	}
	public QualifiedSequenceList getSequenceNames() {
		return sequenceNames;
	}
	public void setSequenceNames(QualifiedSequenceList sequenceNames) {
		if(sequenceNames==null) this.sequenceNames = new QualifiedSequenceList();
		else this.sequenceNames = sequenceNames;
	}
	/**
	 * Loads the genes, transcripts and exons information from the given map file
	 * @param filename Name of the file with the transcriptome description
	 * @throws IOException If the file can not be read
	 */
	public Transcriptome loadMap(String filename) throws IOException {
		InputStream is = null;
		try {
			is = new FileInputStream(filename);
			if(filename.toLowerCase().endsWith(".gz")) {
				is = new ConcatGZIPInputStream(is);
			}
			return loadMap(is);
		} finally {
			if(is!=null) is.close();
		}
		
	}
	/**
	 * Loads the genes, transcripts and exons information from the given map file
	 * @param is Stream with the transcriptome description in GFF3 format
	 * @throws IOException If the file can not be read
	 */
	public Transcriptome loadMap(InputStream is) throws IOException {
		Transcriptome answer = new Transcriptome(sequenceNames);
		Map<String,GFF3GenomicFeature> featuresWithId = new HashMap<String, GFF3GenomicFeature>();
		List<GFF3GenomicFeatureLine> featureLinesWithParentAndNoId = new ArrayList<>();
		GFF3Loader loader = new GFF3Loader(sequenceNames);
		loader.setLog(log);
		List<GFF3GenomicFeatureLine> allLines = loader.loadFeatureLines(is);
		for(GFF3GenomicFeatureLine line:allLines) {
			String id = line.getId();
			if(id!=null) {
				GFF3GenomicFeature feature = featuresWithId.get(id);
				if(feature==null) {
					feature = new GFF3GenomicFeature(line);
					featuresWithId.put(id,feature);
				} else {
					feature.addLine(line);
				}
			} else if(line.hasParents()) {
				featureLinesWithParentAndNoId.add(line);
			}
		}
		List<GFF3GenomicFeatureLine> cdsFeaturesWithoutParent = new ArrayList<>();
		//Assign parents for features with IDs and collect CDS without parents
		for(GFF3GenomicFeature feature:featuresWithId.values()) {
			Set<String> parentIds = feature.getParentIds();
			boolean parentPresent = parentIds.size()>0;
			for(String parentId:parentIds) {
				GFF3GenomicFeature parent = featuresWithId.get(parentId);
				if(parent==null) {
					log.warning("Parent "+parentId+" not found for feature id: "+feature.getId());
					continue;
				}
				parent.addChild(feature);
			}
			if(!parentPresent) {
				String derivesId = feature.getAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_DERIVES_FROM);
				parentPresent = derivesId !=null; 
				if(parentPresent) {
					GFF3GenomicFeature parent = featuresWithId.get(derivesId);
					if(parent==null) {
						log.warning("Parent "+derivesId+" not found for feature id: "+feature.getId());
						continue;
					}
					parent.addChild(feature);
				}
			}
			if(!parentPresent && GFF3GenomicFeatureLine.FEATURE_TYPE_CDS.equals(feature.getType()) && feature.getLines().size()==1) cdsFeaturesWithoutParent.add(feature.getLines().get(0));
		}
		
		//Assign parent for features with parent id and no id
		for(GFF3GenomicFeatureLine featureLine:featureLinesWithParentAndNoId) {
			Set<String> parentIds = featureLine.getParentIds();
			for(String parentId:parentIds) {
				GFF3GenomicFeature parent = featuresWithId.get(parentId);
				if(parent==null) {
					log.warning("Parent  "+parentId+" not found for feature line "+featureLine.getLineNumber()+" at "+featureLine.getSequenceName()+":"+featureLine.getFirst()+"-"+featureLine.getLast()+" Id: "+featureLine.getId());
					continue;
				}
				parent.addChildWithoutId(featureLine);
			}
		}
		int numGenes = 0;
		int numTranscripts = 0;
		for(GFF3GenomicFeature feature:featuresWithId.values()) {
			if(GFF3GenomicFeatureLine.FEATURE_TYPE_GENE.equals(feature.getType()) || GFF3GenomicFeatureLine.FEATURE_TYPE_TRGENE.equals(feature.getType())|| GFF3GenomicFeatureLine.FEATURE_TYPE_PCGENE.equals(feature.getType()) || GFF3GenomicFeatureLine.FEATURE_TYPE_PSEUDOGENE.equals(feature.getType()) ) {
				Gene gene = createGeneFromGeneFeature(feature);
				gene.setOntologyTerms(feature.getOntologyTerms());
				gene.setDatabaseReferences(feature.getDatabaseReferences());
				String note = feature.getAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_NOTE);
				if(loadTextAnnotations && note!=null) gene.addTextFunctionalAnnotation(note);
				boolean transcriptsAdded = false;
				for(GFF3GenomicFeature geneChild:feature.getChildren()) {
					if(loadTextAnnotations) gene.addTextFunctionalAnnotations(geneChild.getProducts());
					Transcript transcript = null;
					List<TranscriptSegment> segments = new ArrayList<TranscriptSegment>();
					if(GFF3GenomicFeatureLine.FEATURE_TYPE_MRNA.equals(geneChild.getType())) {
						transcript = createTranscriptFromMRNAFeature(geneChild);
						if(transcript == null) continue;
						for(GFF3GenomicFeature mrnaChild:geneChild.getChildren()) {
							if(GFF3GenomicFeatureLine.FEATURE_TYPE_POLYPEPTIDE.equals(mrnaChild.getType())) {
								if(loadTextAnnotations) {
									loadPolypeptide(transcript, mrnaChild);
								}
							} else {
								//CDS
								segments.addAll(createFeatureSegments(mrnaChild, transcript));
							}
						}
						if(segments.size()==0) {
							log.warning("CDS not found for mRNA: "+transcript.getId()+". Transcript not loaded");
							continue;
						}
					} else if (GFF3GenomicFeatureLine.FEATURE_TYPE_CDS.equals(geneChild.getType())) {
						//Direct CDs without transcript
						transcript = createTranscriptFromCDSFeature(geneChild);
						if(transcript == null) continue;
						if(transcript.getId()==null) transcript.setId(gene.getId()+"_mRNA");
						segments.addAll(createFeatureSegments(geneChild, transcript));
					}
					if(transcript!=null) {
						transcript.setGene(gene);
						transcript.setTranscriptSegments(segments);
						answer.addTranscript(transcript);
						numTranscripts++;
						transcriptsAdded = true;
					} else {
						
					}
				}
				if(transcriptsAdded) numGenes++;
				else log.warning("No transcripts could be loaded for gene: "+gene.getId()+". Gene not loaded"); 
			}
		}
		//Add orphan CDS
		for (GFF3GenomicFeatureLine cdsFeature: cdsFeaturesWithoutParent) {
			String id = cdsFeature.getId();
			
			Gene gene = new Gene(id+"_gene", id, cdsFeature.getSequenceName(), cdsFeature.getFirst(), cdsFeature.getLast(), cdsFeature.isNegativeStrand() );
			gene.setOntologyTerms(cdsFeature.getOntologyTerms());
			Transcript transcript = new Transcript(id+"_mRNA", cdsFeature.getSequenceName(), cdsFeature.getFirst(), cdsFeature.getLast(), cdsFeature.isNegativeStrand() );
			transcript.setGene(gene);
			List<TranscriptSegment> segments = new ArrayList<TranscriptSegment>(1);
			TranscriptSegment segment = new TranscriptSegment(transcript, cdsFeature.getFirst(), cdsFeature.getLast());
			segment.setStatus(TranscriptSegment.STATUS_CODING);
			segment.setFirstCodonPositionOffset(cdsFeature.getPhase());
			segments.add(segment);
			transcript.setTranscriptSegments(segments);
			answer.addTranscript(transcript);
			numGenes++;
			numTranscripts++;
		}
		
		if(numGenes==0) log.warning("No genes found in transcriptome file. Check that gene features have \""+GFF3GenomicFeatureLine.FEATURE_TYPE_GENE+"\" as the exact feature type");
		if(numTranscripts==0) log.warning("No transcripts found in transcriptome file. Check that transcript features have \""+GFF3GenomicFeatureLine.FEATURE_TYPE_MRNA+"\" as the exact feature type");
		//log.warning("Loaded "+numTranscripts+" transcripts for "+numGenes+" genes");
		return answer;
	}
	private void loadPolypeptide(Transcript transcript, GFF3GenomicFeature feature) {
		Polypeptide p = new Polypeptide(feature.getId(), transcript);
		p.setOntologyTerms(feature.getOntologyTerms());
		p.setProducts(feature.getProducts());
		for(GFF3GenomicFeature child:feature.getChildren()) {
			if(GFF3GenomicFeatureLine.FEATURE_TYPE_PROTEIN_MATCH.equals(child.getType())) {
				p.addPfamTerm(child.getName(), child.getAnnotation(GFF3GenomicFeatureLine.ATTRIBUTE_SIGNATURE_DESC));
			}
		}
		transcript.addPolypeptide(p);
	}
	private Transcript createTranscriptFromCDSFeature(GFF3GenomicFeature feature) {
		String sequenceName = null;
		boolean negativeStrand = false;
		int first=-1;
		int last= -1;
		for(GFF3GenomicFeatureLine line:feature.getLines()) {
			if(sequenceName==null) {
				sequenceName = line.getSequenceName();
				negativeStrand = line.isNegativeStrand();
				first = line.getFirst();
				last = line.getLast();
			} else if (!sequenceName.equals(line.getSequenceName())) {
				log.warning("Can not load CDS with ID "+feature.getId()+". Multiple lines of a CDS feature must have the same sequence name. Expected: "+sequenceName+" given: "+line.getSequenceName());
				return null;
			} else if(negativeStrand!=line.isNegativeStrand()) {
				log.warning("Can not load CDS with ID "+feature.getId()+". Multiple lines of a CDS feature must have the same strand");
				return null;
			} else {
				first = Math.min(first, line.getFirst());
				last = Math.max(last, line.getLast());
			}
		}
		if(sequenceName==null) throw new RuntimeException("Error loading CDS feature with ID: "+feature.getId()+" at  feature lines: "+feature.getLines().size());
		return new Transcript(feature.getId(), sequenceName, first, last, negativeStrand);
	}
	private Gene createGeneFromGeneFeature(GFF3GenomicFeature feature) {
		List<GFF3GenomicFeatureLine> geneLines = feature.getLines();
		if(geneLines.size()>1) {
			//log.warning("Multiple lines found for feature with ID "+feature.getId()+". Features of type "+FEATURE_TYPE_MRNA+" should not have multiple lines. Processing first line only");
			for(GFF3GenomicFeatureLine line: geneLines) log.warning("Multiple lines found for feature with ID "+feature.getId()+". Features of type "+feature.getType()+" should not have multiple lines. Feature at "+line.getSequenceName()+":"+line.getFirst()+"-"+line.getLast()+" Id: "+line.getId()+" Number: "+line.getLineNumber());
			throw new RuntimeException();
		}
		GFF3GenomicFeatureLine geneLine = geneLines.get(0);
		String name = feature.getName();
		if(name ==null) name = feature.getAnnotation("name");
		Gene gene = new Gene(feature.getId(), name, geneLine.getSequenceName(), geneLine.getFirst(), geneLine.getLast(), geneLine.isNegativeStrand() );
		return gene;
	}
	private Transcript createTranscriptFromMRNAFeature(GFF3GenomicFeature feature) {
		List<GFF3GenomicFeatureLine> mrnaLines = feature.getLines();
		if(mrnaLines.size()>1) {
			//log.warning("Multiple lines found for feature with ID "+feature.getId()+". Features of type "+FEATURE_TYPE_MRNA+" should not have multiple lines. Processing first line only");
			for(GFF3GenomicFeatureLine line: mrnaLines) log.warning("Multiple lines found for feature with ID "+feature.getId()+". Features of type "+feature.getType()+" should not have multiple lines. Feature at "+line.getSequenceName()+":"+line.getFirst()+"-"+line.getLast()+" Id: "+line.getId()+" Number: "+line.getLineNumber());
			throw new RuntimeException();
		}
		GFF3GenomicFeatureLine mrnaLine = mrnaLines.get(0);
		Transcript transcript = new Transcript(feature.getId(), mrnaLine.getSequenceName(), mrnaLine.getFirst(), mrnaLine.getLast(), mrnaLine.isNegativeStrand() );
		return transcript;
	}
	private List<TranscriptSegment> createFeatureSegments(GFF3GenomicFeature feature, Transcript transcript) {
		List<TranscriptSegment> segments = new ArrayList<>();
		String featureType = feature.getType();
		for(GFF3GenomicFeatureLine line:feature.getLines()) {
			TranscriptSegment segment = new TranscriptSegment(transcript, line.getFirst(), line.getLast());
			if (GFF3GenomicFeatureLine.FEATURE_TYPE_5PUTR.equals(featureType)) {
				segment.setStatus(TranscriptSegment.STATUS_5P_UTR);
				segments.add(segment);
			} else if (GFF3GenomicFeatureLine.FEATURE_TYPE_3PUTR.equals(featureType)) {
				segment.setStatus(TranscriptSegment.STATUS_3P_UTR);
				segments.add(segment);
			} else if (GFF3GenomicFeatureLine.FEATURE_TYPE_CDS.equals(featureType)) {
				segment.setStatus(TranscriptSegment.STATUS_CODING);
				segment.setFirstCodonPositionOffset(line.getPhase());
				segments.add(segment);
			}
		}
		return segments;
	}

	/**
	 * Loads the transcripts sequences from the given fasta file 
	 * @param filename Name of the fasta file with the transcript sequences
	 * @throws IOException If the file can not be read
	 */
	public void loadSequences(Transcriptome transcriptome, String filename) throws IOException {
		FastaSequencesHandler sequencesHandler = new FastaSequencesHandler();
		sequencesHandler.setSequenceType(DNAMaskedSequence.class);
		List<QualifiedSequence> transcriptSeqs = sequencesHandler.loadSequences(filename);
		
		int n = transcriptSeqs.size();
		for(int i=0;i<n;i++ ) {
			QualifiedSequence sequence = transcriptSeqs.get(i);
			Transcript transcript = transcriptome.getTranscript(sequence.getName());
			if(transcript!=null) {
				transcript.setCDNASequence((DNAMaskedSequence)sequence.getCharacters());
			} else {
				log.warning("Transcript not found with id: "+sequence.getName());
			}
		}
	}
}