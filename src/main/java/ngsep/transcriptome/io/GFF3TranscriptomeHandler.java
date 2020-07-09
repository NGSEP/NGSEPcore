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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.main.io.ConcatGZIPInputStream;
import ngsep.main.io.ParseUtils;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.transcriptome.TranscriptSegment;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.Transcriptome;

public class GFF3TranscriptomeHandler {
	public static final String FEATURE_TYPE_GENE = "gene";
	public static final String FEATURE_TYPE_TRGENE = "transposable_element_gene";
	public static final String FEATURE_TYPE_EXON = "exon";
	public static final String FEATURE_TYPE_MRNA = "mRNA";
	public static final String FEATURE_TYPE_5PUTR = "five_prime_UTR";
	public static final String FEATURE_TYPE_3PUTR = "three_prime_UTR";
	public static final String FEATURE_TYPE_CDS = "CDS";
	
	//Predefined attributes according to the gff3 specification
	public static final String ATTRIBUTE_ID = "ID";
	public static final String ATTRIBUTE_NAME = "Name";
	public static final String ATTRIBUTE_ALIAS = "Alias";
	public static final String ATTRIBUTE_PARENT = "Parent";
	public static final String ATTRIBUTE_TARGET = "Target";
	public static final String ATTRIBUTE_GAP = "Gap";
	public static final String ATTRIBUTE_DERIVES_FROM = "Derives_from";
	public static final String ATTRIBUTE_NOTE = "Note";
	public static final String ATTRIBUTE_DBXREF = "Dbxref";
	public static final String ATTRIBUTE_ONTOLOGY = "Ontology_term";
	public static final String ATTRIBUTE_CIRCULAR = "Is_circular";
	
	
	
	
	
	
	public static final String [] supportedFeatureTypes = {FEATURE_TYPE_GENE,FEATURE_TYPE_CDS,FEATURE_TYPE_MRNA,FEATURE_TYPE_TRGENE,FEATURE_TYPE_5PUTR,FEATURE_TYPE_3PUTR};
	//TODO: Use a file resource
	public static final String [] supportedFeatureTypesSOFAIDs = {"SO:0000704","SO:0000316","SO:0000234","SO:0000111","SO:0000204","SO:0000205"};
	private QualifiedSequenceList sequenceNames;
	
	private Logger log = Logger.getLogger(GFF3TranscriptomeHandler.class.getName());
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
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
		try (BufferedReader in = new BufferedReader(new InputStreamReader(is))) {
			String line=in.readLine();
			if(!line.startsWith("##gff")) throw new IOException("File does not have GFF format");
			line=in.readLine();
			for(int i=1;line!=null;i++) {
				if("##FASTA".equals(line.trim())) break;
				if(line.length()>0 && line.charAt(0)!='#') {
					GFF3GenomicFeatureLine featureLine;
					try {
						featureLine = loadFeatureLine(line,i); 
					} catch (RuntimeException e) {
						e.printStackTrace();
						log.warning("Error loading line "+line+". "+e.getMessage());
						line=in.readLine();
						continue;
					}
					if(featureLine==null) {
						line=in.readLine();
						continue;
					}
					String id = featureLine.getId();
					if(id!=null) {
						GFF3GenomicFeature feature = featuresWithId.get(id);
						if(feature==null) {
							feature = new GFF3GenomicFeature(featureLine);
							featuresWithId.put(id,feature);
						} else {
							feature.addLine(featureLine);
						}
					} else if(featureLine.hasParents()) {
						featureLinesWithParentAndNoId.add(featureLine);
					}
				}
				line=in.readLine();
			}
		}
		//Assign parents for features with IDs
		for(GFF3GenomicFeature feature:featuresWithId.values()) {
			Set<String> parentIds = feature.getParentIds();
			for(String parentId:parentIds) {
				GFF3GenomicFeature parent = featuresWithId.get(parentId);
				if(parent==null) {
					log.warning("Parent "+parentId+" not found for feature id: "+feature.getId());
					continue;
				}
				parent.addChild(feature);
			}
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
			if(FEATURE_TYPE_GENE.equals(feature.getType()) || FEATURE_TYPE_TRGENE.equals(feature.getType())) {
				Gene gene = createGeneFromMRNAFeature(feature);
				numGenes++;
				gene.setOntologyTerms(feature.getOntologyTerms());
				gene.setDatabaseReferences(feature.getDatabaseReferences());
				
				for(GFF3GenomicFeature geneChild:feature.getChildren()) {
					Transcript transcript = null;
					List<TranscriptSegment> segments = new ArrayList<TranscriptSegment>();
					if(FEATURE_TYPE_MRNA.equals(geneChild.getType())) {
						transcript = createTranscriptFromMRNAFeature(geneChild);
						if(transcript == null) continue;
						for(GFF3GenomicFeature mrnaChild:geneChild.getChildren()) {
							segments.addAll(createFeatureSegments(mrnaChild, transcript));
						}
					} else if (FEATURE_TYPE_CDS.equals(geneChild.getType())) {
						//Direct CDs without transcript
						transcript = createTranscriptFromCDSFeature(geneChild);
						if(transcript == null) continue;
						if(transcript.getId()==null) transcript.setId(gene.getId());
						segments.addAll(createFeatureSegments(geneChild, transcript));
					}
					if(transcript!=null) {
						transcript.setGene(gene);
						transcript.setTranscriptSegments(segments);
						answer.addTranscript(transcript);
						numTranscripts++;
					}
					
				}
			}
		}
		if(numGenes==0) log.warning("No genes found in transcriptome file. Check that gene features have \""+FEATURE_TYPE_GENE+"\" as the exact feature type");
		if(numTranscripts==0) log.warning("No transcripts found in transcriptome file. Check that transcript features have \""+FEATURE_TYPE_MRNA+"\" as the exact feature type");
		//log.warning("Loaded "+numTranscripts+" transcripts for "+numGenes+" genes");
		return answer;
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
	private Gene createGeneFromMRNAFeature(GFF3GenomicFeature feature) {
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
			if (FEATURE_TYPE_5PUTR.equals(featureType)) {
				segment.setStatus(TranscriptSegment.STATUS_5P_UTR);
				segments.add(segment);
			} else if (FEATURE_TYPE_3PUTR.equals(featureType)) {
				segment.setStatus(TranscriptSegment.STATUS_3P_UTR);
				segments.add(segment);
			} else if (FEATURE_TYPE_CDS.equals(featureType)) {
				segment.setStatus(TranscriptSegment.STATUS_CODING);
				segment.setFirstCodonPositionOffset(line.getPhase());
				segments.add(segment);
			}
		}
		return segments;
	}
	public GFF3GenomicFeatureLine loadFeatureLine (String line, int lineNumber) {
		String [] items = ParseUtils.parseString(line, '\t');
		//Memory saver for sequenceNames
		QualifiedSequence seq;
		try {
			seq = sequenceNames.addOrLookupName(items[0]);
		} catch (Exception e2) {
			log.warning("Can not load genomic feature at line: "+line+". Unrecognized sequence name. "+items[0]);
			return null;
		}
		String type = loadType(items[2]);
		int first;
		try {
			first = Integer.parseInt(items[3]);
		} catch (NumberFormatException e1) {
			log.warning("Error loading feature at line "+line+". Start coordinate must be a positive integer");
			return null;
		}
		int last;
		try {
			last = Integer.parseInt(items[4]);
		} catch (NumberFormatException e1) {
			log.warning("Error loading feature at line "+line+". End coordinate must be a positive integer");
			return null;
		}
		GFF3GenomicFeatureLine answer = new GFF3GenomicFeatureLine(seq.getName(), first, last,type);
		answer.setSource(items[1]);
		answer.setLineNumber(lineNumber);
		answer.setNegativeStrand(items[6].charAt(0)=='-');
		if(items[7].charAt(0)!='.') {
			try {
				answer.setPhase(Byte.parseByte(items[7]));
			} catch (NumberFormatException e) {
				log.warning("Error loading phase of feature at "+answer.getSequenceName()+":"+answer.getFirst()+"-"+answer.getLast()+". Value: "+items[7]);
			}
		} else if (FEATURE_TYPE_CDS.equals(type)) {
			log.warning("Error loading phase of feature at "+answer.getSequenceName()+":"+answer.getFirst()+"-"+answer.getLast()+". The phase field is required for CDS features");
		}
		String errorMessage = "Error loading annotations of feature at "+answer.getSequenceName()+":"+answer.getFirst()+"-"+answer.getLast()+". ";
		String [] annItems = ParseUtils.parseStringWithText(items[8], ';', '\"');
		for(int i=0;i<annItems.length;i++) {
			if(annItems[i].trim().length()==0) continue;
			int idx = annItems[i].indexOf("=");
			if(idx <=0 || idx == annItems[i].length()-1) log.warning(errorMessage+"Annotation "+annItems[i]+" is not a key-value pair");
			else answer.addAnnotation(annItems[i].substring(0,idx),annItems[i].substring(idx+1));
		}
		return answer;
	}
	/**
	 * Returns the known feature type given the type string. This helps to save memory with feature types
	 * @param typeStr Type to search
	 * @return String If typeStr is a supported SOFA id, it returns the correspoonding name.
	 * otherwise it returns a supported type equals to typeStr or typeStr itself if it is not a known type
	 */
	public String loadType(String typeStr) {
		String type = typeStr;
		int idxType = -1;
		//Search type name
		for(int i=0;i<supportedFeatureTypes.length && idxType == -1;i++) {
			if(supportedFeatureTypes[i].equals(typeStr)) {
				idxType = i;
				type = supportedFeatureTypes[idxType];
			}
		}
		//Search type id
		for(int i=0;i<supportedFeatureTypesSOFAIDs.length && idxType == -1;i++) {
			if(supportedFeatureTypesSOFAIDs[i].equals(typeStr)) {
				idxType = i;
				type = supportedFeatureTypes[idxType];
			}
		}
		return type;
	}
	/**
	 * Loads the transcripts sequences from the given fasta file 
	 * @param filename Name of the fasta file with the transcript sequences
	 * @throws IOException If the file can not be read
	 */
	public void loadSequences(Transcriptome transcriptome, String filename) throws IOException {
		FastaSequencesHandler sequencesHandler = new FastaSequencesHandler();
		sequencesHandler.setSequenceType(DNAMaskedSequence.class);
		QualifiedSequenceList transcriptSeqs = sequencesHandler.loadSequences(filename);
		
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
class GFF3GenomicFeature {
	private List<GFF3GenomicFeatureLine> lines = new ArrayList<>();
	private List<GFF3GenomicFeature> children = new ArrayList<GFF3GenomicFeature>();
	private String id;
	private String type;
	
	public GFF3GenomicFeature(GFF3GenomicFeatureLine line) {
		id = line.getId();
		type = line.getType();
		lines.add(line);
	}
	
	/**
	 * Adds the given line to the feature
	 * @param line New line to add
	 * @throws IOException If the ID of the new line is different from the ID of this feature
	 * @throws IOException If the type of the new line is different from the type of this feature
	 */
	public void addLine(GFF3GenomicFeatureLine line) throws IOException {
		if(id==null) throw new IOException("Error loading line at "+line.getSequenceName()+":"+line.getFirst()+"-"+line.getLast()+". Features with null ids can not have multiple lines");
		if(!id.equals(line.getId())) throw new IOException("Error loading line at "+line.getSequenceName()+":"+line.getFirst()+"-"+line.getLast()+" the type must be the same for all lines with the same id. Expected: "+type+" given: "+line.getType());;
		if(!type.equals(line.getType())) throw new IOException("Error loading line at "+line.getSequenceName()+":"+line.getFirst()+"-"+line.getLast()+" the type must be the same for all lines with the same id. Expected: "+type+" given: "+line.getType());;
		lines.add(line);
	}
	public List<GFF3GenomicFeatureLine> getLines() {
		return lines;
	}
	public void addChild(GFF3GenomicFeature child) {
		children.add(child);
	}
	public void addChildWithoutId(GFF3GenomicFeatureLine featureLine) {
		children.add(new GFF3GenomicFeature(featureLine)); 
		
	}
	public List<GFF3GenomicFeature> getChildren() {
		return children;
	}
	public String getType() {
		return type;
	}
	public String getId() {
		return id;
	}
	public String getName() {
		for(GFF3GenomicFeatureLine line:lines) {
			String name = line.getName(); 
			if(name!=null) return name; 
		}
		return null;
	}
	public String getAnnotation(String attribute) {
		for(GFF3GenomicFeatureLine line:lines) {
			String ann = line.getAnnotation(attribute); 
			if(ann!=null) return ann; 
		}
		return null;
	}
	public boolean hasParents() {
		for(GFF3GenomicFeatureLine line:lines) {
			if(line.hasParents()) return true; 
		}
		return false;
	}
	public Set<String> getParentIds() {
		Set<String> ids = new HashSet<>();
		for(GFF3GenomicFeatureLine line:lines) {
			ids.addAll(line.getParentIds());
		}
		return ids;
	}
	public List<String> getOntologyTerms() {
		return getValuesList(GFF3TranscriptomeHandler.ATTRIBUTE_ONTOLOGY);
	}
	
	public List<String> getDatabaseReferences() {
		return getValuesList(GFF3TranscriptomeHandler.ATTRIBUTE_DBXREF);
	}

	public List<String> getValuesList(String attribute) {
		List<String> values = new ArrayList<>();
		for(GFF3GenomicFeatureLine line:lines) {
			values.addAll(line.getValuesList(attribute));
		}
		return values;
	}
}
class GFF3GenomicFeatureLine implements GenomicRegion {
	private String sequenceName;
	private String source;
	private String type;
	private int first;
	private int last;
	private boolean negativeStrand=false;
	private byte phase = -1;
	private int lineNumber = 0;
	private Map<String,String> annotations = new HashMap<String, String>();
	
	public GFF3GenomicFeatureLine(String sequenceName, int first, int last, String type) {
		super();
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		this.type = type;
	}
	
	public String getSequenceName() {
		return sequenceName;
	}
	
	public int getFirst() {
		return first;
	}
	
	public int getLast() {
		return last;
	}
	
	public String getType() {
		return type;
	}
	
	public String getSource() {
		return source;
	}
	
	public boolean isNegativeStrand() {
		return negativeStrand;
	}
	public Map<String, String> getAnnotations() {
		return annotations;
	}
	public void addAnnotation(String key, String value) {
		annotations.put(key, value);
	}
	@Override
	public int length() {
		return last-first+1;
	}
	@Override
	public boolean isPositiveStrand() {
		return !negativeStrand;
	}
	public String getAnnotation(String tag) {
		return annotations.get(tag);
	}
	public void setSource(String source) {
		this.source = source;
	}
	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}
	public byte getPhase() {
		return phase;
	}
	public void setPhase(byte phase) {
		this.phase = phase;
	}
	public String getId() {
		return annotations.get(GFF3TranscriptomeHandler.ATTRIBUTE_ID);
	}
	public String getName() {
		return annotations.get(GFF3TranscriptomeHandler.ATTRIBUTE_NAME);
	}
	public Set<String> getParentIds() {
		String parentsStr = annotations.get(GFF3TranscriptomeHandler.ATTRIBUTE_PARENT);
		if(parentsStr==null) return new HashSet<>();
		String [] parentIdsArr = parentsStr.split(",");
		Set<String> answer = new HashSet<>();
		for(int i = 0;i<parentIdsArr.length;i++) {
			answer.add(parentIdsArr[i]);
		}
		return answer;
	}
	public boolean hasParents() {
		return annotations.get(GFF3TranscriptomeHandler.ATTRIBUTE_PARENT)!=null;
	}
	public List<String> getOntologyTerms() {
		return getValuesList(GFF3TranscriptomeHandler.ATTRIBUTE_ONTOLOGY);
	}
	public List<String> getDatabaseReferences() {
		return getValuesList(GFF3TranscriptomeHandler.ATTRIBUTE_DBXREF);
	}
	public List<String> getValuesList(String attribute) {
		String parentsStr = annotations.get(attribute);
		if(parentsStr==null) return new ArrayList<>();
		String [] parentIdsArr = parentsStr.split(",");
		return Arrays.asList(parentIdsArr);
	}

	public int getLineNumber() {
		return lineNumber;
	}

	public void setLineNumber(int lineNumber) {
		this.lineNumber = lineNumber;
	}
	
	
}