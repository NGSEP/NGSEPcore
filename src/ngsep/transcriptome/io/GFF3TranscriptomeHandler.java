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
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.main.io.ParseUtils;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.transcriptome.Exon;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.Transcriptome;

public class GFF3TranscriptomeHandler {
	public static final String FEATURE_TYPE_GENE = "gene";
	public static final String FEATURE_TYPE_TRGENE = "transposable_element_gene";
	public static final String FEATURE_TYPE_MRNA = "mRNA";
	public static final String FEATURE_TYPE_5PUTR = "five_prime_UTR";
	public static final String FEATURE_TYPE_5PUTR_INTRON = "five_prime_UTR_intron";
	public static final String FEATURE_TYPE_3PUTR = "three_prime_UTR";
	public static final String FEATURE_TYPE_CDS = "CDS";
	public static final String [] supportedFeatureTypes = {FEATURE_TYPE_CDS,FEATURE_TYPE_5PUTR,FEATURE_TYPE_5PUTR_INTRON,FEATURE_TYPE_GENE,FEATURE_TYPE_MRNA,FEATURE_TYPE_3PUTR,FEATURE_TYPE_TRGENE};
	private QualifiedSequenceList sequenceNames = new QualifiedSequenceList();
	
	private Logger log = Logger.getLogger(GFF3TranscriptomeHandler.class.getName());
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
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
		Transcriptome answer;
		if(sequenceNames.size()==0) answer = new Transcriptome();
		else answer = new Transcriptome(sequenceNames);
		Map<String,GFF3GenomicFeature> features = new TreeMap<String, GFF3GenomicFeature>();
		List<GFF3GenomicFeature> featuresWithParent = new ArrayList<GFF3GenomicFeature>();
		try (FileInputStream fis = new FileInputStream(filename);
			 BufferedReader in = new BufferedReader(new InputStreamReader(fis))) {
			String line=in.readLine();
			if(!line.startsWith("##gff")) throw new IOException("File "+filename+" does not have GFF format");
			line=in.readLine();
			while(line!=null) {
				if("##FASTA".equals(line.trim())) break;
				if(line.length()>0 && line.charAt(0)!='#') {
					GFF3GenomicFeature feature;
					try {
						feature = loadFeature(line); 
					} catch (RuntimeException e) {
						log.warning("Can not load genomic feature at line: "+line+". Unrecognized sequence name. "+e.getMessage());
						line=in.readLine();
						continue;
					}
					int idx = Arrays.binarySearch(supportedFeatureTypes, feature.getType());
					if(idx>=0) {
						//Memory saver for feature types
						feature.setType(supportedFeatureTypes[idx]);
						String id = feature.getId();
						if(id!=null) features.put(id,feature);
						if(feature.getParentId()!=null) {
							featuresWithParent.add(feature);
						}
					}
				}
				line=in.readLine();
			}
		}
		
		//Assign parent for features with parent id
		for(GFF3GenomicFeature feature:featuresWithParent) {
			String parentId = feature.getParentId();
			GFF3GenomicFeature parent = features.get(parentId);
			if(parent!=null) parent.addChild(feature);
			else log.warning("Parent not found for feature with id: "+feature.getId()+" location "+feature.getSequenceName()+":"+feature.getFirst()+"-"+feature.getLast());
		}
		int numGenes = 0;
		int numTranscripts = 0;
		for(GFF3GenomicFeature feature:features.values()) {
			if(FEATURE_TYPE_GENE.equals(feature.getType()) || FEATURE_TYPE_TRGENE.equals(feature.getType())) {
				String name = feature.getAnnotation("Name");
				if(name ==null) name = feature.getAnnotation("name");
				Gene gene = new Gene(feature.getId(), name );
				numGenes++;
				for(GFF3GenomicFeature geneChild:feature.getChildren()) {
					if(FEATURE_TYPE_MRNA.equals(geneChild.getType())) {
						Transcript transcript = new Transcript(geneChild.getId(), geneChild.getSequenceName(), geneChild.getFirst(), geneChild.getLast(), geneChild.isNegativeStrand() );
						transcript.setGene(gene);
						List<Exon> exons = new ArrayList<Exon>();
						for(GFF3GenomicFeature transcriptChild:geneChild.getChildren()) {
							String featureType = transcriptChild.getType();
							if (featureType.startsWith(FEATURE_TYPE_5PUTR)) {
								Exon exon = new Exon(transcript, transcriptChild.getFirst(), transcriptChild.getLast());
								exon.setStatus(Exon.STATUS_5P_UTR);
								exons.add(exon);
							} else if (FEATURE_TYPE_3PUTR.equals(featureType)) {
								Exon exon = new Exon(transcript, transcriptChild.getFirst(), transcriptChild.getLast());
								exon.setStatus(Exon.STATUS_3P_UTR);
								exons.add(exon);
							} else if (FEATURE_TYPE_CDS.equals(featureType)) {
								Exon exon = new Exon(transcript, transcriptChild.getFirst(), transcriptChild.getLast());
								exon.setStatus(Exon.STATUS_CODING);
								exons.add(exon);
							}
						}
						transcript.setExons(exons);
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
	public GFF3GenomicFeature loadFeature (String line) {
		String [] items = ParseUtils.parseString(line, '\t');
		//Memory saver for sequenceNames
		QualifiedSequence seq = sequenceNames.addOrLookupName(items[0]);
		GFF3GenomicFeature answer = new GFF3GenomicFeature(seq.getName(), Integer.parseInt(items[3]), Integer.parseInt(items[4]));
		answer.setType(items[2]);
		answer.setSource(items[1]);
		answer.setNegativeStrand(items[6].charAt(0)=='-');
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
class GFF3GenomicFeature implements GenomicRegion {
	private String sequenceName;
	private String source;
	private String type;
	private int first;
	private int last;
	private boolean negativeStrand;
	private Map<String,String> annotations = new TreeMap<String, String>();
	private List<GFF3GenomicFeature> children = new ArrayList<GFF3GenomicFeature>();
	
	public GFF3GenomicFeature(String sequenceName, int first, int last) {
		super();
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
	}
	public String getSequenceName() {
		return sequenceName;
	}
	
	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}
	
	public String getSource() {
		return source;
	}
	public String getType() {
		return type;
	}
	
	public void setType(String type) {
		this.type = type;
	}
	public int getFirst() {
		return first;
	}
	public int getLast() {
		return last;
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
	public void addChild(GFF3GenomicFeature child) {
		children.add(child);
	}
	public List<GFF3GenomicFeature> getChildren() {
		return children;
	}
	public String getId() {
		return annotations.get("ID");
	}
	public String getParentId() {
		return annotations.get("Parent");
	}
	public void setSource(String source) {
		this.source = source;
	}
	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}
}