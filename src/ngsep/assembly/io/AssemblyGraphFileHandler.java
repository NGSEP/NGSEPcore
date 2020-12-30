package ngsep.assembly.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import ngsep.assembly.AssemblyEdge;
import ngsep.assembly.AssemblyEmbedded;
import ngsep.assembly.AssemblyGraph;
import ngsep.assembly.AssemblyVertex;
import ngsep.main.io.ConcatGZIPInputStream;
import ngsep.sequences.QualifiedSequence;

public class AssemblyGraphFileHandler {
	public static void save(AssemblyGraph graph, String outFileGraph) throws IOException {
		try (OutputStream os = new GZIPOutputStream(new FileOutputStream(outFileGraph));
			 PrintStream out = new PrintStream(os)) {
			List<QualifiedSequence> sequences = graph.getSequences();
			out.println("#SEQUENCES");
			int id=0;
			for(QualifiedSequence seq:sequences) {
				out.println(""+id+"\t"+seq.getName()+"\t"+seq.getLength());
				id++;
			}
			
			out.println("#EMBEDDED");
			for(int i=0;i<sequences.size();i++) {
				List<AssemblyEmbedded> embeddedList = graph.getEmbeddedBySequenceId(i);
				for(AssemblyEmbedded embedded:embeddedList) {
					saveEmbedded(embedded, out);
				}
			}
			
			out.println("#EDGES");
			List<AssemblyEdge> edges = graph.getEdges();
			for(AssemblyEdge edge:edges) {
				saveEdge(edge, out);
			}
			
		}
	}
	
	public static AssemblyGraph load(List<QualifiedSequence> sequences, String graphFilename) throws IOException {
		AssemblyGraph graph = new AssemblyGraph(sequences);
		String line = null;
		try (ConcatGZIPInputStream gzs = new ConcatGZIPInputStream(new FileInputStream(graphFilename));
			 BufferedReader in = new BufferedReader(new InputStreamReader(gzs))) {
			line = in.readLine();
			if(!"#SEQUENCES".equals(line)) throw new IOException("Graph file misses sequence names. First line: "+line);
			int seqId = 0;
			line=in.readLine();
			while(line!=null && seqId<sequences.size() && !line.startsWith("#")) {
				String [] items = line.split("\t");
				QualifiedSequence seq = sequences.get(seqId);
				if(!seq.getName().equals(items[1]))  throw new IOException("Unexpected sequence name at index " +(seqId+2)+". Double check that the graph was built from the given sequences or build again the graph. Expected: "+seq.getName()+" "+seq.getLength()+" loaded: "+line);
				if(seq.getLength()!=Integer.parseInt(items[2])) throw new IOException("Unexpected sequence length at index" +(seqId+2)+". Sequence name: "+seq.getName()+". Double check that the graph was built from the given sequences or build again the graph. Expected: "+seq.getLength()+" loaded: "+items[1]);
				seqId++;
				line=in.readLine();
			}
			if(line==null) throw new IOException("Unexpected end of file reading sequences");
			//if(seqId<sequences.size()) log.info("Missing sequences in graph. Expected: "+sequences.size()+". Loaded: "+seqId);
			if (!line.equals("#EMBEDDED")) throw new IOException("Unexpected line after loading sequences. Expected: #EMBEDDED. Line: "+line);
			line=in.readLine();
			while(line!=null && !line.startsWith("#")) {
				String [] items = line.split("\t");
				AssemblyEmbedded embedded = loadEmbedded(items, sequences);
				graph.addEmbedded(embedded);
				line=in.readLine();
			}
			if(line==null) throw new IOException("Unexpected end of file reading embedded relationships");
			if (!line.equals("#EDGES")) throw new IOException("Unexpected line after loading embedded. Expected: #EDGES. Line: "+line);
			line=in.readLine();
			while(line!=null && !line.startsWith("#")) {
				String [] items = line.split("\t");
				AssemblyEdge edge = loadEdge(items, graph);
				if(!edge.isSameSequenceEdge()) graph.addEdge(edge);
				line=in.readLine();
			}
		} catch (NumberFormatException e) {
			throw new IOException("Error loading number at line: "+line,e);
		}
		graph.updateVertexDegrees();
		return graph;
	}

	
	
	private static void saveEmbedded(AssemblyEmbedded embedded, PrintStream out) {
		int reverse = embedded.isReverse()?1:0;
		out.print(""+embedded.getSequenceId()+"\t"+embedded.getHostId()+"\t"+embedded.getHostStart()+"\t"+embedded.getHostEnd()+"\t"+reverse);
		out.print("\tHSD="+embedded.getHostStartStandardDeviation()+";NSK="+embedded.getNumSharedKmers()+";CSK="+embedded.getCoverageSharedKmers());
		out.print(";WCSK="+embedded.getWeightedCoverageSharedKmers()+";RNSK="+embedded.getRawKmerHits()+";RHSD="+embedded.getRawKmerHitsSubjectStartSD());
		out.println(";HES="+embedded.getHostEvidenceStart()+";HEE="+embedded.getHostEvidenceEnd());
	}
	
	private static AssemblyEmbedded loadEmbedded(String[] items, List<QualifiedSequence> sequences) {
		int embSeqId = Integer.parseInt(items[0]);
		int hostId = Integer.parseInt(items[1]);
		boolean reverse = Integer.parseInt(items[4])==1;
		QualifiedSequence embeddedSeq = sequences.get(embSeqId);
		AssemblyEmbedded embedded = new AssemblyEmbedded(embSeqId, embeddedSeq, reverse, hostId, Integer.parseInt(items[2]), Integer.parseInt(items[3]));
		String [] items2 = items[5].split(";");
		for(String feature:items2) {
			int idxF = feature.indexOf("=");
			String featureName = feature.substring(0,idxF);
			int featureValue = Integer.parseInt(feature.substring(idxF+1));
			if("HSD".equals(featureName)) embedded.setHostStartStandardDeviation(featureValue);
			else if("NSK".equals(featureName)) embedded.setNumSharedKmers(featureValue);	
			else if("CSK".equals(featureName)) embedded.setCoverageSharedKmers(featureValue);
			else if("WCSK".equals(featureName)) embedded.setWeightedCoverageSharedKmers(featureValue);
			else if("RNSK".equals(featureName)) embedded.setRawKmerHits(featureValue);
			else if("RHSD".equals(featureName)) embedded.setRawKmerHitsSubjectStartSD(featureValue);
			else if("HES".equals(featureName)) embedded.setHostEvidenceStart(featureValue);
			else if("HEE".equals(featureName)) embedded.setHostEvidenceEnd(featureValue);

		}
		return embedded;
	}

	private static void saveEdge(AssemblyEdge edge, PrintStream out) {
		out.print(""+edge.getVertex1().getUniqueNumber()+"\t"+edge.getVertex2().getUniqueNumber()+"\t"+edge.getOverlap());
		out.print("\tOSD="+edge.getOverlapStandardDeviation()+";NSK="+edge.getNumSharedKmers()+";CSK="+edge.getCoverageSharedKmers());
		out.print(";WCSK="+edge.getWeightedCoverageSharedKmers()+";RNSK="+edge.getRawKmerHits()+";RSSD="+edge.getRawKmerHitsSubjectStartSD());
		out.println(";AO="+edge.getAverageOverlap()+";MO="+edge.getMedianOverlap()+";FLO="+edge.getFromLimitsOverlap());
	}
	
	private static AssemblyEdge loadEdge(String[] items, AssemblyGraph graph) {
		int v1Idx = Integer.parseInt(items[0]);
		int v2Idx = Integer.parseInt(items[1]);
		int overlap = Integer.parseInt(items[2]);
		AssemblyVertex v1 = graph.getVertexByUniqueId(v1Idx);
		AssemblyVertex v2 = graph.getVertexByUniqueId(v2Idx);
		AssemblyEdge edge;
		//if(v1Idx==49 || v2Idx == 49) System.out.println("Loading edge between: "+v1+" and "+v2+". SeqIdx1: "+v1.getSequenceIndex()+" SeqIdxV2: "+v2.getSequenceIndex());
		if(v1.getSequenceIndex()==v2.getSequenceIndex()) {
			edge = graph.getSameSequenceEdge(v1);
			//if(v1Idx==49 || v2Idx == 49) System.out.println("retrieved edge: "+edge+" total edges vertex 1: "+graph.getEdges(v1).size()+" total edges vertex 2: "+graph.getEdges(v2).size());
		} else {
			edge = new AssemblyEdge(v1, v2, overlap);
			//if(v1Idx==49 || v2Idx == 49) System.out.println("Created edge: "+edge);
		}
		String [] items2 = items[3].split(";");
		for(String feature:items2) {
			int idxF = feature.indexOf("=");
			String featureName = feature.substring(0,idxF);
			int featureValue = Integer.parseInt(feature.substring(idxF+1));
			if("OSD".equals(featureName)) edge.setOverlapStandardDeviation(featureValue);
			else if("NSK".equals(featureName)) edge.setNumSharedKmers(featureValue);	
			else if("CSK".equals(featureName)) edge.setCoverageSharedKmers(featureValue);
			else if("WCSK".equals(featureName)) edge.setWeightedCoverageSharedKmers(featureValue);
			else if("RNSK".equals(featureName)) edge.setRawKmerHits(featureValue);
			else if("RSSD".equals(featureName)) edge.setRawKmerHitsSubjectStartSD(featureValue);
			else if("AO".equals(featureName)) edge.setAverageOverlap(featureValue);
			else if("MO".equals(featureName)) edge.setMedianOverlap(featureValue);
			else if("FLO".equals(featureName)) edge.setFromLimitsOverlap(featureValue);
		}
		return edge;
	}
	
	public static List<QualifiedSequence> loadSequenceNamesFromGraphFile(String graphFilename) throws IOException {
		List<QualifiedSequence> sequenceNames = new ArrayList<QualifiedSequence>();
		String line = null;
		try (ConcatGZIPInputStream gzs = new ConcatGZIPInputStream(new FileInputStream(graphFilename));
			 BufferedReader in = new BufferedReader(new InputStreamReader(gzs))) {
			line = in.readLine();
			if(!"#SEQUENCES".equals(line)) throw new IOException("Graph file misses sequence names. First line: "+line);
			line=in.readLine();
			while(line!=null && !line.startsWith("#")) {
				String [] items = line.split("\t");
				QualifiedSequence seq = new QualifiedSequence(items[1]);
				seq.setLength(Integer.parseInt(items[2]));
				sequenceNames.add(seq);
				line=in.readLine();
			}
		} catch (NumberFormatException e) {
			throw new IOException("Error loading number at line: "+line,e);
		}
		return sequenceNames;
	}
}
