package ngsep.transcriptome.io;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import ngsep.genome.ReferenceGenome;
import ngsep.transcriptome.Codon;
import ngsep.transcriptome.ProteinTranslator;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.Transcriptome;

public class GFF3CombineAnnotations {
	private ReferenceGenome genome;
	public static void main(String[] args) throws Exception {
		GFF3CombineAnnotations instance = new GFF3CombineAnnotations();
		instance.genome = new ReferenceGenome(args[0]);
		String gff1 = args[1];
		String gff2 = args[2];
		instance.process(gff1,gff2,System.out);

	}
	private void process(String gff1, String gff2, PrintStream out) throws IOException {
		GFF3TranscriptomeHandler handler = new GFF3TranscriptomeHandler(genome.getSequencesMetadata());
		Transcriptome t1 = handler.loadMap(gff1);
		t1.fillSequenceTranscripts(genome,null);
		System.err.println("Loaded transcriptome 1 with "+t1.getAllTranscripts().size()+" transcripts");
		Transcriptome t2 = handler.loadMap(gff2);
		t2.fillSequenceTranscripts(genome,null);
		List<Transcript> transcripts = t2.getAllTranscripts();
		System.err.println("Loaded transcriptome 2 with "+transcripts.size()+" transcripts");
		List<Transcript> toAdd = new ArrayList<>();
		List<Transcript> toRemove = new ArrayList<>();
		ProteinTranslator translator = new ProteinTranslator();
		for(Transcript t:transcripts) {
			String protein = t.getProteinSequence(translator);
			if(protein == null) continue;
			List<Transcript> overlappingTranscripts = t1.getTranscripts(t).asList();
			if(overlappingTranscripts.isEmpty()) {
				System.err.println("Transcript: "+t.getId()+" does not overlap with transcripts. Adding transcript");
				toAdd.add(t);
				continue;
			}
			if(!overlapWithCodingSegments(t, overlappingTranscripts)) {
				System.err.println("Transcript: "+t.getId()+" overlaps with "+overlappingTranscripts.size()+" transcripts but does not overlap with their coding segments. Adding transcript");
				toAdd.add(t);
				continue;
			}
			Transcript betterTranscript = null;
			for(Transcript existing:overlappingTranscripts) {
				Codon start = existing.getStartCodon();
				Codon stop = existing.getStopCodon();
				String existingProtein = existing.getProteinSequence(translator);
				//if(t.getGeneId().contains("gene-0.48")) System.err.println("Existing: "+existing.getId()+" coding: "+existing.isCoding()+" start: "+start+" stop: "+stop);
				//if(existing.getId().equals("MSTRG.17.1")) System.err.println("Existing: "+existing.getId()+" coding: "+existing.isCoding()+" start: "+start+" stop: "+stop+" protein length: "+existingProtein.length()+" compared to: "+t.getId()+" length: "+protein.length());
				if(existing.isCoding() && start!=null && start.isStart() && stop!=null && stop.isStop() && existingProtein.length()>= protein.length()) {
					betterTranscript = existing;
					System.err.println("Transcript: "+t.getId()+" filtered by better transcript: "+existing.getId());
					break;
				}
			}
			if(betterTranscript==null) {
				toRemove.addAll(overlappingTranscripts);
				toAdd.add(t);
				System.err.println("Transcript: "+t.getId()+" removed gene: "+overlappingTranscripts.get(0).getGeneId());
			}
		}
		System.err.println("Total transcripts to remove "+toRemove.size());
		for(Transcript t:toRemove) {
			t1.removeTranscript(t);
		}
		System.err.println("Total transcripts to add "+toAdd.size());
		for(Transcript t:toAdd) {
			t1.addTranscript(t);
		}
		GFF3TranscriptomeWriter writer = new GFF3TranscriptomeWriter();
		writer.printTranscriptome(t1, out);
	}
	private boolean overlapWithCodingSegments(Transcript t, List<Transcript> overlappingTranscripts) {
		for(Transcript overlapping:overlappingTranscripts) {
			if(overlapping.getTranscriptSegmentsByAbsolute(t.getFirst(), t.getLast()).size()>0) return true;
		}
		return false;
	}
}
