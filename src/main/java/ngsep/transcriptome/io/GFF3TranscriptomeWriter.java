package ngsep.transcriptome.io;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import ngsep.genome.GenomicRegionImpl;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.TranscriptSegment;
import ngsep.transcriptome.Transcriptome;

public class GFF3TranscriptomeWriter {
	public void printTranscriptome (Transcriptome transcriptome, PrintStream out) {
		List<Transcript> transcripts = transcriptome.getAllTranscripts();
		Set<String> geneIds = new HashSet<>();
		out.println("##gff-version 3");
		for(Transcript t:transcripts) {
			Gene gene = t.getGene();
			String geneId = gene.getId();
			if(!geneIds.contains(geneId)) {
				printGene (gene, out);
				geneIds.add(geneId);
			}
			printTranscript(t, out);
		}
	}

	private void printGene(Gene gene, PrintStream out) {
		String negStr = gene.isNegativeStrand()?"-":"+";
		out.println(gene.getSequenceName()+"\tNGSEP\t"+GFF3TranscriptomeHandler.FEATURE_TYPE_GENE+"\t"+gene.getFirst()+"\t"+gene.getLast()+"\t.\t"+negStr+"\t.\tID="+gene.getId());
	}

	private void printTranscript(Transcript transcript, PrintStream out) {
		String negStr = transcript.isNegativeStrand()?"-":"+";
		//Write mRNA feature
		out.println(transcript.getSequenceName()+"\tNGSEP\t"+GFF3TranscriptomeHandler.FEATURE_TYPE_MRNA+"\t"+transcript.getFirst()+"\t"+transcript.getLast()+"\t.\t"+negStr+"\t.\tID="+transcript.getId()+";Parent="+transcript.getGeneId());
		List<TranscriptSegment> segments = transcript.getTranscriptSegments();
		printRawExons(transcript, segments, negStr, out);
		//int sId = 1;
		for(TranscriptSegment segment:segments) {
			String type = GFF3TranscriptomeHandler.FEATURE_TYPE_TRGENE;
			String phase = ".";
			if(segment.getStatus()==TranscriptSegment.STATUS_5P_UTR) type = GFF3TranscriptomeHandler.FEATURE_TYPE_5PUTR;
			if(segment.getStatus()==TranscriptSegment.STATUS_3P_UTR) type = GFF3TranscriptomeHandler.FEATURE_TYPE_3PUTR;
			if(segment.getStatus()==TranscriptSegment.STATUS_CODING) {
				type = GFF3TranscriptomeHandler.FEATURE_TYPE_CDS;
				phase = ""+segment.getFirstCodonPositionOffset();
			}
			//out.println(transcript.getSequenceName()+"\tNGSEP\t"+type+"\t"+segment.getFirst()+"\t"+segment.getLast()+"\t.\t"+negStr+"\t"+phase+"\tID="+transcript.getId()+"_s"+sId+";Parent="+transcript.getId());
			out.println(transcript.getSequenceName()+"\tNGSEP\t"+type+"\t"+segment.getFirst()+"\t"+segment.getLast()+"\t.\t"+negStr+"\t"+phase+"\tParent="+transcript.getId());
			//sId++;
		}
	}

	private void printRawExons(Transcript t, List<TranscriptSegment> segments, String negStr, PrintStream out) {
		GenomicRegionImpl nextExon = null;
		int exId = 1;
		for(TranscriptSegment segment:segments) {
			if(nextExon == null || nextExon.getLast()+1<segment.getFirst()) {
				if(nextExon!=null) {
					out.println(nextExon.getSequenceName()+"\tNGSEP\t"+GFF3TranscriptomeHandler.FEATURE_TYPE_EXON+"\t"+nextExon.getFirst()+"\t"+nextExon.getLast()+"\t.\t"+negStr+"\t.\tID="+t.getId()+"_e"+exId+";Parent="+t.getId());
					exId++;
				}
				nextExon = new GenomicRegionImpl(segment.getSequenceName(), segment.getFirst(), segment.getLast());
			}
			if(segment.getLast()>nextExon.getLast()) nextExon.setLast(segment.getLast());
		}
		if(nextExon!=null) {
			out.println(nextExon.getSequenceName()+"\tNGSEP\t"+GFF3TranscriptomeHandler.FEATURE_TYPE_EXON+"\t"+nextExon.getFirst()+"\t"+nextExon.getLast()+"\t.\t"+negStr+"\t.\tID="+t.getId()+"_e"+exId+";Parent="+t.getId());
		}
	}
	
}
