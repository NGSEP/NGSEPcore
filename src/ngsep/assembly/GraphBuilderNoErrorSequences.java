package ngsep.assembly;

import java.util.List;

import ngsep.sequences.DNAMaskedSequence;
import static ngsep.assembly.TimeUtilities.timeIt;
import static ngsep.assembly.TimeUtilities.progress;

public class GraphBuilderNoErrorSequences implements GraphBuilder {
	private SimplifiedAssemblyGraph sag;
	private List<CharSequence> sequences;

	@Override
	public SimplifiedAssemblyGraph buildAssemblyGraph(List<CharSequence> sequences) {
		sag = new SimplifiedAssemblyGraph(sequences);
		this.sequences = sequences;
		timeIt("--- buildTheGraph", () -> {
			int total = (sequences.size() * (sequences.size() - 1)) / 2;
			int count = 0;
			for (int i = 0; i < sequences.size(); i++) {
				for (int j = i + 1; j < sequences.size(); j++) {
					progress("--- buildTheGraph", count, total);
					aling(i, false, j, false); // j->i || i'->j'
					aling(i, false, j, true); // j'->i || i'->j
					aling(i, true, j, false); // j->i' || i ->j'
					aling(i, true, j, true); // j'->i' || i ->j
					count++;
				}
			}
			sag.removeAllEmbeddedsIntoGraph();
		});
		return sag;
	}

	private void aling(int id1, boolean r1, int id2, boolean r2) {
		CharSequence a = r1 ? DNAMaskedSequence.getReverseComplement(sequences.get(id1).toString())
				: sequences.get(id1);
		CharSequence b = r2 ? DNAMaskedSequence.getReverseComplement(sequences.get(id2).toString())
				: sequences.get(id2);

		int i = 0;
		int j = b.length() - 1;
		while (j >= 0) {
			if (a.charAt(i) == b.charAt(j)) {
				i++;
				j++;
				if (i >= a.length() || j >= b.length()) {
					decide(id1, r1, id2, r2, j, i + 1);
					break;
				}
			} else {
				j -= i + 1;
				i = 0;
			}
		}
	}

	private void decide(int id1, boolean r1, int id2, boolean r2, int pos, int len) {
		if (pos + len >= sequences.get(id2).length()) {
			// id2^r2->id1^r1
			sag.addEdge((id2 << 1) + (r2 ? 0 : 1), (id1 << 1) + (r1 ? 1 : 0), len, 1);
		} else {
			// id1^r1 into id2^r2
			boolean reversed = r1 ^ r2;
			pos = (r2 && !r1) ? sequences.get(id2).length() - pos - len : pos;
			sag.addEmbedded(id2, id1, pos, reversed, 1);
		}
	}

	public static void main(String[] args) throws Exception {
		System.out.println("¡¡¡ Esta funcion es solo para desarrollo !!!!");
		List<CharSequence> a = Assembler.load(args[0]);
		SimplifiedAssemblyGraph assemblyGraph = (new GraphBuilderNoErrorSequences()).buildAssemblyGraph(a);
		assemblyGraph.save(args[1]);
	}

}
