package ngsep.assembly;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;
import java.util.stream.Collectors;
import ngsep.assembly.GraphQualityAnalyzer.Sequence;
import ngsep.sequences.DNAMaskedSequence;

public class PathQualityAnalizer {

	private SimplifiedAssemblyGraph lec;
	private SimplifiedAssemblyGraph ref;
	private List<Sequence> nams;
	private Map<String, List<Sequence>> mapeds;

	public PathQualityAnalizer(List<Sequence> pNams) throws FileNotFoundException {
		this.nams = pNams;
		Collections.sort(nams, (l1, l2) -> l2.sequence.length() - l1.sequence.length());
		int i = 0;
		for (Sequence seq : nams)
			seq.id = i++;

		this.ref = getGraph(nams);

		AssemblyConfiguration ac = new AssemblyConfiguration();
		GraphBuilderFMIndex builder = new GraphBuilderFMIndex();
		builder.setConfig(ac);
		this.lec = builder
				.buildSimplifiedAssemblyGraph(this.nams.stream().map(a -> a.sequence).collect(Collectors.toList()));

		// can't be sure of the map thanks of the embedes
		for (Entry<Integer, Map<Integer, Embedded>> entry : ref.getEmbbeded().entrySet()) {
			int id = entry.getKey();
			for (Entry<Integer, Embedded> entry2 : entry.getValue().entrySet())
				lec.addEmbedded(id, entry2.getKey(), entry2.getValue().getPos(), entry2.getValue().isReversed(),
						entry2.getValue().getRate());
		}
		for (Entry<Integer, Map<Integer, Embedded>> entry : lec.getEmbbeded().entrySet()) {
			int id = entry.getKey();
			for (Entry<Integer, Embedded> entry2 : entry.getValue().entrySet())
				ref.addEmbedded(id, entry2.getKey(), entry2.getValue().getPos(), entry2.getValue().isReversed(),
						entry2.getValue().getRate());
		}

		lec.removeAllEmbeddedsIntoGraph();
		lec.removeDuplicatedEmbeddes();
		ref.removeAllEmbeddedsIntoGraph();
		ref.removeDuplicatedEmbeddes();

		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		int a = 0;
		for (int j = 0; j < nams.size(); j++) {
			if (!lec.isEmbedded(j))
				map.put(j, a++);
		}
		List<List<Integer>> paths = new ArrayList<List<Integer>>(mapeds.size());
		for (List<Sequence> seqs : mapeds.values()) {
			List<Integer> path = new LinkedList<Integer>();
			for (Sequence seq : seqs) {
				if (!lec.isEmbedded(seq.id)) {
					int id = 2*map.get(seq.id);
					path.add(id + ((seq.rev) ? 1 : 0));
					path.add(id + ((seq.rev) ? 0 : 1));
				}
			}
			paths.add(path);
		}
		System.out.println(paths);

		LayourBuilder builder2 = new LayoutBuilderGreedy();
		AssemblyGraph assemblyGraph = lec.getAssemblyGraph();
		builder2.findPaths(assemblyGraph);
		System.out.println("-------------Graph------------------");
		lec.printInfo();
		emmbededTest();
		
		List<List<Integer>> paths2 = new ArrayList<List<Integer>>(mapeds.size());
		for (List<AssemblyEdge> edge : assemblyGraph.getPaths()) {
			List<Integer> path = new LinkedList<Integer>();
			int m = edge.get(0).getVertex1().getIndex();
			int m1 = edge.get(0).getVertex2().getIndex();
			int n = edge.get(1).getVertex1().getIndex();
			int n1 = edge.get(1).getVertex2().getIndex();
			int curr;
			if (m == n || m == n1) {
				path.add(m1);
				curr = m;
			} else {
				path.add(m);
				curr = m1;
			}
			for (int k = 1; k < edge.size(); k++) {
				path.add(curr);
				m = edge.get(k).getVertex1().getIndex();
				m1 = edge.get(k).getVertex2().getIndex();
				curr = (curr == m) ? m1 : m;
			}
			path.add(curr);
			paths2.add(path);
		}
		System.out.println(paths2);

	}

	
	public void emmbededTest() {
		Set<Integer> refEmb = new HashSet<>();
		for (Map<Integer, Embedded> a : ref.getEmbbeded().values())
			refEmb.addAll(a.keySet());

		Set<Integer> lectEmb = new HashSet<>();
		for (Entry<Integer, Map<Integer, Embedded>> a : lec.getEmbbeded().entrySet())
			lectEmb.addAll(a.getValue().keySet());

		int trueP = 0, falseP = 0;
		for (Entry<Integer, Map<Integer, Embedded>> a : lec.getEmbbeded().entrySet()) {
			for (Entry<Integer, Embedded> b : a.getValue().entrySet()) {
				int i = b.getKey();
				if (refEmb.contains(i))
					trueP++;
				else {
					int h = nams.get(a.getKey()).pos + nams.get(a.getKey()).len;
					int j = nams.get(i).pos + nams.get(i).len;
					int tope = (int) (ref.getSequences().get(i).length() * 0.01);
					if (Math.abs(h - j) <= tope || Math.abs(nams.get(a.getKey()).pos - nams.get(i).pos) <= tope) {
						trueP++;
					} else {
						falseP++;
						
					}
					
				}

				StringBuilder strBuild = new StringBuilder();
				for (int k = 0; k < b.getValue().getPos(); k++)
					strBuild.append('-');
				String seq = (!b.getValue().isReversed())? DNAMaskedSequence.getReverseComplement(ref.getSequences().get(i)):ref.getSequences().get(i).toString();
			}

		}

		int falseN = 0;
		for (int i : refEmb) {
			if (!lectEmb.contains(i))
				falseN++;
		}
		int trueN = (ref.getSequences().size() - ref.amuontOfEmbeddedSequences()) - falseP;
		
		Set<String> refEdg = new HashSet<>();
		for (Entry<Integer, Map<Integer, Alignment>> map : ref.getEdges().entrySet()) {
			int id1 = map.getKey();
			if (!lectEmb.contains(id1 >> 1) && !refEmb.contains(id1 >> 1))
				for (int id2 : map.getValue().keySet())
					if (id1 + 1 < id2 && !lectEmb.contains(id2 >> 1) && !refEmb.contains(id2 >> 1))
						refEdg.add(id1 + "-" + id2);
		}

		Set<String> lecEdg = new HashSet<>();
		for (Entry<Integer, Map<Integer, Alignment>> map : lec.getEdges().entrySet()) {
			int id1 = map.getKey();
			if (!lectEmb.contains(id1 >> 1) && !refEmb.contains(id1 >> 1))
				for (int id2 : map.getValue().keySet())
					if (id1 + 1 < id2 && !lectEmb.contains(id2 >> 1) && !refEmb.contains(id2 >> 1))
						lecEdg.add(id1 + "-" + id2);
		}

		trueP = 0;
		falseP = 0;
		for (String i : lecEdg) {
			if (refEdg.contains(i))
				trueP++;
			else
				falseP++;
		}

		falseN = 0;
		for (String i : refEdg) {
			if (!lecEdg.contains(i))
				falseN++;
		}
		trueN = -1;

		System.out.println("edges (without embeddes in both graphs)");
		System.out.println("false|true");
		System.out.println("neg= " + falseN + "|" + trueN);
		System.out.println("pos= " + falseP + "|" + trueP);
		System.out.println("precision = " + (100 * trueP) / (double) (trueP + falseP) + " %");
		System.out.println("recall = " + (100 * trueP) / (double) (trueP + falseN) + " %");

		int s = 0, l = 0;
		for (Entry<Integer, Map<Integer, Alignment>> map : ref.getEdges().entrySet()) {
			int id1 = map.getKey();
			if (lectEmb.contains(id1) && refEmb.contains(id1))
				for (Entry<Integer, Alignment> a : map.getValue().entrySet()) {
					int id2 = a.getKey();
					if (id1 + 1 < id2 && lectEmb.contains(id2) && refEmb.contains(id2))
						if (!lecEdg.contains(id1 + "-" + id2)) {
							s += a.getValue().getOverlap();
							l++;
						}
				}
		}
		System.out.println("mean overlap (edgeFalsePositive) = " + (s / (double) l));

	}
	
	public SimplifiedAssemblyGraph getGraph(List<Sequence> sequences) throws FileNotFoundException {
		SimplifiedAssemblyGraph sag = new SimplifiedAssemblyGraph(GraphQualityAnalyzer.getSequences(sequences));

		Map<String, List<Sequence>> a = new HashMap<>();
		for (Sequence s : sequences)
			a.computeIfAbsent(s.ref, (x) -> new ArrayList<>()).add(s);

		for (List<Sequence> list : a.values())
			Collections.sort(list, (Sequence x, Sequence y) -> {
				int ans = x.pos - y.pos;
				if (ans != 0)
					return ans;
				return y.len - x.len;
			});

		mapeds = a;
		for (List<Sequence> list : a.values()) {
			for (int i = 0; i < list.size() - 1; i++) {
				Sequence lect1 = list.get(i);
				for (int j = i + 1; j < list.size() && list.get(j).pos < lect1.pos + lect1.len; j++) {
					Sequence lect2 = list.get(j);

					int relativePos = lect2.pos - lect1.pos;
					if (relativePos + lect2.len > lect1.len)
						// lect1 -> lect2
						sag.addEdge((lect1.id << 1) + (lect1.rev ? 0 : 1), (lect2.id << 1) + (lect2.rev ? 1 : 0),
								lect1.len - relativePos, 1);
					else {
						// lect2 into lect1
						boolean reversed = lect1.rev ^ lect2.rev;
						relativePos = (lect1.rev) ? lect1.len - lect2.len - relativePos : relativePos;
						sag.addEmbedded(lect1.id, lect2.id, relativePos, reversed, 1);
					}
				}
			}
		}

		sag.removeAllEmbeddedsIntoGraph();
		return sag;
	}

	public static void main(String[] args) throws FileNotFoundException, ClassNotFoundException, IOException {
		new PathQualityAnalizer(GraphQualityAnalyzer.load(args[0]));
	}

}
