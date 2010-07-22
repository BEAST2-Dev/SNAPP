package snap.util;

import java.util.HashMap;
import java.util.Vector;

import snap.util.TreeFileParser;

public class TreeSetAnalyser {

	/**
	 * move divide y-position of a tree with factor f. Useful to calculate
	 * consensus trees.
	 **/
	static void divideLength(Node node, float f) {
		if (!node.isLeaf()) {
			divideLength(node.m_left, f);
			divideLength(node.m_right, f);
		}
		node.m_fLength /= f;
	}

	/**
	 * add length of branches in src to that of target Useful to calculate
	 * consensus trees. Assumes src and target share same topology
	 */
	static void addLength(Node src, Node target) {
		// assumes same topologies for src and target
		if (!src.isLeaf()) {
			addLength(src.m_left, target.m_left);
			addLength(src.m_right, target.m_right);
		}
		target.m_fLength += src.m_fLength;
	}

	public static void main(String [] args) {
		try {
			Vector<String> sLabels = new Vector<String>();
			String sFile = args[args.length - 1];
			TreeFileParser parser = new TreeFileParser(sLabels, null, null, 0);
			Node [] m_trees = parser.parseFile(sFile);
		
			// count tree topologies
			// first step is find how many different topologies are present
			int [] m_nTopology = new int[m_trees.length];
			HashMap<String, Integer> map = new HashMap<String, Integer>();
			for (int i = 0; i < m_trees.length; i++) {
				Node tree = m_trees[i];
				String sNewick = tree.toShortNewick();
				if (map.containsKey(sNewick)) {
					m_nTopology[i] = map.get(sNewick).intValue();
				} else {
					m_nTopology[i] = map.size();
					map.put(sNewick, map.size());
				}
			}

			// second step is find how many different tree have a particular
			// topology
			int m_nTopologies = map.size();
			int[] nTopologies = new int[map.size()];
			for (int i = 0; i < m_trees.length; i++) {
				nTopologies[m_nTopology[i]]++;
			}

			// sort the trees so that frequently occurring topologies go first in
			// the ordering
			for (int i = 0; i < m_trees.length; i++) {
				for (int j = i + 1; j < m_trees.length; j++) {
					if (nTopologies[m_nTopology[i]] < nTopologies[m_nTopology[j]]
							|| (nTopologies[m_nTopology[i]] == nTopologies[m_nTopology[j]] && m_nTopology[i] > m_nTopology[j])) {
						int h = m_nTopology[j];
						m_nTopology[j] = m_nTopology[i];
						m_nTopology[i] = h;
						Node tree = m_trees[j];
						m_trees[j]=m_trees[i];
						m_trees[i] = tree;
					}

				}
			}

			int i = 0;
			int iOld = 0;
			int iConsTree = 0;
			float [] m_fTreeWeight = new float[m_nTopologies];
			Node [] m_cTrees = new Node[m_nTopologies];
			while (i < m_trees.length) {
				Node tree = m_trees[i].copy();
				Node consensusTree = tree;
				i++;
				while (i < m_trees.length && m_nTopology[i] == m_nTopology[i - 1]) {
					tree = m_trees[i];
					addLength(tree, consensusTree);
					i++;
				}
				divideLength(consensusTree, i - iOld);
				m_fTreeWeight[iConsTree] = (float) (i - iOld + 0.0) / m_trees.length;
				// position nodes of consensus trees
////				positionLeafs(consensusTree);
////				positionRest(consensusTree);
//				float fHeight = positionHeight(consensusTree, 0);
//				offsetHeight(consensusTree, m_fHeight - fHeight);
				m_cTrees[iConsTree] = consensusTree;
				iConsTree++;
				iOld = i;
			}
			int [] m_nTopologyByPopularity = new int[m_trees.length];
			int nColor = 0;
			m_nTopologyByPopularity[0] = 0; 
			for (i = 1; i < m_trees.length; i++) {
				if (m_nTopology[i] != m_nTopology[i-1]) {
					nColor++;
				}
				m_nTopologyByPopularity[i] = nColor; 
			}
			
			System.out.println("nr coverage tree");
			for (i = 0; i < m_nTopologies; i++) {
				System.out.println(i + " " + m_fTreeWeight[i]*100 + "% " + m_cTrees[i].toString(sLabels));
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	} // main
}
