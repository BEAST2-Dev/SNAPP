package snap.util;

import java.util.HashMap;
import java.util.List;
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

	static String getTopology(Node node, List<String> sLabels) {
		if (node.isLeaf()) {
			return sLabels.get(node.getNr()) + "[theta" +node.getNr() + "]:height"+node.getNr(); 
		} else {
			return "(" + getTopology(node.m_left, sLabels) + "," + getTopology(node.m_right, sLabels)  
			+")[theta" +node.getNr() + "]:height"+node.getNr(); 
		}
	}
	
	static String getHeader(Node node) {
		if (node.isLeaf()) {
			return "theta" +node.getNr() + " height"+node.getNr(); 
		} else {
			return getHeader(node.m_left) + " " + getHeader(node.m_right) + " theta" +node.getNr() + " height"+node.getNr();
		}
	}
	static double getHeight(Node node) {
		if (node.isLeaf()) {
			return node.m_fLength;
		} else {
			return getHeight(node.m_left) + node.m_fLength; 
		}
	}
	static String getTheta(String sTheta) {
		return sTheta.replaceAll("theta=", "");
	}
	static String getTreeData(Node node) {
		if (node.isLeaf()) {
			return getTheta(node.m_sMetaData) + " "+(getHeight(node) - node.m_fLength); 
		} else {
			return getTreeData(node.m_left) + " " + getTreeData(node.m_right) + " " + getTheta(node.m_sMetaData) + " "+(getHeight(node) - node.m_fLength);
		}
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
			
			System.out.println("#nr coverage tree");
			int j = 0;
			for (i = 0; i < m_nTopologies; i++) {
				System.out.println("#Tree " + i + ". " + m_fTreeWeight[i]*100 + "% " + getTopology(m_cTrees[i], sLabels));//m_cTrees[i].toString(sLabels));
				System.out.println("nr " + getHeader(m_cTrees[i]));
				boolean bSameTree = true;
				while (bSameTree) {
					System.out.println(getTreeData(m_trees[j]));
					j++;
					bSameTree = (j < m_trees.length) && (m_nTopology[j] == m_nTopology[j-1]);
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	} // main
}

/**
 * 
#Tree 1.  Frequency = 500. 
# ((A[theta = theta0, height = 0],B[theta=theta1, height=0]]):[theta=theta2,height=h0],C:[theta=theta3,height=0]):[theta=theta4,height=height1]
nr        samplenr theta0    ...          theta4        height0        …        height1
1        120        0.2        ….        0.3           1.1                     2.0
2        180        0.1        ….        0.1           1.2                     2.1
3        210        0.3        ….        0.2           1.3                     2.0
...
500      121        0.3        ….        0.2           1.3                     2.0 
**/