package snap.util;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import snap.util.TreeFileParser;

public class TreeSetAnalyser {

	String m_sFileName;
	boolean m_bFirstTopologyOnly = false;
	Vector<String> m_sLabels = new Vector<String>();
	Node [] m_trees;
	int m_nTopologies;
	int [] m_nTopology;
	float [] m_fTreeWeight;
	Node [] m_cTrees;

	List<Double>[] m_heights;
	List<Double>[] m_thetas;

	double m_fBurnIn = 10;

	/**
	 * move divide y-position of a tree with factor f. Useful to calculate
	 * consensus trees.
	 **/
	void divideLength(Node node, float f) {
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
	void addLength(Node src, Node target) {
		// assumes same topologies for src and target
		if (!src.isLeaf()) {
			addLength(src.m_left, target.m_left);
			addLength(src.m_right, target.m_right);
		}
		target.m_fLength += src.m_fLength;
	}

	String getTopology(Node node, List<String> sLabels) {
		if (node.isLeaf()) {
			return sLabels.get(node.getNr()) + "[theta" +node.getNr() + "]:height"+node.getNr(); 
		} else {
			return "(" + getTopology(node.m_left, sLabels) + "," + getTopology(node.m_right, sLabels)  
			+")[theta" +node.getNr() + "]:height"+node.getNr(); 
		}
	}
	
	String getHeader(Node node) {
		if (node.isLeaf()) {
			return "theta" +node.getNr() + "\theight"+node.getNr(); 
		} else {
			return getHeader(node.m_left) + "\t" + getHeader(node.m_right) + "\ttheta" +node.getNr() + "\theight"+node.getNr();
		}
	}
	double getHeight(Node node) {
		if (node.isLeaf()) {
			return node.m_fLength;
		} else {
			return getHeight(node.m_left) + node.m_fLength; 
		}
	}

	double getTheta(String sMetaData) {
		if (sMetaData.contains("theta=")) {
			String sStr = sMetaData.replaceAll("theta=", "");
			double f = Double.parseDouble(sStr);
			return f;
		} else {
			String sStr = sMetaData.replaceAll("coalescenceRate=", "");
			double f = 2.0/Double.parseDouble(sStr);
			return f;
		}
	}
	
	
	String getTreeData(Node node) {
		if (node.isLeaf()) {
			double fTheta = getTheta(node.m_sMetaData);
			double fHeight = (getHeight(node) - node.m_fLength);
			m_thetas[node.getNr()].add(fTheta);
			m_heights[node.getNr()].add(fHeight);
			return fTheta + "\t" +fHeight; 
		} else {
			double fTheta = getTheta(node.m_sMetaData);
			double fHeight = (getHeight(node) - node.m_fLength);
			m_thetas[node.getNr()].add(fTheta);
			m_heights[node.getNr()].add(fHeight);
			return getTreeData(node.m_left) + "\t" + getTreeData(node.m_right) + "\t" + fTheta + "\t"+ fHeight;
		}
	}
	
	@SuppressWarnings("unchecked")
	void initLists(int nNodes) {
		m_heights = new List[nNodes];
		m_thetas = new List[nNodes];
		for (int i = 0; i < nNodes; i++) {
			m_heights[i] = new ArrayList<Double>();
			m_thetas[i] = new ArrayList<Double>();
		}
	}
	
	double mean(List<Double> fList) {
		double fSum = 0;
		for (double f : fList) {
			fSum += f;
		}
		return fSum / fList.size();
	}

	String truncate(double f, int nLength) {
		String sStr = f + "";
		if (sStr.length() > nLength) {
			sStr = sStr.substring(0, nLength);
		}
		return sStr;
	}
	void calcMean(Node node) {
		//List<Double> fHeights = m_heights[node.getNr()];
		List<Double> fThetas = m_thetas[node.getNr()];
		double fMeanTheta = mean(fThetas);
		node.m_sMetaData = "mTheta=" + truncate(fMeanTheta, 6);
		if (!node.isLeaf()) {
			calcMean(node.m_left);
			calcMean(node.m_right);
		}
	}
	String toString(Node node, List<String> sLabels) {
		StringBuffer buf = new StringBuffer();
		if (node.m_left != null) {
			buf.append("(");
			buf.append(toString(node.m_left, sLabels));
			buf.append(',');
			buf.append(toString(node.m_right, sLabels));
			buf.append(")");
		} else {
			buf.append(node.getNr());
		}
		if (node.m_sMetaData != null) {
			buf.append('[');
			buf.append(node.m_sMetaData);
			buf.append(']');
		}
		buf.append(":" + truncate(node.m_fLength, 6));
		return buf.toString();
	}
	
	
	void printUsageAndExit() {
		System.out.println("Usage: " + getClass().getName() + " [-n] <tree set file>\n");
		System.out.println("Prints tree from tree set in order of popularity of topology.\n" +
				"-n : only the most popular topology is shown in the output.\n\n" +
				"On stdout, it prints header lines and a table with individual heights and thetas, like this: " +
				"#Tree 1.  Frequency = 500. \n" +
				"# ((A[theta = theta0, height = 0],B[theta=theta1, height=0]]):[theta=theta2,height=h0],C:[theta=theta3,height=0]):[theta=theta4,height=height1]\n" +
				"nr        samplenr theta0    ...          theta4        height0        …        height1\n" +
				"1        120        0.2        ….        0.3           1.1                     2.0\n" +
				"2        180        0.1        ….        0.1           1.2                     2.1\n" +
				"3        210        0.3        ….        0.2           1.3                     2.0\n" +
				"...\n" +
				"500      121        0.3        ….        0.2           1.3                     2.0\n" + 
				"\n\n\n" +
				"On stderr, it prints out the mean tree for each topology, like this:\n" +
				"#Tree 0. 99.7003% (((0[mTheta=0.0851]:0.1596,1[mTheta=0.0976]:0.1596)[mTheta=0.0108]:0.2024,2[mTheta=0.0895]:0.3621)[mTheta=0.0108]:0.3167,3[mTheta=0.0486]:0.6788)[mTheta=0.0107]:0.0#Tree 1. 0.1998002% ((A[theta0]:height0,(B[theta1]:height1,C[theta2]:height2)[theta4]:height4)[theta5]:height5,D[theta3]:height3)[theta6]:height6\n" +
				"#Tree 1. 0.1998002% ((0[mTheta=0.0514]:1.4090,(1[mTheta=0.0938]:0.6033,2[mTheta=0.0821]:0.6033)[mTheta=0.0455]:0.8056)[mTheta=0.0171]:2.7178,3[mTheta=0.0066]:4.1268)[mTheta=0.0263]:0.0\n" +
				"#Tree 2. 0.0999001% (0[mTheta=0.0036]:0.0862,((1[mTheta=0.0200]:0.0189,2[mTheta=0.0037]:0.0189)[mTheta=0.0125]:0.0275,3[mTheta=0.0012]:0.0465)[mTheta=0.0071]:0.0397)[mTheta=0.0188]:0.0\n" +
				"");
		System.exit(0);
	}

	void parseArgs(String [] args) {
		int i = 0;
		try {
			while (i < args.length) {
				int iOld = i;
				if (i < args.length) {
					if (args[i].equals("")) {
						i += 1;
					} else if (args[i].equals("-n")) {
						m_bFirstTopologyOnly = true;
						i += 1;
					}
					if (i == iOld) {
						if (i == args.length-1) {
							m_sFileName = args[i];
							i++;
						} else {
							throw new Exception("Wrong argument");
						}
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error parsing command line arguments: " + Arrays.toString(args) + "\nArguments ignored\n\n");
			printUsageAndExit();
		}
		if (m_sFileName == null) {
			printUsageAndExit();
		}
	} // parseArgs
	

	public void analyse(String [] args) {
		parseArgs(args);
		try {
			TreeFileParser parser = new TreeFileParser(m_sLabels, null, null, 0);
			Node [] trees = parser.parseFile(m_sFileName);
			int nBase = (int) (trees.length * m_fBurnIn / 100.0);
			m_trees = new Node[trees.length - nBase];
			System.arraycopy(trees, nBase, m_trees, 0, m_trees.length);
		
			// count tree topologies
			// first step is find how many different topologies are present
			m_nTopology = new int[m_trees.length];
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
			m_nTopologies = map.size();
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
			m_fTreeWeight = new float[m_nTopologies];
			m_cTrees = new Node[m_nTopologies];
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
			
			produceOutput();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	void produceOutput() throws Exception {
		int nNodes = m_sLabels.size() * 2 -1;
		System.out.println("#nr coverage tree");
		int j = 0;
		for (int i = 0; i < m_nTopologies; i++) {
			System.out.println("#Tree " + i + ". " + m_fTreeWeight[i]*100 + "% " + getTopology(m_cTrees[i], m_sLabels));//m_cTrees[i].toString(sLabels));
			System.out.println("nr\t" + getHeader(m_cTrees[i]));
			initLists(nNodes);
			boolean bSameTree = true;
			while (bSameTree) {
				System.out.println(j + "\t" +getTreeData(m_trees[j]));
				j++;
				bSameTree = (j < m_trees.length) && (m_nTopology[j] == m_nTopology[j-1]);
			}
			calcMean(m_cTrees[i]);
			System.err.println("#Tree " + i + ". " + m_fTreeWeight[i]*100 + "% " + toString(m_cTrees[i], m_sLabels));
			if (m_bFirstTopologyOnly) {
				return;
			}
		}
	}
	
	public static void main(String [] args) {
		TreeSetAnalyser analyser = new TreeSetAnalyser();
		analyser.analyse(args);
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