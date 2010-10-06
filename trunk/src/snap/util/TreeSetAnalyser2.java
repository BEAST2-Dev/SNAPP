package snap.util;


import java.util.Arrays;
import java.util.Vector;

import snap.util.TreeFileParser;

public class TreeSetAnalyser2 extends TreeSetAnalyser {
	/** tree to compare the tree set with, if provided **/
	Node m_originalTree = null;

	void printUsageAndExit() {
		System.out.println("Usage: " + getClass().getName() + " [-tree <newick tree>] <tree set file>\n");
		System.out.println("Analyses the tree set, and compares with a tree if provided (e.g. the original used to simulate data from)\n" +
				"-tree <newick tree>: tree in newick format on command line\n" +
				"<tree set file>: name of the file containing tree set in NEXUS format");
		System.out.println("Outputs:\n"+
				"(1) Size of the 95% credible set.\n"+
				"(2) Probabilities of the top 20 (?) trees\n"+
				"(3) Whether or not the true tree is in the credible set."
				);
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
					} else if (args[i].equals("-tree")) {
						String sTree = args[i+1];
						Vector<String> sLabels = null;
						TreeFileParser parser = new TreeFileParser(sLabels, null, null, 0);
						m_originalTree = parser.parseNewick(sTree);
						i += 2;
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
	
	void produceOutput() {
			// find #topologies with 95% coverage
			double fSum = 0;
			int k = 0;
			while (fSum < 0.95) {
				fSum += m_fTreeWeight[k++];
			}
			System.out.println("95% HPD contains " + k +" topologies, out of a total of " + m_nTopologies + " topologies");
			
			// print top 20, or top 95% coverage tree topologies, whichever is smaller
			System.out.println("#nr coverage tree");
			for (int i = 0; i < Math.min(k, 20); i++) {
				System.out.println("Tree " + i + ". " + m_fTreeWeight[i]*100 + "% " + getTopology(m_cTrees[i], sLabels));
			}
			
			// check if standard tree is in credible set
			int iTopology = -1;
			if (m_originalTree == null) {
				// no tree provided, TODO print message perhaps?
			} else {
				String sTopology = getTopology(m_originalTree, sLabels);
				boolean bContained = false;
				for (int j = 0; j < k && !bContained; j++) {
					if (sTopology.equals(getTopology(m_cTrees[j], sLabels))) {
						bContained = true;
						iTopology = j;
					}
				}
				System.out.println("Original tree is " + (bContained?"part of":"not in") + " the 95% HPD");
			}
			
			// collect theta and height information of 
			if (iTopology < 0) {
				// cannot do this
			} else {
				int nNodes = sLabels.size() * 2 -1;
				System.out.println("#nr coverage tree");
				int j = 0;
				for (int i = 0; i < m_nTopologies; i++) {
					if (i == iTopology) {
						System.out.println("#Tree " + i + ". " + m_fTreeWeight[i]*100 + "% " + getTopology(m_cTrees[i], sLabels));//m_cTrees[i].toString(sLabels));
						System.out.println("nr\t" + getHeader(m_cTrees[i]));
					}
					initLists(nNodes);
					boolean bSameTree = true;
					while (bSameTree) {
						if (i == iTopology) {
							System.out.println(j + "\t" +getTreeData(m_trees[j]));
						}
						j++;
						bSameTree = (j < m_trees.length) && (m_nTopology[j] == m_nTopology[j-1]);
					}
					calcMean(m_cTrees[i]);
				}
			}
			
			
	}
	
	public static void main(String [] args) {
		TreeSetAnalyser2 analyser = new TreeSetAnalyser2();
		analyser.analyse(args);
	} // main
}

