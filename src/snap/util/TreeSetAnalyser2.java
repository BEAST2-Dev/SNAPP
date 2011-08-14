package snap.util;


import java.util.Arrays;
import java.util.List;

import snap.util.TreeFileParser;

public class TreeSetAnalyser2 extends TreeSetAnalyser {
	/** tree to compare the tree set with, if provided **/
	String m_sOriginalTree = null;

	public static void printUsage() {
		System.out.println("Usage: " + TreeSetAnalyser2.class.getName() + " [-tree <newick tree>] [-b <burnin percentage>] <tree set file>\n");
		System.out.println("Analyses the tree set, and compares with a tree if provided (e.g. the original used to simulate data from)\n" +
				"-tree <newick tree>: tree in newick format on command line\n" +
				"-b <burnin percentage>: percentage (so a number between 0 and 100) of tree at the start that are discarded, default 10%\n" +
				"<tree set file>: name of the file containing tree set in NEXUS format");
		System.out.println("Outputs:\n"+
				"(1) Size of the 95% credible set.\n"+
				"(2) Probabilities of the top 20 (?) trees\n"+
				"(3) Whether or not the true tree is in the credible set."
				);
	}

	@Override
	boolean parseArgs(String [] args) {
		int i = 0;
		try {
			while (i < args.length) {
				int iOld = i;
				if (i < args.length) {
					if (args[i].equals("")) {
						i += 1;
					} else if (args[i].equals("-tree")) {
						m_sOriginalTree = args[i+1];
						i += 2;
					} else if (args[i].equals("-b")) {
						m_fBurnIn = Double.parseDouble(args[i+1]);
						m_fBurnIn = Math.max(m_fBurnIn,0);
						m_fBurnIn = Math.min(m_fBurnIn,100);
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
			printUsage();
			return false;
		}
		if (m_sFileName == null) {
			printUsage();
			return false;
		}
		return true;
	} // parseArgs
	
	String getShortTopology(Node node, List<String> sLabels) {
		if (node.isLeaf()) {
			return sLabels.get(node.getNr()); 
		} else {
			return "(" + getShortTopology(node.m_left, sLabels) + "," + getShortTopology(node.m_right, sLabels)  
			+")"; 
		}
	}

		
		
	@Override
	void produceOutput() throws Exception {
			// find #topologies with 95% coverage
			double fSum = 0;
			int k = 0;
			while (fSum < 0.95) {
				fSum += m_fTreeWeight[k++];
			}
			System.out.println("95% HPD contains " + k +" topologies, out of a total of " + m_nTopologies + " topologies");
			
			//Standard error is condensed info only.
			System.err.print(m_sFileName+"\t"+k+"\t"+m_nTopologies);
		
			// print top 20, or top 95% coverage tree topologies, whichever is smaller
			System.out.println("#nr coverage tree");
			for (int i = 0; i < Math.min(k, 20); i++) {
				System.out.println("Tree " + (i+1) + ": " + m_fTreeWeight[i]*100 + "% " + getShortTopology(m_cTrees[i], m_sLabels));
			}
			
			// check if standard tree is in credible set
			int iTopology = -1;
			Node originalTree = null;
			if (m_sOriginalTree == null) {
				// no tree provided, TODO print message perhaps?
			} else {
				TreeFileParser parser = new TreeFileParser(m_sLabels, null, null, 0);
				originalTree = parser.parseNewick(m_sOriginalTree);
				originalTree.sort();
				originalTree.labelInternalNodes(m_sLabels.size());
				
				String sTopology = getShortTopology(originalTree, m_sLabels);
				boolean bContained = false;
				for (int j = 0; j < k && !bContained; j++) {
					if (sTopology.equals(getShortTopology(m_cTrees[j], m_sLabels))) {
						bContained = true;
						iTopology = j;
					}
				}
				System.out.println("Original tree " + sTopology + "is " + (bContained?"part of":"not in") + " the 95% HPD");
				System.err.print("\t"+(bContained?"Y":"N"));
			}
			
			// collect theta and height information of trees with same topology as 'original tree' 
			if (iTopology < 0) {
				System.err.println();
			} else {
				int nNodes = m_sLabels.size() * 2 -1;
				initLists(nNodes);
				System.out.println("#nr coverage tree");
				System.out.println("#Original tree " + getTopology(originalTree, m_sLabels));
				System.out.println("Original\t" +getTreeData(originalTree));
				int j = 0;
				for (int i = 0; i < m_nTopologies; i++) {
					if (i == iTopology) {
						System.out.println("#Tree " + i + ". " + m_fTreeWeight[i]*100 + "% " + getTopology(m_cTrees[i], m_sLabels));//m_cTrees[i].toString(sLabels));
						System.out.println("nr\t" + getHeader(m_cTrees[i]));
					}
					//initLists(nNodes);
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
				// RRB: at this point m_heights and m_thetas (both List<Double>[]) contain
				// arrays of height and theta information. The first entry is for the original
				// tree, so m_height.get(0)[0] is the height of node 0 in the original tree
				// and entries 1 and upward are from the sample.
				// Taxa are numbered 0,...,n and internal nodes n+1,...,2n-1 with node 2n-1 the root. 
				// TODO: do something useful with this...
				
				//Compute posterior bias and variance of root height and root population size.
				
				List<Double> root_heights = m_heights[m_cTrees[iTopology].getNr()];
				List<Double> root_thetas = m_thetas[m_cTrees[iTopology].getNr()];

				
				int ntrees = root_heights.size(); 
				
				Double sum = 0.0;
				Double sumSquared = 0.0;
				for(int i=1;i<ntrees;i++) {
					Double r = root_heights.get(i);
					sum+=r;
					sumSquared+=r*r;
				}
				System.out.println("Root Height | Root Height Bias | Root Height StdDev\n");
				Double mean = sum/(ntrees-1.0);
				Double bias = mean - root_heights.get(0);
				Double var = sumSquared/(ntrees-1.0) - mean*mean;
				System.out.println(""+mean+"\t"+bias+"\t"+Math.sqrt(var));
				//System.err.print("\t"+root_heights.get(0)+"\t"+bias+"\t"+Math.sqrt(var));

				sum = 0.0;
				sumSquared = 0.0;
				for(int i=1;i<ntrees;i++) {
					Double r = root_thetas.get(i);
					sum+=r;
					sumSquared+=r*r;
				}
				System.out.println("Root Theta|Root Theta Bias | Root Theta StdDev\n");
				mean = sum/(ntrees-1.0);
				bias = mean - root_thetas.get(0);
				var = sumSquared/(ntrees-1.0) - (mean*mean);
				System.out.println(""+mean + "\t"+bias+"\t"+Math.sqrt(var));
				//System.err.println("\t"+root_thetas.get(0) + "\t"+bias+"\t"+Math.sqrt(var));

			}
			
			
	}
	
	public static void main(String [] args) {
		TreeSetAnalyser2 analyser = new TreeSetAnalyser2();
		if (analyser.parseArgs(args)) {
			analyser.run();
		}
	} // main
}

