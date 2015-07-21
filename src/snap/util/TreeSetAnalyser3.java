package snap.util;


//import java.io.PrintStream;
//import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import snap.util.TreeFileParser;

public class TreeSetAnalyser3 extends TreeSetAnalyser {
	/** tree to compare the tree set with, if provided **/
	String m_sOriginalTree = null;
	
	public static void printUsage() {
		System.out.println("Usage: " + TreeSetAnalyser3.class.getName() + " [-tree <newick tree>] [-b <burnin percentage>] <tree set file>\n");
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
	
//    final static int MAX_LAG = 2000;
//	/** values from which the ESS is calculated **/
//	List<Double> m_trace;
//	/** sum of trace, excluding burn-in **/
//	double m_fSum = 0;
//	/** keep track of sums of trace(i)*trace(i_+ lag) for all lags, excluding burn-in  **/
//    List<Double> m_fSquareLaggedSums;
//	/** return effective sample size of a sample 
//	 * @param trace:  values from which the ESS is calculated **/
//	public double ESS(List<Double> trace) {
//		m_fSum = 0;
//		m_trace = new ArrayList<Double>();
//		m_fSquareLaggedSums = new ArrayList<Double>();
//		double fESS = 0;
//		for (int i = 1; i < trace.size(); i++) {
//			fESS = incrementalESS(trace.get(i));
//		}
//		return fESS;
//	}
//	double incrementalESS(double fNewValue) {
//		m_trace.add(fNewValue);
//		m_fSum += fNewValue;
//		
//        final int nTotalSamples = m_trace.size();
//
//        // take no burn in, it is already taken out of the input
//        //final int iStart = 0;
//        final int nSamples = nTotalSamples - 0;
//        final int nMaxLag = Math.min(nSamples, MAX_LAG);
//
//        // calculate mean
//        final double fMean = m_fSum / nSamples;
//
//        while (m_fSquareLaggedSums.size()< nMaxLag) {
//        	m_fSquareLaggedSums.add(0.0);
//        }
//    
//        // calculate auto correlation for selected lag times
//        double[] fAutoCorrelation = new double[nMaxLag];
//        // fSum1 = \sum_{iStart ... nTotalSamples-iLag-1} trace
//    	double fSum1 = m_fSum;
//        // fSum1 = \sum_{iStart+iLag ... nTotalSamples-1} trace
//    	double fSum2 = m_fSum;
//        for (int iLag = 0; iLag < nMaxLag; iLag++) {
//            m_fSquareLaggedSums.set(iLag, m_fSquareLaggedSums.get(iLag) + m_trace.get(nTotalSamples - iLag - 1) * m_trace.get(nTotalSamples - 1));
//            // The following line is the same approximation as in Tracer 
//            // (valid since fMean *(nSamples - iLag), fSum1, and fSum2 are approximately the same)
//            // though a more accurate estimate would be
//            // fAutoCorrelation[iLag] = m_fSquareLaggedSums.get(iLag) - fSum1 * fSum2
//            fAutoCorrelation[iLag] = m_fSquareLaggedSums.get(iLag) - (fSum1 + fSum2) * fMean + fMean * fMean * (nSamples - iLag);
//            fAutoCorrelation[iLag] /= ((double) (nSamples - iLag));
//        	fSum1 -= m_trace.get(nTotalSamples - 1 - iLag);
//        	fSum2 -= m_trace.get(0 + iLag);
//        }
//
//        double integralOfACFunctionTimes2 = 0.0;
//        for (int iLag = 0; iLag < nMaxLag; iLag++) {
//            if (iLag == 0) {
//                integralOfACFunctionTimes2 = fAutoCorrelation[0];
//            } else if (iLag % 2 == 0) {
//                // fancy stopping criterion - see main comment
//                if (fAutoCorrelation[iLag - 1] + fAutoCorrelation[iLag] > 0) {
//                    integralOfACFunctionTimes2 += 2.0 * (fAutoCorrelation[iLag - 1] + fAutoCorrelation[iLag]);
//                } else {
//                    // stop
//                    break;
//                }
//            }
//        }
//
//        // auto correlation time
//        final double fACT = integralOfACFunctionTimes2 / fAutoCorrelation[0];
//
//        // effective sample size
//        final double fESS = nSamples / fACT;
//        return fESS;
//    } // ESS
	
	@Override
	public boolean parseArgs(String [] args) {
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
		System.out.println("95% HPD contains " + k +" topologies,... out of a total of " + m_nTopologies + " topologies");
		
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
			
			//Write true value, posterior mean value, posterior std dev for each height and each theta.
			
			//Heights first.
			System.out.println("\tPProb\t"+m_fTreeWeight[iTopology]);
			System.out.println("\tHeight \toriginal\tmean\tstderr of mean");
			for(int node = m_sLabels.size(); node<nNodes;node++) {
				List<Double> heights = m_heights[node];
				int n = heights.size()-1;
								
				//System.out.print("\t"+heights.get(1));
				
				Double sum = 0.0;
				Double sumSquared = 0.0;
				for(int i=1;i<=n;i++) {
					Double r = heights.get(i);
					sum+=r;
					sumSquared+=r*r;
				}
				
				Double mean = sum/n;
				Double var = sumSquared/(n - 1.0)  - sum*sum/(n*(n-1.0));
				
				double fESS = beast.core.util.ESS.calcESS(heights);
				
				//System.err.print("\t"+heights.get(0)+"\t"+mean+"\t"+Math.sqrt(var));
				System.out.println("\tNode"+node+"\t"+heights.get(0)+"\t"+mean+"\t"+Math.sqrt(var/fESS));
			}
			
			//Now Thetas
			System.out.println("\tTheta \toriginal\tmean\tstderr of mean");
			for(int node = 0; node<nNodes;node++) {
				List<Double> thetas = m_thetas[node];
				int n = thetas.size()-1;
				
				//System.out.print("\t"+thetas.get(1));

				
				
				Double sum = 0.0;
				Double sumSquared = 0.0;
				for(int i=1;i<=n;i++) {
					Double r = thetas.get(i);
					sum+=r;
					sumSquared+=r*r;
				}
				
				Double mean = sum/n;
				Double var = sumSquared/(n - 1.0)  - sum*sum/(n*(n-1.0));
				
				double fESS = beast.core.util.ESS.calcESS(thetas);
				
				//System.err.print("\t"+thetas.get(0)+"\t"+mean+"\t"+Math.sqrt(var));
				System.out.println("\ttheta"+node+"\t"+thetas.get(0)+"\t"+mean+"\t"+Math.sqrt(var/fESS));// + "\t" + var);
			}
			
			System.err.println();
			System.out.println();
			
			
		}
		
		
	}
	
	public static void main(String [] args) {
		TreeSetAnalyser3 analyser = new TreeSetAnalyser3();
		if (analyser.parseArgs(args)) {
			analyser.run();
		}
	} // main
}

