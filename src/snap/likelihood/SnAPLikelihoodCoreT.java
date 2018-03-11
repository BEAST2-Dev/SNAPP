
/*
 * File SnAPLikelihoodCoreT.java
 *
 * Copyright (C) 2010 Remco Bouckaert, David Bryant remco@cs.auckland.ac.nz
 *
 * This file is part of SnAP.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * SnAP is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  SnAP is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SnAP; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */
package snap.likelihood;





import snap.Data;
import snap.NodeData;

import java.util.Arrays;
import java.util.concurrent.Future; 

import beast.app.BeastMCMC;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Node;



/** threaded version of SSSLikelihoodCore **/ 
public class SnAPLikelihoodCoreT  extends SnAPLikelihoodCore {
	SiteProbabilityCalculatorT m_siteProbabilityCalculatorT;
	public SnAPLikelihoodCoreT(Node root, Alignment data) {
		super(root, data);
		m_siteProbabilityCalculatorT = new SiteProbabilityCalculatorT();
	}
	
	static boolean stopRequested = false;
	Future<?> [] m_future;
	private class SSSRunnable implements Runnable {
		int m_iStart;
		int m_iStep;
		int m_iMax; 
		Data m_data;
		NodeData m_root;
		double m_u;
		double m_v;
		double m_rate;
		boolean m_bMutationOnlyAtRoot;
		boolean m_bHasDominantMarkers;
		boolean m_bUseCache;
		Double [] m_coalescenceRate;
		
	  SSSRunnable(int iStart, int iStep, int iMax, Data data, Double [] coalescenceRate, NodeData root, double u, double v, double rate, boolean bMutationOnlyAtRoot, boolean bHasDominantMarkers, boolean bUseCache) {
	    m_iStart = iStart;
	    m_iStep = iStep;
	    m_iMax = iMax;
	    m_data = data;
	    m_root = root;
	    m_u = u;
	    m_v = v;
	    m_rate = rate;
		m_bMutationOnlyAtRoot = bMutationOnlyAtRoot;
		m_bHasDominantMarkers = bHasDominantMarkers;
	    m_bUseCache = bUseCache;
	    m_coalescenceRate = coalescenceRate;
	  }
	  public void run() {
		  int iThread = m_iStart;
	    for (int id = m_iStart; id < m_iMax; id+= m_iStep) {
			//if (id>0 && id%100 == 0)
			//	System.err.print(id + " ");
			double siteL=0.0;

			try {
				int [] thisSite = m_data.getPattern(id);
				int [] lineageCounts = m_data.getPatternLineagCounts(id);
				//siteL =  SiteProbabilityCalculatorT.computeSiteLikelihood(m_root, m_u, m_v, thisSite, m_bUseCache, false, m_iStart);
				siteL =  m_siteProbabilityCalculatorT.computeSiteLikelihood(m_root, m_u, m_v, m_rate, m_coalescenceRate, thisSite, lineageCounts, m_bMutationOnlyAtRoot, m_bHasDominantMarkers, m_bUseCache, false, 0);

			}
			catch (Exception ex) {
				ex.printStackTrace();
				System.exit(1);
			}
//			if (siteL==0.0) {
//				return -10e100;
//			}
			patternProb[id] = siteL;
	    }
	  }
	}
	
	double [] patternProb;
	
	/**
	 Compute Likelihood of the allele counts
	 
	 @param root  The tree. Uses branch lengths and coalescenceRate values stored on this tree.
	 @param u  Mutation rate from red to green
	 @param v Mutation rate from green to red
	 @param sampleSize  Number of samples taken at each species (index by id field of the NodeData)
	 @param redCount  Each entry is a different marker... for each marker the number of red alleles in each species.
	 @param siteProbs  Vector of probabilities (logL) for each site.
	 * @throws Exception 
	 **/
	@Override
	public double [] computeLogLikelihood(NodeData root, double u, double v, double rate, 
			int [] sampleSizes, 
			Data data, 
			Double [] coalescenceRate,
			boolean bMutationOnlyAtRoot,	
			boolean bHasDominantMarkers,							  
			boolean bUseCache,
			boolean dprint /*= false*/) throws Exception
	{
		m_lineageCountCalculator.computeCountProbabilities(root,sampleSizes,coalescenceRate,bHasDominantMarkers,dprint);
		//dprint = true;
			
			//TODO: Partial subtree updates over all sites.
			
			
			//double forwardLogL = 0.0;
			int numPatterns = data.getPatternCount();

			//Temporarily store pattern probabilities... used for numerical checks.
			patternProb = new double[numPatterns];
			Arrays.fill(patternProb, -1);
			int nThreads = BeastMCMC.m_nThreads;
			m_siteProbabilityCalculatorT.clearCache(root.getNodeCount(), data.getMaxStateCount(), nThreads);

			m_future = new Future<?>[nThreads];
			for (int i = 0; i < nThreads; i++) {
				m_future[i] = BeastMCMC.g_exec.submit(new SSSRunnable(i, nThreads, numPatterns, data, coalescenceRate, root.copy(), u, v, rate, bMutationOnlyAtRoot, bHasDominantMarkers, bUseCache));
			}

			// correction for constant sites
//			int [] thisSite = new int[data.m_sTaxaNames.size()];
//			double P0 =  SiteProbabilityCalculator.computeSiteLikelihood(root,u,v,thisSite, false, false);
//			for (int i = 0; i < thisSite.length; i++) {
//				thisSite[i] = data.m_nStateCounts.elementAt(i);
//			}
//			double P1 =  SiteProbabilityCalculator.computeSiteLikelihood(root,u,v,thisSite, false, false);
//			forwardLogL-=(double) data.getSiteCount() * Math.log(1.0 - P0 - P1);
//			System.err.println(numPatterns + " " + forwardLogL);
			
			// wait for the other thread to finish
			for (int i = 0; i < nThreads; i++) {
				m_future[i].get();
			}

			return patternProb;
			
//			for(int id = 0; id < numPatterns-(bUsenNonPolymorphic ? 0 : 2); id++) {
//				forwardLogL+=(double)data.getPatternWeight(id) * Math.log(patternProb[id]);
//			}
//			if (!bUsenNonPolymorphic) {
//				// correction for constant sites
//				double P0 =  patternProb[numPatterns - 2];
//				double P1 =  patternProb[numPatterns -1 ];
//				forwardLogL-=(double) data.getSiteCount() * Math.log(1.0 - P0 - P1);
//			}
//			//System.err.println(numPatterns + " " + forwardLogL);
//			return forwardLogL;
	} // computeLogLikelihood
	

} // class TreeLikelihood
