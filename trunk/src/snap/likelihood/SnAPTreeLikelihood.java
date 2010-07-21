
/*
 * File SnAPTreeLikelihood.java
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



import java.util.List;
import java.util.Random;

import beast.app.BeastMCMC;
import beast.evolution.alignment.Alignment;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Distribution;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

import snap.NodeData;
import snap.likelihood.SnAPLikelihoodCore;


@Description("Implements a tree Likelihood Function for Single Site Sorted-sequences on a tree.") 
@Citation("David Bryant, Remco Bouckaert, Noah Rosenberg. Inferring species trees directly from SNP and AFLP data: full coalescent analysis without those pesky gene trees. arXiv:0910.4193v1. http://arxiv.org/abs/0910.4193")
public class SnAPTreeLikelihood extends Distribution {

	Alignment m_data;
	int[] m_nSampleSizes;
	SnAPLikelihoodCore m_core;
	
	public Input<Alignment> m_pData = new Input<Alignment>("data", "set of alignments");
	public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations");
	public Input<RealParameter> m_pU = new Input<RealParameter>("mutationRateU", "mutation rate from red to green?");
	public Input<RealParameter> m_pV = new Input<RealParameter>("mutationRateV", "mutation rate from green to red?");
	
	
    /**
     * Constructor.
     */
	public SnAPTreeLikelihood() {
    }
    
    
    //public SSSTreeLikelihood(Data patternList, State state) {
    @Override
    public void initAndValidate() {
    	m_data = m_pData.get();
    	if ( BeastMCMC.m_nThreads == 1) {
    		m_core = new SnAPLikelihoodCore(m_pTree.get().getRoot(), m_pData.get());
    	} else {
    		m_core = new SnAPLikelihoodCoreT(m_pTree.get().getRoot(), m_pData.get());
    	}
    	Integer [] nSampleSizes = m_data.m_nStateCounts.toArray(new Integer[0]);
    	m_nSampleSizes = new int[nSampleSizes.length];
    	for (int i = 0; i < nSampleSizes.length; i++) {
    		m_nSampleSizes[i] = nSampleSizes[i];
    	}
    	
//    	try {
//    		m_nU = state.getParameterIndex("u");
//    		m_nV = state.getParameterIndex("v");
//    	} catch (Exception e) {
//			e.printStackTrace();
//		}
    }

   int m_nTreeID = -1;
    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLogP() {
    	try {
	    	NodeData root = (NodeData) m_pTree.get().getRoot();
	    	double u = m_pU.get().getValue();
	    	double v  = m_pV.get().getValue();
			boolean useCache = true;
			boolean dprint = false;
			logP = m_core.computeLogLikelihood(root, u, v, 
	    			m_nSampleSizes, 
	    			m_data, 
	    			useCache,
	    			dprint /*= false*/);
			return logP;
    	} catch (Exception e) {
			e.printStackTrace();
			return 0;
		}
    } // calculateLogLikelihood

	@Override
	public void store(int nSample) {
    	super.store(nSample);
    	m_core.m_bReuseCache = true;
    }

	@Override
    public void restore(int nSample) {
    	super.restore(nSample);
    	m_core.m_bReuseCache = false;
    }


	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
	}    


} // class SSSTreeLikelihood
