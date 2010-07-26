
/*
 * File SnAPPrior.java
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

import beast.core.Description;
import beast.core.Input;
import beast.core.Distribution;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.tree.Node;

@Description("Standard prior for SnAP analysis, consisting of a Yule prior on the tree " +
		"(parameterized by lambda) " +
		"and gamma distribution over the theta values " +
		"(with parameters alpha and beta). " +
		"Thetas are represented by the gamma parameter where values are theta=2/gamma")
public class SnAPPrior extends Distribution {
//	public Input<Parameter> m_pU = new Input<Parameter>("mutationRateU", "mutation rate from red to green?");
//	public Input<Parameter> m_pV = new Input<Parameter>("mutationRateV", "mutation rate from green to red?");
	public Input<RealParameter> m_pAlpha = new Input<RealParameter>("alpha", "prior parameter -- see docs for details");
	public Input<RealParameter> m_pBeta = new Input<RealParameter>("beta", "prior parameter -- see docs for details");
	public Input<RealParameter> m_pLambda = new Input<RealParameter>("lambda", "parameter for Yule birth process");
	public Input<RealParameter> m_pGamma = new Input<RealParameter>("gamma", "Populations sizes for the nodes in the tree");
	public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations");
	
	public SnAPPrior() {
	}
	
	public void initAndValidate() {
	}
	

	@Override
	public double calculateLogP() throws Exception {
		logP = 0.0;

		//Mutation rates
		//	uniform pi_0 on [0,1]. Rate equals one.
//		if (state.getValue(m_pU) < 0.0 || state.getValue(m_pV) < 0.0 )
//			return Double.NEGATIVE_INFINITY; //zero prior probability.


		//branch lengths / tree height
		//Assume a yule prior with given birthrate lambda.
		//computeHeights(state.tree);
		//state.m_trees[m_nTreeID].calculateHeightsFromLengths();
		Tree tree = (Tree) m_pTree.get();
		double heightsum = tree.getRoot().getHeight(); 
		heightsum += heightSum(tree.getRoot()); 

		int nspecies = (tree.getNodeCount() + 1) / 2;
		double lambda = m_pLambda.get().getValue();
		double alpha = m_pAlpha.get().getValue();
		double beta = m_pBeta.get().getValue();
		// note: nspecies-2, not nspecies-1
		// SHOULD BE -1 not -2?!?
		logP += (nspecies-1)*Math.log(lambda) - lambda*heightsum;

		//Gamma values in tree
		RealParameter gamma = m_pGamma.get();
		//double [] gamma = state.getParameter(m_pGamma).getValues();
		//	We assume that theta has a gamma (alpha,beta) distribution, so that
		//the gamma parameter has 2/gamma(alpha,beta) distribution
		for (int iNode = 0; iNode < gamma.getDimension(); iNode++) {
			double x = 2.0/gamma.getValue(iNode);
			logP += (alpha - 1.0)*Math.log(x) - (beta * x);
		}
		return logP;
	} // calculateLogLikelihood

	double heightSum(Node node) throws Exception {
		if (node.isLeaf()) {
			return 0;
		} else {
			double h = node.getHeight(); 
			h +=
				heightSum(node.m_left) +
				heightSum(node.m_right);
			return h;
		}
	} // heightSum


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

	
} // class SSSPrior 
