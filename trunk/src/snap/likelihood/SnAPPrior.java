
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

import java.util.Iterator;
import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Input;
import beast.core.Distribution;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.tree.Node;

@Description("Standard prior for SnAP analysis, consisting of a Yule prior on the tree " +
        "(parameterized by lambda) " +
        "and gamma distribution over the theta values " +
        "(with parameters alpha and beta). " +
        "Thetas are represented by the gamma parameter where values are theta=2/gamma")
public class SnAPPrior extends Distribution {
    public Input<RealParameter> m_pAlpha = new Input<RealParameter>("alpha", "prior parameter -- see docs for details", Validate.REQUIRED);
    public Input<RealParameter> m_pBeta = new Input<RealParameter>("beta", "prior parameter -- see docs for details", Validate.REQUIRED);
    public Input<RealParameter> m_pLambda = new Input<RealParameter>("lambda", "parameter for Yule birth process", Validate.REQUIRED);
    public Input<RealParameter> m_pGamma = new Input<RealParameter>("gamma", "Populations sizes for the nodes in the tree", Validate.REQUIRED);
    public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    	// nothing to initialise
    }


    @Override
    public double calculateLogP() throws Exception {
        logP = 0.0;

        Tree tree = m_pTree.get();
        double heightsum = tree.getRoot().getHeight();
        heightsum += heightSum(tree.getRoot());
        int nspecies = (tree.getNodeCount() + 1) / 2;
        double lambda = m_pLambda.get().getValue();
        double alpha = m_pAlpha.get().getValue();
        double beta = m_pBeta.get().getValue();

        double mu = 0.0; //Death process

        if (mu==0.0) {
        	//Yule Model
            logP += (nspecies-1)*Math.log(lambda) - lambda*heightsum;
        } else {
            //Birth death model.    See Thompson 1975, pg 56
            List<Double> allHeights = getSortedHeights(tree.getRoot());
            Iterator<Double> p = allHeights.iterator();
            double x1 = p.next();
            double p0n= p0(x1,lambda,mu);

            for(int n=2; n<nspecies;n++) {
                double xn = p.next();
                logP += Math.log(mu*p1(xn,lambda,mu)/p0n);
            }

        }
        //Gamma values in tree
        RealParameter gamma = m_pGamma.get();
        //double [] gamma = state.getParameter(m_pGamma).getValues();
        //	We assume that theta has a gamma (alpha,beta) distribution, so that
        //the gamma parameter has 2/gamma(alpha,beta) distribution
		//Because of the transform of variables, we need to include the appropriate
		//Jacobian.... that is, if f(theta) is the GAMMA distribution, we want
		// f(2/gamma) / (gamma^2). This has been independently tested in mathlab.
		
        for (int iNode = 0; iNode < gamma.getDimension(); iNode++) {
			double g = gamma.getValue(iNode);
			logP += -(alpha+1.0)*Math.log(g) - beta*(2.0/g);
			
            //double x = 2.0/gamma.getValue(iNode);
            //logP += (alpha - 1.0)*Math.log(x) - (beta * x);
        }
        return logP;
    } // calculateLogLikelihood

    double heightSum(Node node) throws Exception {
        if (node.isLeaf()) {
            return 0;
        } else {
            double h = node.getHeight();
            h += heightSum(node.m_left);
            if (node.m_right != null) {
            	h += heightSum(node.m_right);
            }
            return h;
        }
    } // heightSum

    //Returns a list of branching times in the tree, sorted in an decreasing sequence. First one is
    //the height of the mrca of the tree.
    List<Double> getSortedHeights(Node node) throws Exception {
        return null;
    }

    private double p0(double t, double lambda, double mu) {
        return mu*(1-Math.exp(-(lambda-mu)*t))/(lambda - mu*Math.exp(-(lambda-mu)*t));
    }

    private double p1(double t, double lambda, double mu) {
        double denominator = (lambda - mu*Math.exp(-(lambda-mu)*t));
        return (lambda - mu)*(lambda - mu) * Math.exp(-(lambda-mu)*t) / (denominator * denominator);
    }

	@Override public List<String> getArguments() {return null;}
	@Override public List<String> getConditions() {return null;}
	@Override public void sample(State state, Random random) {};
} // class SSSPrior 
