
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

import snap.NodeData;
import snap.distribution.ChiSquareNoncentralDist;

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
        "Thetas are represented by the coalescenceRate parameter where values are theta=2/coalescenceRate")
public class SnAPPrior extends Distribution {
    public Input<RealParameter> m_pAlpha = new Input<RealParameter>("alpha", "prior parameter -- see docs for details", Validate.REQUIRED);
    public Input<RealParameter> m_pBeta = new Input<RealParameter>("beta", "prior parameter -- see docs for details", Validate.REQUIRED);
    public Input<RealParameter> m_pKappa = new Input<RealParameter>("kappa", "prior parameter -- see docs for details");
    public Input<RealParameter> m_pCoalescenceRate = new Input<RealParameter>("coalescenceRate", "Populations sizes for the nodes in the tree", Validate.REQUIRED);
    public Input<RealParameter> m_pLambda = new Input<RealParameter>("lambda", "parameter for Yule birth process");//, Validate.REQUIRED);
    public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations"); //, Validate.REQUIRED);

    @Override
    public void initAndValidate() throws Exception {
    	if (m_pKappa.get() == null) {
    		System.err.println("WARNING: kappa parameter not set for SnAPPrior. using default value of 1.0");
    		m_pKappa.setValue(new RealParameter("1.0"), this);
    	}
    }


    @Override
    public double calculateLogP() throws Exception {
        logP = 0.0;

        double alpha = m_pAlpha.get().getValue();
        double beta = m_pBeta.get().getValue();
        double kappa = m_pKappa.get().getValue();
        
        Tree tree = m_pTree.get();
        double heightsum = tree.getRoot().getHeight();
        heightsum += heightSum(tree.getRoot());
        int nspecies = (tree.getNodeCount() + 1) / 2;
        double lambda = m_pLambda.get().getValue();

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
        RealParameter coalescenceRate = m_pCoalescenceRate.get();

		
		int PRIORCHOICE = 2;
		
		if (PRIORCHOICE == 0) {
			//Assume independent gamma distributions for thetas.
			
			//We assume that 2/r has a gamma(alpha,beta) distribution. That means that r has density proportional to
			// 1/(r^2)  * GAMMA(2/r|alpha,beta)
			//which has log (alpha - 1.0)*Math.log(2.0/r) - (beta *(2.0/ r)) - 2*log(r), which in turn simplifies to the expr. below (w/ consts)
		
			for (int iNode = 0; iNode < coalescenceRate.getDimension(); iNode++) {
				double r = coalescenceRate.getValue(iNode);
				logP += -(alpha + 1.0)*Math.log(r) - 2.0* beta / r;
			}
		} else if (PRIORCHOICE == 1) {
			//Assume independent inverse gamma distributions for thetas
			
			//We assume that (2/r) has an inverse gamma (alpha,beta) distribution. That means that r has density proportional to
			// 1/(r^2) * INVGAMMA(2/r|alpha,beta)
			//which has logarithm:
			// -(alpha+1) log(2/r) - beta*r/2 - 2 log(r) = C + (alpha-1) log(r) - beta*r / 2.
			for (int iNode = 0; iNode < coalescenceRate.getDimension(); iNode++) {
				double r = coalescenceRate.getValue(iNode);
				logP += (alpha - 1.0)*Math.log(r) - 0.5* beta * r;
			}
		} else if (PRIORCHOICE == 2) {
	        //> let x be the root.
	        //> r = rate for node x.
	        double r = ((NodeData)tree.getRoot()).coalescenceRate();
	        //> logP += (alpha - 1.0)*Math.log(r) - 0.5* beta * r;
	        logP += (alpha - 1.0)*Math.log(r) - 0.5* beta * r;

	        
	        //> for all nodes x in a pre-order traversal
	        //> 	let y be parent of node x. 
	        //> 	let t be length of branch connecting x with y
	        //> 	r = rate for node x. 
	        //> 	r0 = rate for node y.
	        //> 	df = 2*alpha
	        //> 	nc = 2*(2/r0) * beta * exp(-kappa * t) / (1 - exp(-kappa*t))
	        //> 	p = ChiSquareNoncentralDist(df,nc) -> density (2/r)
	        //> 	logP +=log(p) - 2 * log (r)         [this second part is the Jacobian for the transform r <-> theta]
	        //> end        

	        for (int iNode = 0; iNode < tree.getNodeCount(); iNode++) {
	        	NodeData node = (NodeData) tree.getNode(iNode);
	        	if (!node.isRoot()) {
	        		NodeData parent = (NodeData) node.getParent();
	        		double t = parent.getHeight() - node.getHeight();
	        		r = node.coalescenceRate();
	        		double r0 = parent.coalescenceRate();
	        		double df = 2 * alpha;
	        		double nc = 2*(2.0/r0) * beta * Math.exp(-kappa * t) / (1 - Math.exp(-kappa*t));
	        		double p = ChiSquareNoncentralDist.density(df, nc, 2/r);
	        		logP += Math.log(p) - 2.0 * Math.log(r);
	        	}
	        }
	        
			
		} else {
			//Assume that rate has uniform distribution on [[0,1000]
			for (int iNode = 0; iNode < coalescenceRate.getDimension(); iNode++) {
				double r = coalescenceRate.getValue(iNode);
				if (r>1000.0 || r<0.0)
					return -Double.NEGATIVE_INFINITY;
			}
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
