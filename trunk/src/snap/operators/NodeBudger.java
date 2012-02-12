
/*
 * File NodeBudger.java
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
package snap.operators;



import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

@Description("Moves internal node height without changing the tree topology. " +
		"So the range is limited by the height of the parent node and the height " +
		"of the highest child.")
public class NodeBudger extends NodeSwapper {
	public Input<Double> m_pWindowSize = new Input<Double>("size", "Relative size of the window in which to move the root node");
	double m_fWindowSize;
	int m_nInternalNodeCount = -1;
	int m_nLeafNodeCount = -1;

	@Override
	public void initAndValidate() {
		m_nInternalNodeCount = m_pTree.get().getInternalNodeCount();
		m_nLeafNodeCount = m_pTree.get().getLeafNodeCount();
		m_fWindowSize = m_pWindowSize.get();
	}
	
	@Override
	public double proposal() {
		double hastingsRatio = 1.0;
		Tree tree = m_pTree.get(this);
		Node [] nodes = tree.getNodesAsArray();

		//Choose a random node internal node 
		// that is a node with number [m_nNodeCount/2 + 1, .... , m_nNodeCount - 1] 
		int whichNode = m_nLeafNodeCount + Randomizer.nextInt(m_nInternalNodeCount);
		Node p = nodes[whichNode];

		if (p.isRoot()){
//			System.err.println("*********************Root budger********************");
//			
//			System.exit(0);			
			// RRB: budging the root node leads to very long calculation times
			// so we reject its move. The root time can still be changed through
			// the scale operator, so the root time is not necessarily fixed.
			
			//DJB Given the problems with mixing heights, and the bugs fixed since RRB's comments, have decided
			// to implement a root move, using the description in Drummond et al, 2002.
			double beta = 1.0/m_fWindowSize;
			
			beta = 2.0; //TEMPORARY FIX!!!!
			
			if (beta < 1.0)
				beta = 1.0;
			
			
			if (beta==1.0) //No move possible, return a reject. 
				return Double.NEGATIVE_INFINITY;
			double maxc = p.getLeft().getHeight();
			if (p.getRight() != null) {
				maxc = Math.max(maxc, p.getRight().getHeight());
			}
			
			/**
			 Let h be the height, and maxc be the height of the nearest child. This move moves h to 
				h' = maxc + u * (h - maxc)
			 where u is uniform [1/beta, beta].
			 
			   Note that the reverse move is h = maxc + u'(h' - maxc) = maxc + u'(maxc + u*(h - maxc) - maxc) = maxc + u'*u*(h-maxc)
			 so that u' = 1/u. Hastings ration is therefore the Jacobian of
			 
			 u	h
			 0  -u^{-2}  
			 
			 which is 1/u.
			 
			 
			 **/
			
			double u = Randomizer.nextDouble()*(beta - 1.0/beta)+1.0/beta;
			double oldh = p.getHeight();
			double h2 = maxc + u*(p.getHeight() - maxc);
			p.setHeight(maxc + u*(p.getHeight() - maxc));
			//System.out.println("oldh = "+oldh+"\tbeta = "+beta+"\t u = "+u+"\t maxc = "+maxc+" \tp.getHeight = "+p.getHeight()+" proposed new height = "+h2 +"="+p.getHeight());
			
			//return Double.NEGATIVE_INFINITY;
			return -Math.log(u);
		}

		//Find shortest branch to any child.
		double minb = p.getLeft().getHeight();
		if (p.getRight() != null) {
			minb = Math.max(minb, p.getRight().getHeight());
		}

		double range = p.getParent().getHeight() - minb;

		double move = minb + Randomizer.nextDouble()*range;

		p.setHeight(move);

		return Math.log(hastingsRatio);
	}

	/** automatic parameter tuning **/
	@Override
	public void optimize(double logAlpha) {
		Double fDelta = calcDelta(logAlpha);
		fDelta += Math.log(m_fWindowSize);
		m_fWindowSize = Math.exp(fDelta);
		if (m_fWindowSize > 1.0) {
			m_fWindowSize = 1.0;
		}
		//System.err.println("NodeBudger " + m_fWindowSize);
    }
	
} // class NodeBudger
