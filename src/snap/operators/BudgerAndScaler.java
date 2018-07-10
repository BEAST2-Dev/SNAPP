
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
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import snap.GammaParameter;
import snap.NodeData;

@Description("Moves internal node height without changing the tree topology. " +
		"So the range is limited by the height of the parent node and the height " +
		"of the highest child. The gamma value is also scaled so that gamma*height is constant for all three branches")
public class BudgerAndScaler extends NodeSwapper {
	public Input<GammaParameter> m_gamma = new Input<GammaParameter>("gamma","array representation of gamma parameter", Validate.REQUIRED);
	int m_nNodeCount = -1;

	@Override
	public void initAndValidate() {
		m_nNodeCount = m_pTree.get().getNodeCount();
	}
	
	@Override
	public double proposal() {
		double hastingsRatio = 1.0;
		Tree tree = m_pTree.get(this);
		GammaParameter gamma = m_gamma.get(this);
		Node [] nodes = tree.getNodesAsArray();

		//Choose a random node internal node 
		int whichNode = m_nNodeCount/2 + 1 + Randomizer.nextInt(m_nNodeCount/2 - 1);
		NodeData p = (NodeData) nodes[whichNode];

		if (p.isRoot()){
			// RRB: budging the root node leads to very long calculation times
			// so we reject its move. The root time can still be changed through
			// the scale operator, so the root time is not necessarily fixed.
			return Double.NEGATIVE_INFINITY;
		}

		//Find shortest branch to any child.
		double minb = 10e10;
		minb = Math.max(p.getLeft().getHeight(), p.getRight().getHeight());

		double range = p.getParent().getHeight() - minb;

		/* Hastings ratio calculations
		 Variables are height h of node. g1, g2 and g3, gamma values of edges to left, right and parent.  u = uniform(h - minb,parent.height).
		 
		 If new height is h', we want
			(h' - left.height)*g1' = (h - left.height)*g1, so 
			g1' = [(h - left.height)/(h' - left.height)]*g1
			g2' = [(h - right.height)/(h' - right.height)]*g2
			g3' = [(h - par.height)/(h' - par.height)]*g3.
		 
		 (h,thetaL,thetaR,thetaP,u) -> (h+u,[(h - left.height)/(h + u - left.height)]*g1,[(h - right.height)/(h + u - right.height)]*g2, g3' = [(h - par.height)/(h+u - par.height)]*g3,-u)
		 
		 Jacobian is det of
		
		 1		0		0		0	1
		 ?		[(h - left.height)/(h + u - left.height)]	0	0	?
		 ?		0		[(h - right.height)/(h + u - right.height)]		0	?
		 ?		0		0	[(h - par.height)/(h+u - par.height)]		?
		 0		0		0		0		-1
		 
		 The ? indicating terms that don't contribute to the determinant. The Jacobian is therefore
		 
		- [(h - left.height)/(h + u - left.height)]*[(h - right.height)/(h + u - right.height)]* [(h - par.height)/(h+u - par.height)]
		 
		 */
		
		double h = p.getHeight();
		int iNode = p.getNr();
		int iLeft = p.getLeft().getNr();
		int iRight = p.getRight().getNr();
		NodeData left = (NodeData)p.getLeft(); 
		NodeData right = (NodeData)p.getRight();
		
		double g1 = gamma.getValue(iLeft);
		double h1 = left.getHeight();
		double g2 = gamma.getValue(iRight);
		double h2 = right.getHeight();

		double g3 = gamma.getValue(iNode);
		double h3 = p.getParent().getHeight();

		double u = Randomizer.nextDouble()*range;
		
		if (outsideBounds(((h3-h)/(h3-h-u) * g3), gamma) ||
			outsideBounds(((h-h1)/(h+u-h1) * g1), gamma) ||
			outsideBounds(((h-h2)/(h+u-h2) * g2), gamma)) {
			return Double.NEGATIVE_INFINITY;
		}
		
		gamma.setValue(iNode, ((h3-h)/(h3-h-u) * g3));
		gamma.setValue(iLeft, ((h-h1)/(h+u-h1) * g1));
		gamma.setValue(iRight, ((h-h2)/(h+u-h2) * g2));
		
		p.setHeight(h+u);
		
		hastingsRatio = (h-h1)*(h-h2)*(h3-h) / ( (h+u-h1)*(h+u-h2)*(h3-h-u));
		
		return Math.log(hastingsRatio);
	}

	private boolean outsideBounds(double d, RealParameter realParameter) {
		if (d < realParameter.getLower() || d > realParameter.getUpper()) {
			return true;
		}
		return false;
	}

	
} // class BudgerAndScaler
