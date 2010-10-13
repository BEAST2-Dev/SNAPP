
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
	public Input<Double> m_pWindowSize = new Input<Double>("size", "Relative size of the window in which to move the node");
	double m_fWindowSize;
	int m_nNodeCount = -1;

	@Override
	public void initAndValidate() {
		m_nNodeCount = m_pTree.get().getNodeCount();
		m_fWindowSize = 1;//m_pWindowSize.get();
	}
	
	@Override
	public double proposal() {
		double hastingsRatio = 1.0;
		Tree tree = m_pTree.get(this);
		Node [] nodes = tree.getNodesAsArray();

		//Choose a random node internal node 
		int whichNode = m_nNodeCount/2 + 1 + Randomizer.nextInt(m_nNodeCount/2 - 1);
		Node p = nodes[whichNode];

		if (p.isRoot()){
			// RRB: budging the root node leads to very long calculation times
			// so we reject its move. The root time can still be changed through
			// the scale operator, so the root time is not necessarily fixed.
			return Double.NEGATIVE_INFINITY;
		}

		//Find shortest branch to any child.
		double minb = 10e10;
		minb = Math.max(p.m_left.getHeight(), p.m_right.getHeight());

		double range = p.getParent().getHeight() - minb;

		double move = minb + m_fWindowSize * Randomizer.nextDouble()*range;

		p.setHeight(move);

		return Math.log(hastingsRatio);
	}

	/** automatic parameter tuning **/
	@Override
	public void optimize(double logAlpha) {
		Double fDelta = calcDelta(logAlpha);
		fDelta += Math.log(m_fWindowSize);
		m_fWindowSize = Math.exp(fDelta);
		if (m_fWindowSize > 1) {
			m_fWindowSize = 1;
		}
    }
	
} // class NodeBudger
