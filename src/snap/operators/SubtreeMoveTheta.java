
/*
 * File GammaMoveSubtree.java
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


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("Scales all theta values in a subtree by the same amount")
public class SubtreeMoveTheta extends Operator {
	public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations");
	public Input<RealParameter> m_coalescenceRate = new Input<RealParameter>("coalescenceRate", "population sizes");
	public Input<Double> m_pScale = new Input<Double>("scale", "scale of move");

	double m_fScale;

	@Override
	public void initAndValidate() {
		m_fScale = m_pScale.get();
	}
	
	@Override
	public double proposal() {

		Tree tree = m_pTree.get(); //Species tree
		Node [] nodes = tree.getNodesAsArray();
		int nnodes = nodes.length; //Number of nodes in the species tree
		int subtreeRootId = Randomizer.nextInt(nnodes);
		Node subtreeRoot = nodes[subtreeRootId];
		
		RealParameter coalescenceRates = m_coalescenceRate.get();
		double scale = Math.exp(m_fScale*(2.0*Randomizer.nextDouble() - 1.0));
		
		int numNodes = scaleSubtreeRates(subtreeRoot, coalescenceRates,scale);
		
		return numNodes * Math.log(scale); //Hastings ratio
	}
 
 	private int scaleSubtreeRates(Node v, RealParameter coalescenceRates, double scale) {
 		int nNodes = 1;
 		if (!v.isLeaf()) {
 			nNodes+=scaleSubtreeRates(v.getLeft(),coalescenceRates,scale);
 			nNodes+=scaleSubtreeRates(v.getRight(),coalescenceRates,scale);
 		}
 		coalescenceRates.setValue(v.getNr(), coalescenceRates.getValue(v.getNr())*scale);
 		return nNodes;
 	}
 	
 	
	/** automatic parameter tuning **/
	@Override
	public void optimize(double logAlpha) {
		Double fDelta = calcDelta(logAlpha);
		fDelta += Math.log(m_fScale);
		m_fScale = Math.exp(fDelta);
    }
	
} // class GammaMover 
