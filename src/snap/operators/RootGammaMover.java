
/*
 * File GammaMover.java
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
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

@Description("Scales value in gamma parameter associated with root node.")
public class RootGammaMover extends Operator {
	public Input<RealParameter> m_coalescenceRate = new Input<RealParameter>("coalescenceRate", "population sizes", Validate.REQUIRED);
	public Input<Double> m_pScale = new Input<Double>("scale", "scale of move");
	public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations", Validate.REQUIRED);

	double m_fScale;
	Tree tree;

	@Override
	public void initAndValidate() {
		m_fScale = m_pScale.get();
		tree = m_pTree.get();
	}
	
	@Override
	public double proposal() {
		
		RealParameter coalescenceRate = m_coalescenceRate.get(this);

		// We want whichNode to be equal to the root node.	
		int whichNode = tree.getRoot().getNr(); 
		
		double scale = Math.exp(m_fScale*(2.0*Randomizer.nextDouble() - 1.0));
		coalescenceRate.setValue(whichNode, coalescenceRate.getValue(whichNode)*scale);
		return Math.log(scale);
	}

	/** automatic parameter tuning **/
	@Override
	public void optimize(double logAlpha) {
		Double fDelta = calcDelta(logAlpha);
		fDelta += Math.log(m_fScale);
		m_fScale = Math.exp(fDelta);
    }
	
} // class GammaMover 
