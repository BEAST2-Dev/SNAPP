
/*
 * File RateMixer.java
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
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

@Description("Moves length of branch and coalescence rate of branch in the opposit direction.")
public class RateMixer extends Operator {

	public Input<Double> m_pScale = new Input<Double>("scaleFactors", "scaling factor: larger means more bold proposals");
	public Input<RealParameter> m_coalescenceRate = new Input<RealParameter>("coalescenceRate", "population sizes");
	public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations");
	
	@Override
	public void initAndValidate() {
		m_fMixScale = m_pScale.get();
	}

	double m_fMixScale;

	@Override
	public double proposal() { // throws Exception {
		RealParameter coalescenceRate = m_coalescenceRate.get(this);

		double scale = (m_fMixScale + (Randomizer.nextDouble() * ((1.0 / m_fMixScale) - m_fMixScale)));

		
		for (int i = 0; i < coalescenceRate.getDimension(); i++) {
			coalescenceRate.setValue(i, coalescenceRate.getValue(i) * scale);
		}

		Tree tree = m_pTree.get(this);
		int nInternalNodes = 0;
		try {
			nInternalNodes = tree.scale(1/scale);
		} catch (Exception e) {
			return Double.NEGATIVE_INFINITY;
		}
		
		
		
		//For the Hastings ratio. If there are n leaves, then there are n-1 heights and 2n-1 branches (including one above root).
		// Hastings ratio is therefore (1/u)^(n-1) * u^(2n-1) * u^{-2} = u^(-n+1+2n-1-2} = u^{n-2}, which is not the same as below! 
		//Here n = nInternalNodes+1;
		
		
		return Math.log(scale) * (nInternalNodes-1);
	}

	/** automatic parameter tuning **/
	@Override
	public void optimize(double logAlpha) {
		Double fDelta = calcDelta(logAlpha);
		fDelta += Math.log(m_fMixScale);
		m_fMixScale = Math.exp(fDelta);
    }
	
	@Override
    public double getCoercableParameterValue() {
    	return m_fMixScale;
    }

	@Override
    public void setCoercableParameterValue(final double value) {
		m_fMixScale = value;
    }

} // class RateMixer
