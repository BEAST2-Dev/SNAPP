
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


import snap.GammaParameter;
import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

@Description("Moves length of branch and gamma of branch in the oposite direction.")
public class RateMixer extends Operator {

	public Input<Double> m_pScale = new Input<Double>("scaleFactors", "scaling factor: larger means more bold proposals");
	public Input<GammaParameter> m_pGamma = new Input<GammaParameter>("gamma", "population sizes");
	public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations");
	
	@Override
	public void initAndValidate() {
		m_fMixGamma = m_pScale.get();
	}

	double m_fMixGamma;

	@Override
	public double proposal() { // throws Exception {
		GammaParameter gamma = m_pGamma.get(this);

		//double scale = Math.exp(m_fMixGamma*(2.0*Randomizer.nextDouble() - 1.0));
        double scale = (m_fMixGamma + (Randomizer.nextDouble() * ((1.0 / m_fMixGamma) - m_fMixGamma)));
		//state.mulValues(scale, gamma);
		for (int i = 0; i < gamma.getDimension(); i++) {
			gamma.setValue(i, gamma.getValue(i) * scale);
		}
		//gamma.mulValues(scale);
		Tree tree = m_pTree.get(this);
		int nInternalNodes = 0;
		try {
			nInternalNodes = tree.scale(1/scale);
		} catch (Exception e) {
			return Double.NEGATIVE_INFINITY;
		}

		return Math.log(scale) * (-2 * gamma.getDimension());// - nInternalNodes - 0);
	}

	/** automatic parameter tuning **/
	@Override
	public void optimize(double logAlpha) {
		Double fDelta = calcDelta(logAlpha);
		fDelta += Math.log(m_fMixGamma);
		m_fMixGamma = Math.exp(fDelta);
    }
}
