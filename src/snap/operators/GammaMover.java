
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


import snap.GammaParameter;
import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.util.Randomizer;

@Description("Scales single value in gamma parameter.")
public class GammaMover extends Operator {
	public Input<GammaParameter> m_pGamma = new Input<GammaParameter>("gamma", "population sizes");
	public Input<Double> m_pScaleGamma = new Input<Double>("pGammaMove", "scale of move");

	double m_fScale;

	@Override
	public void initAndValidate() {
		m_fScale = m_pScaleGamma.get();
	}
	
	@Override
	public double proposal() {
		GammaParameter gamma = m_pGamma.get(this);
		int whichNode = Randomizer.nextInt(gamma.getDimension());
		
		double scale = Math.exp(m_fScale*(2.0*Randomizer.nextDouble() - 1.0));
		gamma.setValue(whichNode, gamma.getValue(whichNode)*scale);
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
