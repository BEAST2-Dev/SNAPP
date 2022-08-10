
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


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("Scales single value in gamma parameter.")
public class GammaMover extends Operator {
	public Input<RealParameter> m_coalescenceRate = new Input<RealParameter>("coalescenceRate", "population sizes");
	public Input<Double> m_pScale = new Input<Double>("scale", "scale of move");

	double m_fScale;

	@Override
	public void initAndValidate() {
		m_fScale = m_pScale.get();
	}
	
	@Override
	public double proposal() {
		RealParameter coalescenceRate = m_coalescenceRate.get();
		int whichNode = Randomizer.nextInt(coalescenceRate.getDimension());
		
		double scale = Math.exp(m_fScale*(2.0*Randomizer.nextDouble() - 1.0));
		double newValue = coalescenceRate.getValue(whichNode)*scale;
		if (newValue < coalescenceRate.getLower() || newValue > coalescenceRate.getUpper()) {
			return Double.NEGATIVE_INFINITY;
		}

		coalescenceRate.setValue(whichNode, newValue);
		return Math.log(scale);
	}

	/** automatic parameter tuning **/
	@Override
	public void optimize(double logAlpha) {
		Double fDelta = calcDelta(logAlpha);
		fDelta += Math.log(m_fScale);
		m_fScale = Math.exp(fDelta);
    }
	
    @Override
    public double getCoercableParameterValue() {
    	return m_fScale;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
		m_fScale = value;
    }
} // class GammaMover 
