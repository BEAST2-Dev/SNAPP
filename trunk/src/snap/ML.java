
/*
 * File MCMC.java
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
package snap;


import java.util.ArrayList;
import java.util.List;

import beast.core.*;

@Description("Maximum likelihood search by hill climbing.")
public class ML extends beast.core.MCMC {
	
	public Input<Integer> m_oStateBurnIn = new Input<Integer>("stateBurnin", "Number of burn in samples taken on the prior only to determine initial state", new Integer(0));
	public Input<Distribution> m_stateDistribution = new Input<Distribution>("stateDistribution", "Uncertainty from which the initial state is sampled for 'stateBurnin' of samples. Must be specified if stateBurnin is larger than zero.");
	
	List<Integer> m_gammas;
	
	@Override
	public void initAndValidate() throws Exception {
		super.initAndValidate();
		if (m_oStateBurnIn.get() > 0) {
			if (m_stateDistribution.get() == null) {
				throw new Exception("stateBurnin is larger than zero, but stateUncertainty not specified");
			}
		}
		
		m_gammas = new ArrayList<Integer>();
		for (int i = 0 ; i < state.getNrOfStateNodes(); i++) {
			StateNode stateNode = state.getStateNode(i);
			if (stateNode instanceof GammaParameter) {
				m_gammas.add(stateNode.getIndex());
			}
		}
	} // init


	public void run() throws Exception {
		long tStart = System.currentTimeMillis();
		state.setEverythingDirty(true);
		
		int nBurnIn = m_oBurnIn.get();
		int nChainLength = m_oChainLength.get();
		//int nStateBurnin = m_oStateBurnIn.get();

		System.err.println("Start state:");
		System.err.println(state.toString());
		
		// do the sampling
		double logAlpha = 0;

//		// Sample initial state from 'prior' = stateUncertainty
//		// For this to work, stateBurnin must be specified as the number of samples to be 
//		// taken from the 'prior', and the stateUncertainty must contain the 'prior'.
//		if (nStateBurnin > 0) {
//			System.err.println("Sampling state from prior for " + nStateBurnin + " samples...");
//			prepare();
//			double fOldLogLikelihood = m_stateDistribution.get().calculateLogP();
//			for (int iSample = -nBurnIn - nStateBurnin; iSample <= -nBurnIn; iSample++) {
//				
//				//State proposedState = state.copy();
//	        	state.store();
//	        	posteriorInput.get().store(iSample);
//				Operator operator = operatorSet.selectOperator();
//				double fLogHastingsRatio = operator.proposal();
//				if (fLogHastingsRatio != Double.NEGATIVE_INFINITY) {
//					storeCachables(iSample);
//					m_stateDistribution.get().store(iSample);
//					
//					prepare();
//					double fNewLogLikelihood = m_stateDistribution.get().calculateLogP();
//					logAlpha = fNewLogLikelihood -fOldLogLikelihood + fLogHastingsRatio; //CHECK HASTINGS
//		            if (logAlpha>=0) { // || Randomizer.nextDouble() < Math.exp(logAlpha)) {
//						// accept
//						fOldLogLikelihood = fNewLogLikelihood;
//						//state = proposedState;
//						state.setDirty(false);
//					} else {
//						m_stateDistribution.get().restore(iSample);
//			        	state.restore();
//			        	posteriorInput.get().restore(iSample);
//	                    restoreCachables(iSample);
//					}
//				}
//			}
//			// store initial state
//			System.err.println("Post Start state:");
//			//System.err.println(m_state.toString(m_data.m_sTaxaNames));
//			System.err.println(state.toString());
//		}
//
//		
		// Go into the main loop
		boolean bDebug = false;
		state.setEverythingDirty(true);
		double fOldLogLikelihood = calc();
		System.err.println("Start likelihood: = " + fOldLogLikelihood);
		for (int iSample = -nBurnIn; iSample <= nChainLength; iSample++) {
			
			//State proposedState = state.copy();
        	state.store(iSample);
			Operator operator = operatorSet.selectOperator();
			if (iSample == 24) {
				int h = 3;
				h++;
				//proposedState.makeDirty(State.IS_GORED);
			}
			double fLogHastingsRatio = operator.proposal();
			if (fLogHastingsRatio != Double.NEGATIVE_INFINITY) {
				//System.out.print("store ");
				state.storeCalculationNodes();
//				proposedState.makeDirty(State.IS_GORED);
				//System.out.print(operator.getName()+ "\n");
				//System.err.println(proposedState.toString());
				//if (bDebug) {
					//System.out.print(operator.getName()+ "\n");
					//System.err.println(proposedState.toString());
					//state.validate();
				//}			
				
				
				
				double fNewLogLikelihood = calc();
				logAlpha = fNewLogLikelihood -fOldLogLikelihood;// + fLogHastingsRatio; //CHECK HASTINGS
	            if (logAlpha>=0) {// || Randomizer.nextDouble() < Math.exp(logAlpha)) {
					// accept
					fOldLogLikelihood = fNewLogLikelihood;
					state.setEverythingDirty(false);
					if (iSample>=0) {
						operator.accept();
					}
				} else {
					// reject
					if (iSample>=0) {
						operator.reject();
					}
					state.restore();
                    state.setEverythingDirty(false);
                    state.restoreCalculationNodes();
					//System.out.println("restore ");
				}
			} else {
				// operation failed
				if (iSample>0) {
					operator.reject();
				}
			}
			log(iSample);
			
			if (bDebug) {
				//state.validate();
				state.setEverythingDirty(true);
				//System.err.println(m_state.toString());
				double fLogLikelihood = calc();
				if (Math.abs(fLogLikelihood - fOldLogLikelihood) > 1e-10) {
					throw new Exception("Likelihood incorrectly calculated: " + fOldLogLikelihood + " != " + fLogLikelihood);
				}
				if (iSample > NR_OF_DEBUG_SAMPLES) {
					bDebug = false;
				}
			} else {
				operator.optimize(logAlpha);
			}
		}
		operatorSet.showOperatorRates(System.out);
		long tEnd = System.currentTimeMillis();
		System.out.println("Total calculation time: " + (tEnd - tStart)/1000.0 + " seconds");
		close();
	} // run;

	/** calculate log likelihood for posterior **/
	protected double calc() throws Exception {
		prepare();
		// recalculates all likelihoods from scratch
		double fLogP = 0;
		fLogP = posteriorInput.get().calculateLogP();
		return fLogP;
	} // calc

	/** synchronise gamma values in tree with parameter **/
	void prepare() {
		for (Integer iGamma : m_gammas) {
			GammaParameter gamma = (GammaParameter) state.getStateNode(iGamma);
			gamma.prepare();
		}
	} // prepare
} // class MCMC
