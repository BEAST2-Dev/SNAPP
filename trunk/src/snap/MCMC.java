
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


import beast.core.*;
import beast.evolution.tree.Tree;
import beast.util.*;

@Description("Allow sampling from the prior.")
public class MCMC extends beast.core.MCMC {
	public Input<Integer> m_oStateBurnIn = new Input<Integer>("stateBurnin", "Number of burn in samples taken on the prior only to determine initial state", 0);
	public Input<Distribution> m_stateDistribution = new Input<Distribution>("stateDistribution", "Uncertainty from which the initial state is sampled for 'stateBurnin' of samples. Must be specified if stateBurnin is larger than zero.");
	public Input<Integer> m_killAfterXSeconds = new Input<Integer>("killAfter", "Number of seconds after which the job is killed. When negative (default) this is ignored", -1);
	
	@Override
	public void initAndValidate() throws Exception {
		super.initAndValidate();
		if (m_oStateBurnIn.get() > 0) {
			if (m_stateDistribution.get() == null) {
				throw new Exception("stateBurnin is larger than zero, but stateUncertainty not specified");
			}
		}
	} // init


	@Override
	public void run() throws Exception {
        // initialises log so that log file headers are written, etc.
        for (Logger log : m_loggers.get()) {
            log.init();
        }
    	// set up state (again). Other plugins may have manipulated the
    	// StateNodes, e.g. set up bounds or dimensions
    	state.initAndValidate();
    	// also, initialise state with the file name to store and set-up whether to resume from file
    	state.setStateFileName(m_sStateFile + "." + Randomizer.getSeed());
		int nBurnIn = m_oBurnIn.get();
		int nChainLength = m_oChainLength.get();
		int nStateBurnin = m_oStateBurnIn.get();
        if (m_bRestoreFromFile) {
        	state.restoreFromFile();
        	nBurnIn = 0;
        	//prepare();
        }
		long tStart = System.currentTimeMillis();
		state.setEverythingDirty(true);
		
		int nKillAfterXSeconds = m_killAfterXSeconds.get();

		System.err.println("Start state:");
		System.err.println(state.toString());
		
		// do the sampling
		double logAlpha = 0;		
		
		
		// Sample initial state from 'prior' = stateUncertainty
		// For this to work, stateBurnin must be specified as the number of samples to be 
		// taken from the 'prior', and the stateUncertainty must contain the 'prior'.
		if (!m_bRestoreFromFile && nStateBurnin > 0) {
			System.err.println("Sampling state from prior for " + nStateBurnin + " samples...");
			//prepare();
			
	        state.store(-1);
	        state.setEverythingDirty(true);
	        state.checkCalculationNodesDirtiness();
			double fOldLogLikelihood = m_stateDistribution.get().calculateLogP();
	        state.setEverythingDirty(false);
	        state.acceptCalculationNodes();
			
			for (int iSample = -nBurnIn - nStateBurnin; iSample <= -nBurnIn; iSample++) {
				
				//State proposedState = state.copy();
	        	state.store(iSample);
	        	posteriorInput.get().store();
				Operator operator = operatorSchedule.selectOperator();
				double fLogHastingsRatio = operator.proposal();
				if (fLogHastingsRatio != Double.NEGATIVE_INFINITY) {
					state.storeCalculationNodes();
	                state.checkCalculationNodesDirtiness();
					//m_stateDistribution.get().store();
					
					//prepare();
					double fNewLogLikelihood = m_stateDistribution.get().calculateLogP();
					
					logAlpha = fNewLogLikelihood -fOldLogLikelihood + fLogHastingsRatio; //CHECK HASTINGS
		            if (logAlpha>=0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
						// accept
						fOldLogLikelihood = fNewLogLikelihood;
	                    state.acceptCalculationNodes();
						//state = proposedState;
						state.setEverythingDirty(false);
					} else {
						//m_stateDistribution.get().restore();
			        	//posteriorInput.get().restore();
			        	state.restore();
			        	state.restoreCalculationNodes();
					}
		            
		            Tree tree = (Tree)state.stateNodeInput.get().get(0);
		            if (tree.getRoot().getHeight() > 2) {
		            	int k = 3;
		            	k++;
		            }
		            
				} else {
	                state.restore();
				}
				
			}
			// store initial state
			System.err.println("Post Start state:");
			//System.err.println(m_state.toString(m_data.m_sTaxaNames));
			System.err.println(state.toString());
		}

		
		// Go into the main loop
		boolean bDebug = true;
		state.setEverythingDirty(true);
		Distribution posterior = posteriorInput.get();
		
		
		double fOldLogLikelihood = robustlyCalcPosterior(posterior);//posterior.calculateLogP();
		
		System.err.println("Start likelihood: = " + fOldLogLikelihood);
		
//		int iGammaParameter = -1;
//		for (int i = 0; i < state.getNrOfStateNodes(); i++) {
//			if (state.getStateNode(i) instanceof GammaParameter) {
//				iGammaParameter = i; 
//			}
//		}
		
		
		boolean bStayAlive = true;
		for (int iSample = -nBurnIn; iSample <= nChainLength && bStayAlive; iSample++) {
			
			//State proposedState = state.copy();
        	state.store(iSample);

			Operator operator = operatorSchedule.selectOperator();
			//System.err.println(operator.getClass().getName() + " " + operator.getID());
			double fLogHastingsRatio = operator.proposal();
			if (fLogHastingsRatio != Double.NEGATIVE_INFINITY) {
				state.storeCalculationNodes();
                state.checkCalculationNodesDirtiness();
				
				double fNewLogLikelihood = posterior.calculateLogP();
				
				logAlpha = fNewLogLikelihood -fOldLogLikelihood + fLogHastingsRatio; //CHECK HASTINGS
	            if (!Double.isNaN(fNewLogLikelihood) && (logAlpha>=0 || Randomizer.nextDouble() < Math.exp(logAlpha))) {
					// accept
					fOldLogLikelihood = fNewLogLikelihood;
					state.setEverythingDirty(false);
                    state.acceptCalculationNodes();
					if (iSample>=0) {
						operator.accept();
					}
					//System.err.println("ACCEPT");
				} else {
					// reject
					if (iSample>=0) {
						operator.reject();
					}
					state.restore();
                    state.setEverythingDirty(false);
                    state.restoreCalculationNodes();
                    // needed for logging tree theta's correctly
//                    ((GammaParameter)state.getStateNode(iGammaParameter)).prepare();
					//System.err.println("REJECT");
				}
			} else {
				// operation failed
				if (iSample>0) {
					operator.reject();
				}
				state.restore();
			}
            log(iSample);
            
			
			if (bDebug || iSample % 10000 == 0) {
				state.store(iSample);
				state.setEverythingDirty(true);
				//System.err.println(m_state.toString());
				double fLogLikelihood = posterior.calculateLogP();
				state.setEverythingDirty(false);
				if (Math.abs(fLogLikelihood - fOldLogLikelihood) > 1e-6) {
					throw new Exception("Likelihood incorrectly calculated: " + fOldLogLikelihood + " != " + fLogLikelihood);
				}
				if (iSample > NR_OF_DEBUG_SAMPLES) {
					bDebug = false;
				}
			} else {
				operator.optimize(logAlpha);
			}
			if (nKillAfterXSeconds > 0) {
				long tEnd = System.currentTimeMillis();
				if ((tEnd - tStart)/1000 >= nKillAfterXSeconds) {
					System.out.println("\n\nTime is up! Going to terminate now.\n");
					bStayAlive = false;
				}
			}
		}
		operatorSchedule.showOperatorRates(System.out);
		long tEnd = System.currentTimeMillis();
		System.out.println("Total calculation time: " + (tEnd - tStart)/1000.0 + " seconds");
		close();

	
		state.setEverythingDirty(true);
        System.err.println("End likelihood: " + fOldLogLikelihood);
        System.err.println("End state:");
		System.err.println(state.toString());
		state.storeToFile(nChainLength);
	} // run;

} // class MCMC
