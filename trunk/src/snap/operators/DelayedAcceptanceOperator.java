package snap.operators;

import java.util.List;

import snap.likelihood.SnapSubstitutionModel;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.util.CompoundDistribution;
import beast.core.Operator;
import beast.core.StateNode;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.alignment.distance.HammingDistance;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;


@Description("An operator that uses an approximate likelihood and filters out proposales with low acceptance " +
		"based on the approximate likleihood")
public class DelayedAcceptanceOperator extends Operator {
	public Input<Operator> operatorInput = new Input<Operator>("operator", "Operator for proposing moves", Validate.REQUIRED);
	public Input<State> stateInput = new Input<State>("state", "state object for which we do proposals", Validate.REQUIRED);

	// TODO: do we need the complete posterior? If so, we have to distinguish between prior and likelihood
	public Input<Distribution> posteriorInput = new Input<Distribution>("posterior", "Distribution to be approximated", Validate.REQUIRED);

	
    public Input<SiteModelInterface> siteModelInput = new Input<SiteModelInterface>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
	
	
	Operator operator;
	Distribution prior = null;
	// data and tree used in approx treelikelihood
	Alignment data = null;
	Tree tree = null;
	
	// TODO: take site model in account in approximation
	SiteModel.Base siteModel;
	SnapSubstitutionModel substitutionmodel;
	
    State state;
    
    // empirical distance and variance between taxa
    double [][] distance;
    double var;

	
	public void initAndValidate() {
		operator = operatorInput.get();
    	siteModel = (SiteModel.Base) siteModelInput.get();
    	substitutionmodel = ((SnapSubstitutionModel)siteModel.substModelInput.get());
		
        if (posteriorInput.get() instanceof CompoundDistribution) {
            final CompoundDistribution posterior = (CompoundDistribution) posteriorInput.get();
            final List<Distribution> distrs = posterior.pDistributions.get();
            final int nDistr = distrs.size();
            for (int i = 0; i < nDistr; i++) {
                final Distribution distr = distrs.get(i);
                final String sID = distr.getID();
                if (sID != null && sID.equals("prior")) {
                	prior = distr;
                } else if (sID != null && sID.equals("likelihood")) {
                    if (distr instanceof CompoundDistribution) {
                    	final List<Distribution> distrs2 = ((CompoundDistribution)distr).pDistributions.get();
                    	if (distrs2.size() != 1) {
                    		throw new RuntimeException("DelayedAcceptanceOperator: Expected only one distribution in likelihood");
                    	}
                    	Distribution distr2 = distrs2.get(0);
                    	if (!(distr2 instanceof TreeLikelihood)) {
                    		TreeLikelihood tl = (TreeLikelihood) distr2;
                    		data = tl.dataInput.get();
                    		tree = (Tree) tl.treeInput.get();
                    	}
                    } else {
                		throw new RuntimeException("DelayedAcceptanceOperator: Expected likelihood to be a CompoundDistribution");
                    }
	            } else {
	        		throw new RuntimeException("DelayedAcceptanceOperator: Expected likelihood or prior in posterior and nothing else");
	            }
            }
        } else {
            throw new RuntimeException("DelayedAcceptanceOperator: Expected a CompoundDistribution as posterior input");
        }

        
        if (prior == null) {
        	throw new RuntimeException("DelayedAcceptanceOperator: could not identify prior in posterior input");
        }
        if (data == null || tree == null) {
        	throw new RuntimeException("DelayedAcceptanceOperator: could not identify data or tree in treelikelihood in posterior input");
        }
        
        
        calcDistanceAndVariance();
    }


    private void calcDistanceAndVariance() {
    	// set up distance matrix
    	
    	// TODO: verify this is the correct distance
		HammingDistance d = new HammingDistance();
		d.setPatterns(data);
		int nrOfTaxa = data.getNrTaxa();
		distance = new double[nrOfTaxa][nrOfTaxa];
		for (int i = 0; i < nrOfTaxa; i++) {
			for (int j = i+ 1; j < nrOfTaxa; j++) {
				distance[i][j] = d.pairwiseDistance(i,  j);
				distance[j][i] = distance[i][j];
			}
		}
		// TODO: normalise?
		
		// estimate variance
		// TODO: should be variance for each entry [i,j]?
		double v = 0;
		double mean = 0;
		for (int i = 0; i < nrOfTaxa; i++) {
			for (int j = i+ 1; j < nrOfTaxa; j++) {
				mean += distance[i][j];
				v += distance[i][j] * distance[i][j];
			}
		}
		// TODO: should we take the diagonal in account (as is implemented) or not?
		v = v * 2 / (nrOfTaxa * nrOfTaxa);
		mean = mean * 2 / (nrOfTaxa * nrOfTaxa);
		var = v - mean * mean; 
	}

	@Override
    public double proposal()  {
    	try {
	    	double oldApproxLogLikelihood = evaluate();
	    	double logHastingsRatio = operator.proposal();
	    	if (logHastingsRatio == Double.NEGATIVE_INFINITY) {
	    		// instant reject
	    		return Double.NEGATIVE_INFINITY;
	    	}

            state.storeCalculationNodes();
            state.checkCalculationNodesDirtiness();
	    	double newApproxLogLikelihood = evaluate();

	    	double logAlpha = newApproxLogLikelihood - oldApproxLogLikelihood + logHastingsRatio; //CHECK HASTINGS
	        if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
	        	// accept
	        	
	        	// TODO: do we need to restore 
	        	//state.restore();
	            // TODO: can we get the state nr?
	            //state.store(-1);

	        	// TODO: this is probably not the correct correction of the HR
	            logHastingsRatio += newApproxLogLikelihood - oldApproxLogLikelihood;
	        	
	        } else {
	        	// reject;
	            state.restore();
	            // TODO: can we get the state nr?
	            state.store(-1);
	        	return Double.NEGATIVE_INFINITY;
	        }
	    	return logHastingsRatio;
    	} catch (Exception e) {
    		System.err.println("DelayedAcceptanceOperator.proposal " + e.getMessage());
            state.restore();
            // TODO: can we get the state nr?
            state.store(-1);
    		return Double.NEGATIVE_INFINITY;
    	}
    }

	private double evaluate() throws Exception {
		double logP = prior.calculateLogP();
		logP  += approxLikelihood();
		return logP;
	}
	
	
	/* 
	 * calculate approximate treelikelihood for tree & data 
	 * */
	private double approxLikelihood() {
    	Double [] coalescenceRate = substitutionmodel.m_pCoalescenceRate.get().getValues();
    	double u = substitutionmodel.m_pU.get().getValue();
    	double v  = substitutionmodel.m_pV.get().getValue();

    	
    	// TODO: do the real work here
    	double mu[][] = calcMu();
    	
    	// TODO: fill in this constant
    	double K = 0;
    	
    	double approxL = 0;
    	int nrOfTaxa = distance.length;
		for (int i = 0; i < nrOfTaxa; i++) {
			for (int j = i+ 1; j < nrOfTaxa; j++) {
				// TODO: should be variance for each entry [i,j]?
				approxL += (distance[i][j] - mu[i][j]) * (distance[i][j] - mu[i][j]) / var;   
			}
		}

		return K  + approxL;
	}
	
	private double[][] calcMu() {
		// TODO calculate estimates of distance between taxa based on the 
		// tree and other parameters
		return null;
	}

	// BEAST infrastructure stuff
	@Override
	public void accept() {
		// TODO do we want operator tuning? If not, comment out this line as well as in reject() below 
		operator.accept();
	}
	
	@Override
	public void reject() {
		operator.reject();
	}
	
	@Override
	public List<StateNode> listStateNodes() throws Exception {
		return operator.listStateNodes();
	}

} // class DelayedAcceptanceOperator