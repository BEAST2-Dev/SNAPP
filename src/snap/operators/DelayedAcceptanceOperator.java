package snap.operators;



import java.util.List;

import snap.likelihood.SnAPTreeLikelihood;
import snap.likelihood.SnapSubstitutionModel;

import beast.core.Description;
import beast.core.Input;
import beast.core.OperatorSchedule;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.StateNode;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

//Things to measure:
//	time of full calculation over all samples
//  time of approx calculation
//  rate of rejection rate
//  final acceptance rate




@Description("An operator that uses an approximate likelihood and filters out proposales with low acceptance " +
		"based on the approximate likleihood")
public class DelayedAcceptanceOperator extends Operator {
	public Input<Operator> operatorInput = new Input<Operator>("operator", "Operator for proposing moves", Validate.REQUIRED);

	public Input<Double> tuningInput = new Input<Double>("tuning", "tuning parameter for approx likelihood");

	
	public Input<State> stateInput = new Input<State>("state", "state object for which we do proposals", Validate.REQUIRED);

//	public Input<TreeInterface> treeInput = new Input<TreeInterface>("tree", "tree used for approximate likelihood", Validate.REQUIRED);
//	public Input<Alignment> dataInput = new Input<Alignment>("data", "alignment used for approximate likelihood", Validate.REQUIRED);
//    public Input<SiteModelInterface> siteModelInput = new Input<SiteModelInterface>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
    public Input<SnAPTreeLikelihood> treeLikelihoodInput = new Input<SnAPTreeLikelihood>("treelikelihood", "SNAPP tree likelihood for the tree", Validate.REQUIRED);
	
    
    public Input<ApproximateLikelihood> approxLiklihoodInput = new Input<ApproximateLikelihood>("approxLikelihood", "Approximate SNAPP tree likelihood for the tree", Validate.REQUIRED);
    
    
    
	
	Operator operator;
	
	// TODO: take site model in account in approximation?
	// TODO: take clock model in account in approximation?
	SiteModel.Base siteModel;
	SnapSubstitutionModel substitutionmodel;
	SnAPTreeLikelihood treelikelihood;
	ApproximateLikelihood approxLiklihood;
	
    State state;
    
	TreeInterface tree = null;

	
	public void initAndValidate() {
		operator = operatorInput.get();
		treelikelihood = treeLikelihoodInput.get();
    	siteModel = (SiteModel.Base) treelikelihood.siteModelInput.get();
    	substitutionmodel = ((SnapSubstitutionModel)siteModel.substModelInput.get());
		
    	tree = treelikelihood.treeInput.get();
    	state = stateInput.get();
    	
    	approxLiklihood = approxLiklihoodInput.get();

/*
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
*/
        
        
        
    }


	@Override
    public double proposal()  {

		double probVariableSites = treelikelihood.getProbVariableSites();
		
    	try {
    		Node oldRoot = tree.getRoot().copy();
        	Double [] oldCoalescenceRate = substitutionmodel.m_pCoalescenceRate.get().getValues();
        	double oldU = substitutionmodel.m_pU.get().getValue();
        	double oldV  = substitutionmodel.m_pV.get().getValue();

        	double oldApproxLogLikelihood = approxLiklihood.approxPosterior(oldRoot, oldCoalescenceRate, oldU, oldV);
	    	double oldPrior = approxLiklihood.getPriorValue();
	    	
	    	double logHastingsRatio = operator.proposal();
            
			Node newRoot = tree.getRoot();
	    	Double [] newCoalescenceRate = substitutionmodel.m_pCoalescenceRate.get().getValues();
	    	double newU = substitutionmodel.m_pU.get().getValue();
	    	double newV = substitutionmodel.m_pV.get().getValue();

	    	
	    	// could skip till after checking logHR == -infinity
	    	// but since the proposal can return -infinity in two different cases
	    	// (if slave-operator return -infinity, OR if proposal is rejected)
	    	// that would make it difficult to deal with in the main MCMC loop
	    	// to distinguish between those two cases.
	    	state.storeCalculationNodes();
            state.checkCalculationNodesDirtiness();

            if (logHastingsRatio == Double.NEGATIVE_INFINITY) {
	    		// instant reject
	    		// need to store state, so it is restored properly
	    		return Double.NEGATIVE_INFINITY;
	    	}

	    	double newApproxLogLikelihood = approxLiklihood.approxPosterior(newRoot, newCoalescenceRate, newU, newV);
	    	double newPrior = approxLiklihood.getPriorValue();

	    	double logAlpha = newApproxLogLikelihood + newPrior - oldApproxLogLikelihood - oldPrior + logHastingsRatio; //CHECK HASTINGS
	        if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
	        	// accept
	        	
	        	// TODO: do we need to restore 
	        	//state.restore();
	            // TODO: can we get the state nr?
	            //state.store(-1);

	        	// reset the HR
	        	// no need to worry about HR of slave-operator
	        	// note HR contains priors, that cancel out in MCMC loop
	        	
	        	if (probVariableSites == 1.0) {
		            logHastingsRatio = oldApproxLogLikelihood - newApproxLogLikelihood;
	        	} else {
	        		double operatorLogHastingsRatio = logHastingsRatio;
	        		// we need to correct for non-constant site probability in the newly proposed site
		    		if (logAlpha < 0) {
		    			logHastingsRatio = -logAlpha;
		    		} else {
		    			logHastingsRatio = 0;	
		    		}
//System.err.print(operator.getName() + " " + oldApproxLogLikelihood + " " +  newApproxLogLikelihood + " " + oldPrior + " " + newPrior + " ");			    	
		    		probVariableSites = treelikelihood.getNewProbVariableSites();
		    		// no need to recalculate: no parameter changed in between
			    	//oldApproxLogLikelihood = oldPrior + approxLiklihood.approxLikelihood(oldRoot, oldCoalescenceRate, oldU, oldV);
			    	//newApproxLogLikelihood = newPrior + approxLiklihood.approxLikelihood(newRoot, newCoalescenceRate, newU, newV);
			    	double ratio = oldApproxLogLikelihood - newApproxLogLikelihood - operatorLogHastingsRatio;
			    	if (ratio < 0) {
			    		logHastingsRatio += ratio;
			    	}
//System.err.println(" " + logHastingsRatio + " ");			    	
 			    	return operatorLogHastingsRatio + logHastingsRatio;
	        	}
	        	
	        	
	        } else {
	        	// reject;
	            //state.restore();
	            // TODO: can we get the state nr?
	            //state.store(-1);
	        	return Double.NEGATIVE_INFINITY;
	        }
	    	return logHastingsRatio;
    	} catch (Exception e) {
    		e.printStackTrace();
    		System.err.println("DelayedAcceptanceOperator.proposal " + e.getMessage());
            state.restore();
            // TODO: can we get the state nr?
            state.store(-1);
    		return Double.NEGATIVE_INFINITY;
    	}
    }

	
	

	// BEAST infrastructure stuff
	@Override
    public void setOperatorSchedule(final OperatorSchedule operatorSchedule) {
        operator.setOperatorSchedule(operatorSchedule);
        super.setOperatorSchedule(operatorSchedule);
    }
	
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
	public void optimize(double logAlpha) {
		//operator.optimize(logAlpha);
	}
	
	@Override
	public List<StateNode> listStateNodes() {
		return operator.listStateNodes();
	}

	@Override
	public boolean requiresStateInitialisation() {
		return false;
	}
	
	@Override
	public String toString() {
		String s = operator.toString();
		String label = operator.getName();
		int i = label.indexOf("operator");
		if (i >= 0) {
			label = label.substring(0, i) + "DAOprtr" + label.substring(i + 8);
		} else {
			label = this.getClass().getName();
		}
		s = s.replaceAll(operator.getName(), label);
		return s;
	}
	
} // class DelayedAcceptanceOperator