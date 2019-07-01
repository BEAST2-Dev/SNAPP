package snap.likelihood;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import snap.distribution.ChiSquareNoncentralDist;

@Description("Standard rate prior for SnAP analysis, which is a " +
        "gamma distribution over the theta values " +
        "(with parameters alpha and beta). " +
        "Thetas are represented by the coalescenceRate parameter where values are theta=2/coalescenceRate")
public class RatePrior extends Distribution {
    public Input<RealParameter> m_pAlpha = new Input<RealParameter>("alpha", "Alpha parameter for the gamma prior on population size (theta) values", Validate.REQUIRED);
    public Input<RealParameter> m_pBeta = new Input<RealParameter>("beta", "Beta parameter for the gamma prior on population size (theta) values", Validate.REQUIRED);
    public Input<RealParameter> m_pKappa = new Input<RealParameter>("kappa", "prior parameter -- see docs for details");
    public Input<RealParameter> m_pCoalescenceRate = new Input<RealParameter>("coalescenceRate", "Populations sizes for the nodes in the tree", Validate.REQUIRED);
    public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations"); //, Validate.REQUIRED);

    
   enum Priors {
        gamma,inverseGamma,CIR,uniform
    }

   public Input<Priors> m_priors = new Input<Priors>("rateprior", "prior on rates. " +
           "This can be " + Arrays.toString(Priors.values()) + " (default 'gamma')", Priors.gamma, Priors.values());

	int PRIORCHOICE = 0;

    @Override
    public void initAndValidate() {
    	if (m_pKappa.get() == null) {
    		System.err.println("WARNING: kappa parameter not set for SnAPPrior. using default value of 1.0");
    		m_pKappa.setValue(new RealParameter(new Double[]{1.0}), this);
    	}

    	// determine rate prior
    	switch (m_priors.get()) {
	    	case gamma: PRIORCHOICE =0;break;
	    	case inverseGamma: PRIORCHOICE =1;break;
	    	case CIR: PRIORCHOICE =2;break;
	    	case uniform: PRIORCHOICE =3;break;
    	}
    	System.err.println("Rate prior = " + Priors.values()[PRIORCHOICE] + "");
    }


    @Override
    public double calculateLogP() {
        logP = 0.0;

        double alpha = m_pAlpha.get().getValue();
        double beta = m_pBeta.get().getValue();
        double kappa = m_pKappa.get().getValue();
        
        if (outsideBounds(m_pAlpha.get()) || outsideBounds(m_pBeta.get()) || outsideBounds(m_pKappa.get()) || outsideBounds(m_pCoalescenceRate.get())) {
        	logP = Double.NEGATIVE_INFINITY;
        	return logP;
        }
        
        Tree tree = m_pTree.get();
        
        //Gamma values in tree
        RealParameter coalescenceRate = m_pCoalescenceRate.get();

		
		
		if (PRIORCHOICE == 0) {
			//Assume independent gamma distributions for thetas.
			
			//We assume that 2/r has a gamma(alpha,beta) distribution. That means that r has density proportional to
			// 1/(r^2)  * GAMMA(2/r|alpha,beta)
			//which has log (alpha - 1.0)*Math.log(2.0/r) - (beta *(2.0/ r)) - 2*log(r), which in turn simplifies to the expr. below (w/ consts)
		
			for (int iNode = 0; iNode < coalescenceRate.getDimension(); iNode++) {
				double r = coalescenceRate.getValue(iNode);
				logP += -(alpha + 1.0)*Math.log(r) - 2.0* beta / r;
			}
		} else if (PRIORCHOICE == 1) {
			//Assume independent inverse gamma distributions for thetas
			
			//We assume that (2/r) has an inverse gamma (alpha,beta) distribution. That means that r has density proportional to
			// 1/(r^2) * INVGAMMA(2/r|alpha,beta)
			//which has logarithm:
			// -(alpha+1) log(2/r) - beta*r/2 - 2 log(r) = C + (alpha-1) log(r) - beta*r / 2.
			for (int iNode = 0; iNode < coalescenceRate.getDimension(); iNode++) {
				double r = coalescenceRate.getValue(iNode);
				logP += (alpha - 1.0)*Math.log(r) - 0.5* beta * r;
			}
		} else if (PRIORCHOICE == 2) {
	        //> let x be the root.
	        //> r = rate for node x.
    		
			//> let x be the root.
	        //> r = rate for node x.
    		double r = coalescenceRate.getArrayValue(tree.getRoot().getNr());
	        //TEMPORARY!!!! logP += -(alpha + 1.0)*Math.log(r) - 2.0* beta / r;
			
			
			logP += -(alpha*2 + 1.0)*Math.log(r) - 2.0* beta / r;
			
			
			
			/*
			 
			 The CIR process (Cox, Ingersoll and Ross 1985. Econometrica) has SDE
			 dr = \kappa (\theta - r) dt + \sigma \sqrt{r} dz_1
			 has a stationary distribution that is gamma with 
				alpha = 2 \kappa \theta / \sigma^2
			 and 
				beta = 2 \kappa / \sigma^2
			 The correlation between time 0 and time t is  exp(-kappa t).
			 
			 
			 Let c = 2 \kappa / (sigma^2 * (1 - exp(-kappa t)))
			 
			 If we condition on rate r0 at time 0, the distribution of 2*c*rt is non-central
			 chi squared with 
				2q+2 = 4\kappa \theta / sigma^2  
			 degrees of freedom and parameter of non-centrality
				2u = c r_0 exp(-kappa t)
			 
			 
			 Converting these into our set of parameters (alpha, beta, kappa) we have
			 \theta = \alpha / \beta
			 \sigma^2 = 2 \kappa / \beta
			 so
			 c = \beta/(1 - exp(-kappa t))
			 
			 df = 2q+2 = 2 \alpha
			 
			 nc= 2u = 2 \beta r_0 \frac{exp(-kappa t)}{1-\exp(-kappa t)}
			 
			 We are applying the CIR process to THETA = 2/r
			 
			*/
			
	        
			/* Priors for the hyperparameters */
			logP -= Math.log(kappa);
			
				
			
			
			
	        Node [] nodes = tree.getNodesAsArray();
	        
	        for (int iNode = 0; iNode < tree.getNodeCount(); iNode++) {
	        	Node node = nodes[iNode];
	        	if (!node.isRoot()) {
	        		Node parent = node.getParent();
	        		double t = parent.getHeight() - node.getHeight();
	        		r = coalescenceRate.getArrayValue(node.getNr());
					
					
					double r0 = coalescenceRate.getArrayValue(parent.getNr());
										
	        		//double r0 = parent.coalescenceRate(); <=== I think that this might be buggy.
					
					double theta0 = 2.0/r0;			
					
					
					double ekt = Math.exp(-kappa*t);
					double c = beta / (1.0 - ekt);
					double df = 2*alpha;
					double nc = 2.0 * beta * theta0 * ekt / (1.0 - ekt);
					double theta = 2.0/r;
					
					//System.err.println("kappa,t,beta,ekt,df,nc = "+kappa+", "+t+", "+beta+", "+ekt+", "+df+", "+nc);
					
					double p = ChiSquareNoncentralDist.density(df, nc, 2*c*theta);
					
	        		logP += Math.log(p) - 2.0 * Math.log(r);
	        	}
	        }
			
			
       		if (Double.isInfinite(logP)) {
       			// take care of numeric instability
       			logP = Double.NEGATIVE_INFINITY;
       		}
	        
			
		} else {
			//Assume that rate has uniform distribution on [[0,1000]
			for (int iNode = 0; iNode < coalescenceRate.getDimension(); iNode++) {
				double r = coalescenceRate.getValue(iNode);
				if (r>1000000000.0 || r<0.0)
					return Double.NEGATIVE_INFINITY;
			}
		}
		
		
		
		
		
		
        return logP;
    } // calculateLogLikelihood

    private boolean outsideBounds(RealParameter realParameter) {
    	if (realParameter == null) {
    		return false;
    	}
    	for (int i = 0; i < realParameter.getDimension(); i++) {
    		double d = realParameter.getArrayValue(i);
    		if (d < realParameter.getLower() || d > realParameter.getUpper()) {
    			return true;
    		}
    	}
		return false;
	}

    
	@Override
	public List<String> getArguments() {return null;}

	@Override
	public List<String> getConditions() {return null;}

	@Override
	public void sample(State state, Random random) {}

}
