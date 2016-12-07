package snap.operators;



import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import snap.Data;
import snap.distribution.GammaDist;
import snap.likelihood.SnAPTreeLikelihood;
import snap.likelihood.SnapSubstitutionModel;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.OperatorSchedule;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.Operator;
import beast.core.StateNode;
import beast.core.Loggable;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

//import java.util.Arrays;
//
//import beast.core.Citation;
//import beast.core.Description;
//import beast.core.Input;
//import beast.core.Input.Validate;
//import beast.core.parameter.RealParameter;
//import beast.evolution.tree.TreeInterface;
import static org.apache.commons.math.special.Gamma.logGamma;


//Things to measure:
//	time of full calculation over all samples
//  time of approx calculation
//  rate of rejection rate
//  final acceptance rate


@Description("An operator that uses an approximate likelihood based on genetic distances together with sampling of constant sites ")
public class DistDAOperator extends Operator implements Loggable{
	public Input<Operator> operatorInput = new Input<Operator>("operator", "Operator for proposing moves", Validate.REQUIRED);
	public Input<Distribution> priorInput = new Input<Distribution>("prior", "prior used when likelihood is approximated", Validate.REQUIRED);	
	public Input<State> stateInput = new Input<State>("state", "state object for which we do proposals", Validate.REQUIRED);
    public Input<SnAPTreeLikelihood> treeLikelihoodInput = new Input<SnAPTreeLikelihood>("treelikelihood", "SNAPP tree likelihood for the tree", Validate.REQUIRED);
	public Input<IntegerParameter> totalCountInput = new Input<IntegerParameter>("totalSiteCount", "Sampled count of the total number of sites", Validate.REQUIRED);

	Operator operator; //Operator for proposing trees + parameters.
	
	// TODO: take site model in account in approximation?
	// TODO: take clock model in account in approximation?
	
	double priorValue;
	Distribution prior = null;
	
	SiteModel.Base siteModel;
	SnapSubstitutionModel substitutionmodel;
	SnAPTreeLikelihood treelikelihood;
	ApproximateLikelihood approxLiklihood;
	
	Alignment data = null;
    

	//IntegerParameter ascSiteCount;
	int totalSiteCount;
	
	int numSampledSites;
	boolean resampleConstant = true; //First time through, X needs to be resampled from P(X|theta)
	
	double [][] baseDist; //Distances computed from sampled sites only.
	double [][] baseVar; //Variances computed from sampled sites only.
	double [][] baseW; //Summed weights of sampled sites for which both taxa are present
	
    State state;   
	TreeInterface tree = null;
	int numProposals = 0;
	int numInitialPass = 0;
	
	public void initAndValidate() {
		operator = operatorInput.get();
		treelikelihood = treeLikelihoodInput.get();
    	siteModel = (SiteModel.Base) treelikelihood.siteModelInput.get();
    	substitutionmodel = ((SnapSubstitutionModel)siteModel.substModelInput.get());
    	prior = priorInput.get();
    	
    	tree = treelikelihood.treeInput.get();
    	state = stateInput.get();
    	data = (Data)treelikelihood.dataInput.get();
    	numSampledSites = data.getSiteCount();
   		//ascSiteCount = treelikelihood.ascSiteCountInput.get();
        totalSiteCount = treelikelihood.totalSiteCountInput.get().getValue();
        
   		//Compute empirical distances from the data, using only those sites in the original data (not sampled)
		calcBaseDistAndCovariance();		
		
		
		System.err.println("DistDAOperator is in development and should not be used!");

		
    }

	public void init(PrintStream out) { 
		out.print("acceptance rates: "+operator.getID()+"\t");
	}

    public void log(int sample, PrintStream out) {
    	out.print(""+((double)numInitialPass/numProposals) + "\t");
	}
    
    public void close(PrintStream out) {
	// nothing to do
	}


	@Override
    public double proposal()  {
		
    	try {
    		int z = numSampledSites;
    				
    		numProposals++;
    		
    		//Save old state
    		Node oldRoot = tree.getRoot().copy();
        	Double [] oldCoalescenceRate = substitutionmodel.m_pCoalescenceRate.get().getValues();
        	double oldU = substitutionmodel.m_pU.get().getValue();
        	double oldV  = substitutionmodel.m_pV.get().getValue();        	 	
        	double oldPrior = prior.getCurrentLogP();
	    	//Integer[]  oldAscSiteCounts = ascSiteCount.getValues();	  
	    	int oldTotalSiteCount = totalSiteCount;
	    	double oldConstSiteProb = 1.0 - treelikelihood.getNewProbVariableSites();	   
	    	double[][] dThetaOld = computeExpectedDist(oldRoot, oldCoalescenceRate, oldU, oldV);

	    	if (resampleConstant) {
	    		//First time through the algorithm we need to sample the initial number of constant sites.
	            oldTotalSiteCount = sampleNegBinomial(numSampledSites,oldConstSiteProb); 
	            resampleConstant = false;
	       	}
	    	
//	    	
//	    	double[] oldConstSiteProbs = new double[2];
//	    	oldConstSiteProbs[0]  = treelikelihood.getSitesProbs(-2);
//	    	oldConstSiteProbs[1] = treelikelihood.getSitesProbs(-1);
	    		    	
	    	//Propose new state
	    	double logHastingsRatio = operator.proposal();
	    	
	    	
	    	//Extract parameters from the new state
  			Node newRoot = tree.getRoot();
	    	Double [] newCoalescenceRate = substitutionmodel.m_pCoalescenceRate.get().getValues();
	    	double newU = substitutionmodel.m_pU.get().getValue();
	    	double newV = substitutionmodel.m_pV.get().getValue();
	    	state.storeCalculationNodes();
            state.checkCalculationNodesDirtiness();
         
            if (logHastingsRatio == Double.NEGATIVE_INFINITY) {
	    		// instant reject
	    		// need to store state, so it is restored properly
	    		return Double.NEGATIVE_INFINITY;
	    	}

            
            // Propose a new count of constant sites.
	    	double newConstSiteProb = 1.0 - treelikelihood.getNewProbVariableSites();	    			
	    	int newTotalSiteCount = sampleNegBinomial(numSampledSites,newConstSiteProb); 
	    	

	    	//private double calculateLogPhi(int z, int N, int Nstar, double[][] mu) {	
		
	    	double oldPhi = calculateLogPhi(z,oldTotalSiteCount,oldTotalSiteCount,dThetaOld);
	    	double[][] dThetaNew = computeExpectedDist(newRoot, newCoalescenceRate, newU, newV);
	    	double newPhi = calculateLogPhi(z,newTotalSiteCount,oldTotalSiteCount,dThetaNew);
	    	
	    	double newPrior = prior.getCurrentLogP();
	    	double oldLogNegBin = logNegBinomialDensity(oldTotalSiteCount,numSampledSites,1.0-oldConstSiteProb);
	    	double newLogNegBin = logNegBinomialDensity(newTotalSiteCount,numSampledSites,1.0-newConstSiteProb);
	    	
	    	double logAlpha = logHastingsRatio -Math.log(newTotalSiteCount) + Math.log(oldTotalSiteCount);
	    	logAlpha += newPhi - oldPhi + newPrior - oldPrior + oldLogNegBin - newLogNegBin;	    	
	    	logAlpha = Math.min(logAlpha,0.0);
	    		    	
            
            if (Randomizer.nextDouble() < Math.exp(logAlpha)) {
	        	// accept. Now compute the appropriate Hasting's ratio.
            	numInitialPass++;
            	double oldPhiBack = calculateLogPhi(z,oldTotalSiteCount,newTotalSiteCount,dThetaOld);
            	double newPhiBack = calculateLogPhi(z,newTotalSiteCount,newTotalSiteCount,dThetaNew);
            	double logAlphaBack = -logHastingsRatio +Math.log(newTotalSiteCount) - Math.log(oldTotalSiteCount);
            	logAlphaBack += oldPhiBack - newPhiBack + oldPrior - newPrior + newLogNegBin - oldLogNegBin;
            	logAlphaBack = Math.min(logAlphaBack,0.0);
            	
            	logHastingsRatio += logAlphaBack - logAlpha;
            } else {
            	//Reject
            	return Double.NEGATIVE_INFINITY;
            }
            return logHastingsRatio;
    	}  catch (Exception e) {
    		e.printStackTrace();
    		System.err.println("DelayedAcceptanceOperator.proposal " + e.getMessage());
            state.restore();
            // TODO: can we get the state nr?
            state.store(-1);
    		return Double.NEGATIVE_INFINITY;
    	}
    }

	//TODO: Move these to a utility class
	private void calcBaseDistAndCovariance() {
		
		// Computes empirical distances from the SNP data, together with an estimate of 
    	//sampling variances and the total summed weight of sites for which there is a state for both taxa.
		//Note: these are computed from the input (sampled) data set only, and ignore any ascertainment 
		//process that went on.
		//Note: formula differ slightly from Gordon's thesis since with a single population we assume sampling
		//of individuals without replacement (as this is what the coalescent formula provides)
    	
		Distance d = new Distance.Base() {
			@Override
			public double pairwiseDistance(int taxon1, int taxon2) {
				double Kxy = 0;
				double d = 0;
				snap.Data _data = (snap.Data) data;
				for (int k = 0; k < _data.getPatternCount(); k++) {
					int [] lineageCounts = _data.getPatternLineagCounts(k);
					int [] sitePattern = _data.getPattern(k);
					double rkx = sitePattern[taxon1];
					double rky = sitePattern[taxon2];
					int nkx = lineageCounts[taxon1];
					int nky = lineageCounts[taxon2];
					double weight = _data.getPatternWeight(k);
					if (weight > 0 && nkx > 0 && nky > 0) {
						if (taxon1 != taxon2) {
							d += weight * (rkx * (nky - rky) + rky * (nkx - rkx))/(nkx * nky);
						} else {
							d += weight * (2.0 * rkx * (nkx - rkx))/(nkx * (nkx-1));
						}
						Kxy += weight;
					}
				}
				
//              Here we assume that the sampled ascertained sites are *not* taken into account.
//				These come into play when we compute the likelihood.
				return d / Kxy;
			}
		};
		
		//We compute variance estimates based on pg 384 of Nei and Roychoudhury (1973).
		//Note:- these formula are based on the assumption that *all* sites are present. We need to account for
		//ascertained sites separately.
		//TODO: Check these for sampling individuals without replacement in a single species/population.
		Distance v = new Distance.Base() {
			@Override
			public double pairwiseDistance(int taxon1, int taxon2) {
				double Kxy = 0;
				double v = 0;
				snap.Data _data = (snap.Data) data;
				for (int k = 0; k < _data.getPatternCount(); k++) {
					int [] lineageCounts = _data.getPatternLineagCounts(k);
					int [] sitePattern = _data.getPattern(k);
					int nki = lineageCounts[taxon1];
					int nkj = lineageCounts[taxon2];
					double weight = _data.getPatternWeight(k);
					if (weight > 0 && nki > 0 && nkj > 0) {
						
						double gk0i = (double)sitePattern[taxon1]/nki;
						double gk1i = 1.0 - gk0i;
						double gk0j = (double)sitePattern[taxon2]/nkj;
						double gk1j = 1.0 - gk0j;
						double vk;
						
						if (taxon1 != taxon2) {
							vk = 1.0/(nki*nkj) * ((1-nki-nkj)*(gk0i*gk0j + gk1i*gk1j)*(gk0i*gk0j + gk1i*gk1j) 
									+ (nki - 1)*(gk0i*gk0i*gk0j + gk1i*gk1i*gk1j) + (nkj-1)*(gk0j*gk0j*gk0i + gk1j*gk1j*gk1i) 
									+ gk0i*gk0j + gk1i*gk1j );						
						} else {
							vk = 2*(nki - 1.0)/(nki*nki*nki) * ((3.0 - 2.0*nki)*(gk0i*gk0i + gk1i*gk1i)*(gk0i*gk0i + gk1i*gk1i) 
									+ gk0i*gk0i + gk1i*gk1i + 2.0*(nki-2.0)*(gk0i*gk0i*gk0i + gk1i*gk1i*gk1i));
						}
						v += weight * vk;
						Kxy += weight;
					}
				}
				v = v / (Kxy * Kxy);
				return v;			
				}
		};
		
		Distance w = new Distance.Base() {		
			/* Returns summed weight of all sites for which both pop1 and pop2 have at least one sample.
			 * 
			 */
			@Override
			public double pairwiseDistance(int taxon1, int taxon2) {
				double Kxy = 0;
				snap.Data _data = (snap.Data) data;
				for (int k = 0; k < _data.getPatternCount(); k++) {
					int [] lineageCounts = _data.getPatternLineagCounts(k);
					int nkx = lineageCounts[taxon1];
					int nky = lineageCounts[taxon2];
					double weight = _data.getPatternWeight(k);
					if (weight > 0 && nkx > 0 && nky > 0) {
						Kxy += weight;
					}
				}
				return Kxy;				
			}
		};
		
		
		
		((Distance.Base)d).setPatterns(data);
		((Distance.Base)v).setPatterns(data);
		((Distance.Base)w).setPatterns(data);

		
		int nrOfTaxa = data.getTaxonCount();
		baseDist = new double[nrOfTaxa][nrOfTaxa];
		baseVar = new double[nrOfTaxa][nrOfTaxa];
		baseW = new double[nrOfTaxa][nrOfTaxa];
		
		for (int i = 0; i < nrOfTaxa; i++) {
			for (int j = i; j < nrOfTaxa; j++) {
				baseDist[i][j] = d.pairwiseDistance(i,  j);
				baseDist[j][i] = baseDist[i][j];
				baseVar[i][j] = v.pairwiseDistance(i, j);
				baseVar[j][i] = baseVar[i][j];
				baseW[i][j] = w.pairwiseDistance(i, j);
				baseW[j][i] = baseW[i][j];
			} 
		}	
	}
	

	private int sampleNegBinomial(int n,double p) {
		/* Samples from a negative binomial 
		 * P[X = N|n,p] = \binom{N-1}{z-1} p^z (1-p)^(N-z)
		 * 
		 * Use the fact that X ~ Poisson(lambda) where lambda is Gamma with 
		 * shape = n and rate \theta = (1-p)/p. 
		 * TODO Check these formulas.
		 */
		double m = GammaDist.inverseF(n, (1.0-p)/p, 6, Randomizer.nextDouble());
		int x =  (int) Randomizer.nextPoisson(m); 
		return x+n;
	}
	
	//TODO: Move these into a 'DistanceBasedUtilities' class.
	private double[][] computeExpectedDist(Node root, Double [] coalescenceRate, double u, double v) {
		// calculate expected differences between individuals based on the 
		// tree and other parameters.
		//
		//Note - these differ slightly from formula in Gordon's thesis since we don't correct for 

		// 1. eval moment generating functions
		double [] M = new double[tree.getNodeCount()];
		calcMomentGeneratingFunction(M, root, -2*(u+v), coalescenceRate);

		// 2. calc approx distances, store result in u
		int ntax = data.getTaxonCount();
		double [][] mu = new double[ntax][ntax];
		//double probVariableSites = treelikelihood.getProbVariableSites();
		calcApproxDistance(mu, M, root, u, v);	
		return mu;		
	}

	
	/**
	 * calculate moment generating function for each node in the tree based on SNAPP parameters, 
	 * evaluate the mgf at x and store the result in M.
	 */
	private void calcMomentGeneratingFunction(double[] M, Node node, double x, Double [] coalescenceRate) {
		if (node.isRoot()) {
			double lambda = coalescenceRate[node.getNr()];
			M[node.getNr()] = lambda / (lambda - x);
		} else {
			double lambda_i = coalescenceRate[node.getNr()];
			int i = node.getNr();
			double ti = node.getLength();
			int parent = node.getParent().getNr();			
			M[i] = lambda_i/(lambda_i-x) * (1 - Math.exp((x - lambda_i)*ti)) + Math.exp((x - lambda_i)*ti)*M[parent];			
		}
		if (!node.isLeaf()) {
			calcMomentGeneratingFunction(M, node.getLeft(), x, coalescenceRate);
			calcMomentGeneratingFunction(M, node.getRight(), x, coalescenceRate);
		}
	}
	
	private List<Node> calcApproxDistance(double[][] mu, double[] M, final Node node,
			final double u, final double v) {
		int x = node.getNr();
		double t = node.getHeight();
		double pi0 = v/(u+v);
		double pi1 = 1.0 - pi0;
		
		
		if (node.isLeaf()) {
			mu[x][x] = 2.0 * pi0 * pi1 * (1.0 - M[x]);  
			List<Node> list = new ArrayList<Node>();
			list.add(node);
			return list;
		} else {
			List<Node> left = calcApproxDistance(mu, M, node.getLeft(), pi0, pi1);
			List<Node> right = calcApproxDistance(mu, M, node.getRight(), pi0, pi1);
			//Note: implicit assumption that the population tree is ultrametric.
			double mu_ij = 2.0 * pi0 * pi1 * (1.0 - Math.exp(-2.0 * (u+v) * t) * M[x]);			
			for (Node lNode : left) {
				for (Node rNode : right) {
					int i = lNode.getNr();
					int j = rNode.getNr();
					mu[i][j] = mu_ij;
					mu[j][i] = mu_ij;
				}
			}
			left.addAll(right);
			return left;
		}
	}
	
	private double calculateLogPhi(int z, int N, int Nstar, double[][] mu) {	
		double logL = 0.0;
		
		//TODO: Need to deal with weighting better in the case of missing sites.
		
		int ntax = baseDist.length;
		for (int i=0;i<ntax;i++) {
			for (int j=i;j<ntax;j++) {
				double mu_ij = mu[i][j];
				double d_ij = ((double)z/N)*baseDist[i][j];
				double v_ij = ((double)z/Nstar)*((double)z/Nstar) * baseVar[i][j];
				logL += -0.5 * (mu_ij - d_ij) * (mu_ij - d_ij) /v_ij;
			}
		}
		return logL;		
	}

	private double logNegBinomialDensity(int N,int z,double p) {
		/* Computes log( \binom{N-1}{z-1} p^z (1-p)^(N-z) )
		 *
		 */
		double logp = logGamma(N) - logGamma(z) - logGamma(N-z+1);
		logp += z*Math.log(p) + (N-z)*Math.log(1.0-p);				
		return logp;
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