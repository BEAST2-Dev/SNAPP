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
//  time of full calculation over all samples
//  time of approx calculation
//  rate of rejection rate
//  final acceptance rate


@Description("An operator that uses an approximate likelihood based on genetic distances together with sampling of constant sites ")
public class DistDAOperator extends Operator implements Loggable{
	public Input<Operator> operatorInput = new Input<Operator>("operator", "Operator for proposing moves", Validate.REQUIRED);
	public Input<Distribution> priorInput = new Input<Distribution>("prior", "prior used when likelihood is approximated", Validate.REQUIRED);  
	public Input<State> stateInput = new Input<State>("state", "state object for which we do proposals", Validate.REQUIRED);
	public Input<SnAPTreeLikelihood> treeLikelihoodInput = new Input<SnAPTreeLikelihood>("treelikelihood", "SNAPP tree likelihood for the tree", Validate.REQUIRED);
	public Input<Boolean> skipDAInput = new Input<Boolean>("skipDA","skip delayed acceptance step but output approx likelihood", Validate.REQUIRED);
			
	
	
	//public Input<IntegerParameter> totalCountInput = new Input<IntegerParameter>("totalSiteCount", "Sampled count of the total number of sites", Validate.REQUIRED);

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
	double oldPstar;
	double newPstar;
	double oldPrior;
	double newPrior;  //TODO Make these local variables after debugging
  

	//IntegerParameter ascSiteCount;
	//IntegerParameter totalSiteCountParameter = null;
	//int totalSiteCount;

	int numSampledSites;

	double [][] baseDist; //Distances computed from sampled sites only.
	double [][] baseVar; //Variances computed from sampled sites only.
	double [][] baseW; //Summed weights of sampled sites for which both taxa are present

	State state;   
	TreeInterface tree = null;
	int numProposals = 0;
	int numInitialPass = 0;
	boolean skipDA = false;

	public void initAndValidate() {
		operator = operatorInput.get();
		treelikelihood = treeLikelihoodInput.get();
		siteModel = (SiteModel.Base) treelikelihood.siteModelInput.get();
		substitutionmodel = ((SnapSubstitutionModel)siteModel.substModelInput.get());
		prior = priorInput.get();
		skipDA = skipDAInput.get().booleanValue();

		tree = treelikelihood.treeInput.get();
		state = stateInput.get();
		data = (Data)treelikelihood.dataInput.get();
		numSampledSites = data.getSiteCount();

		//Compute empirical distances from the data, using only those sites in the original data (not sampled)
		calcBaseDistAndCovariance();    


		System.err.println("DistDAOperator is in development and should not be used!");


	}

	public void init(PrintStream out) { 
		out.print("approxLogPost\t");
	}

	public void log(int sample, PrintStream out) {
		out.print(""+ oldPstar  + "\t");
		oldPstar = 0.0;
	}

	public void close(PrintStream out) {
		// nothing to do
	}


	@Override
	public double proposal()  {

		try {
			//int z = numSampledSites;			
			numProposals++;

			//Save old state
			Node oldRoot = tree.getRoot().copy();
			Double [] oldCoalescenceRate = substitutionmodel.m_pCoalescenceRate.get().getValues();
			double oldU = substitutionmodel.m_pU.get().getValue();
			double oldV  = substitutionmodel.m_pV.get().getValue();           
			double oldPrior = prior.getCurrentLogP();
			//Integer[]  oldAscSiteCounts = ascSiteCount.getValues();   
			double oldConstSiteProb = 1.0 - treelikelihood.getNewProbVariableSites();    
			double[][] dThetaOld = computeExpectedDist(oldRoot, oldCoalescenceRate, oldU, oldV);
		

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
            
            oldPstar = calculateLogPstar(1.0-oldConstSiteProb,dThetaOld,oldPrior);
            
            if (skipDA) 
            	return logHastingsRatio;
            
			double[][] dThetaNew = computeExpectedDist(newRoot, newCoalescenceRate, newU, newV);
			double newPrior = prior.getCurrentLogP();
			newPstar = calculateLogPstar(1.0-newConstSiteProb,dThetaNew,newPrior);
		
			
			double logAlpha = Math.min(0.0, logHastingsRatio + newPstar - oldPstar);
			
			if (Randomizer.nextDouble() < Math.exp(logAlpha)) {
				// accept. Now compute the appropriate Hasting's ratio.
				numInitialPass++;


				logHastingsRatio =  oldPstar -  newPstar;
			} else {
				//Reject
				return Double.NEGATIVE_INFINITY;
			}
			
			
			return logHastingsRatio;
		}  catch (Exception e) {
			e.printStackTrace();
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
					int nkx = lineageCounts[taxon1];
					int nky = lineageCounts[taxon2];
					double weight = _data.getPatternWeight(k);
					if (weight > 0 && nkx > 0 && nky > 0) {

						double gk0i = (double)sitePattern[taxon1]/nkx;
						double gk1i = 1.0 - gk0i;
						double gk0j = (double)sitePattern[taxon2]/nky;
						double gk1j = 1.0 - gk0j;
						
						d+=weight*(gk0i*gk1j + gk1i * gk0j);
					
						Kxy += weight;
					}
				}

				//              Here we assume that the sampled ascertained sites are *not* taken into account.
				//        These come into play when we compute the likelihood.
				return d / Kxy;
			}
		};
		
		

		
		//We compute empirical variance estimates 
		Distance v = new Distance.Base() {
			@Override
			public double pairwiseDistance(int taxon1, int taxon2) {
				double Kxy = 0;
				double v = 0;
				double d = 0;
				double d2 = 0;
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
						
						d+=weight*(gk0i*gk1j + gk1i * gk0j);
						d2+=weight*((gk0i*gk1j + gk1i * gk0j)*(gk0i*gk1j + gk1i * gk0j));
					
						Kxy += weight;
					}
				}
				double vsite = d2/Kxy - (d/Kxy)*(d/Kxy);
				v = vsite/Kxy;
                //v = v/Kxy; //Compute per site variance
				return v;     
			}
		};
		

		((Distance.Base)d).setPatterns(data);
		((Distance.Base)v).setPatterns(data);


		int nrOfTaxa = data.getTaxonCount();
		baseDist = new double[nrOfTaxa][nrOfTaxa];
		baseVar = new double[nrOfTaxa][nrOfTaxa];
		

		for (int i = 0; i < nrOfTaxa; i++) {
			for (int j = i; j < nrOfTaxa; j++) {
				baseDist[i][j] = d.pairwiseDistance(i,  j);
				baseDist[j][i] = baseDist[i][j];
				baseVar[i][j] = v.pairwiseDistance(i, j);
				baseVar[j][i] = baseVar[i][j];
			} 
		} 
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

	
	
	private double calculateLogPstar(double p, double[][] mu, double logPrior) {
        int ntax = baseDist.length;
        double logZ = ntax*(ntax+1.0)/4.0 * Math.log(2.0*Math.PI);      
        double logL = 0.0;
        for (int i=0;i<ntax;i++) {
            for (int j=i;j<ntax;j++) {
                double mu_ij = mu[i][j]/p;
                double d_ij = baseDist[i][j];
                double v_ij = baseVar[i][j];
                logL += -0.5 * (mu_ij - d_ij) * (mu_ij - d_ij) /v_ij;
                logZ += 0.5 * Math.log(v_ij);
                //System.out.println(""+((mu[i][j]-d_ij)));
            }
        }   
        
        return -logZ + logL + logPrior;
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