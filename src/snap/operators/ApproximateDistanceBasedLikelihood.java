package snap.operators;

import java.util.ArrayList;
import java.util.List;

import snap.likelihood.SnAPTreeLikelihood;
import snap.likelihood.SnapSubstitutionModel;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

@Description("Approximate Likelihood based on distances between lineages")
public class ApproximateDistanceBasedLikelihood extends BEASTObject implements ApproximateLikelihood {
	public Input<Distribution> priorInput = new Input<Distribution>("prior", "prior used when likelihood is approximated", Validate.REQUIRED);
	public Input<SnAPTreeLikelihood> treeLikelihoodInput = new Input<SnAPTreeLikelihood>("treelikelihood", "SNAPP tree likelihood for the tree", Validate.REQUIRED);

    public boolean useMatLabFormulae = false;
    /** probability that sites are variable for the state of the tree before operating on it **/
    //double probVariableSites = 1.0;

    // empirical distance and variance between taxa
    double [][] distance;
    double [][] var;


	// log of constant in approximate tree likelihood
	double K;

	double priorValue;

	Distribution prior = null;
	// data and tree used in approximate tree likelihood
	Alignment data = null;
	SnAPTreeLikelihood treelikelihood;
	TreeInterface tree = null;
	IntegerParameter ascSiteCount;

	@Override
	public void initAndValidate() {
		treelikelihood = treeLikelihoodInput.get();
		SiteModel.Base siteModel = (SiteModel.Base) treelikelihood.siteModelInput.get();
		SnapSubstitutionModel substitutionmodel = ((SnapSubstitutionModel)siteModel.substModelInput.get());
		
    	prior = priorInput.get();
    	tree = treelikelihood.treeInput.get();
    	data = treelikelihood.dataInput.get();
		ascSiteCount = treelikelihood.ascSiteCountInput.get();

		
        if (prior == null) {
        	throw new RuntimeException("DelayedAcceptanceOperator: could not identify prior in posterior input");
        }
        if (data == null || tree == null) {
        	throw new RuntimeException("DelayedAcceptanceOperator: could not identify data or tree in treelikelihood in posterior input");
        }

        calcDistanceAndVariance();
	}


    private void calcDistanceAndVariance() {
    	// Compute and estimate of the genetic distances from the SNP data, together with an estimate of 
    	//sampling variances.
    	
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
							d += weight * (2.0 * rkx * (nkx - rkx))/(nkx * nkx);
						}
						Kxy += weight;
					}
				}
				
				//TODO NEED TO ADD ascertained site counts to Kxy.
				if (ascSiteCount != null) {
					Kxy += ascSiteCount.getValue(0) + ascSiteCount.getValue(1);
				}
				return d / Kxy;
			}
		};
		
		Distance v = new Distance.Base() {
			@Override
			public double pairwiseDistance(int taxon1, int taxon2) {
				double Kxy = 0;
				double v = 0;
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
						double v1 = (rkx * rky + (nkx - rkx) * (nky - rky))/(nkx * nky);
						double v2 = (rkx * rkx * rky + (nkx - rkx) * (nkx - rkx) * (nky - rky))/(nkx * nkx * nky);
						if (taxon1 != taxon2) {
							double v3 = (rkx * rky * rky + (nkx - rkx) * (nky - rky) * (nky - rky))/(nkx * nky * nky);
							v += weight * ((1 - nkx - nky) * v1 * v1 + (nkx - 1) * v2 + (nky - 1) * v3 + v1) /(nkx * nky);
						} else {
							v += weight * 2.0 * ((nkx - 1)/(nkx * nkx)) * ((3.0 - 2.0 * nkx) * v1 * v1 + 2.0 * (nkx - 2.0) * v2 + v1);
						}
						Kxy += weight;
					}
					
					//TODO NEED TO ADD ascertained site counts to Kxy.	
					if (ascSiteCount != null) {
						Kxy += ascSiteCount.getValue(0) + ascSiteCount.getValue(1);
					}
				}
				v = v / (Kxy * Kxy);
				//if (true) return 0.00005;
				return v;			
				}
		};
		
		((Distance.Base)d).setPatterns(data);
		((Distance.Base)v).setPatterns(data);
		
		int nrOfTaxa = data.getNrTaxa();
		distance = new double[nrOfTaxa][nrOfTaxa];
		var = new double[nrOfTaxa][nrOfTaxa];
		
		for (int i = 0; i < nrOfTaxa; i++) {
			for (int j = i; j < nrOfTaxa; j++) {
				distance[i][j] = d.pairwiseDistance(i,  j);
				distance[j][i] = distance[i][j];
				var[i][j] = v.pairwiseDistance(i, j);
				var[j][i] = var[i][j];
			} 
		}
		
		
		// calculate log of constant of approximate likelihood
		K = -Math.log(2.0 * Math.PI) * nrOfTaxa * (nrOfTaxa + 1.0)/4.0;
		for (int i = 0; i < nrOfTaxa; i++) {
			for (int j = i; j < nrOfTaxa; j++) {
				K -= Math.sqrt(var[i][j]);
			}
		}
	}

    
    @Override
    public double getPriorValue() {
    	return priorValue;
    }
    
	/** calculate approximate posterior **/
    @Override
	public double approxPosterior(Node root, Double [] coalescenceRate, double u, double v) throws Exception {
		priorValue = prior.calculateLogP();
		double logP = priorValue;
		logP  += approxLikelihood(root, coalescenceRate, u, v);
		return logP;
	}


	/**
	 * calculate approximate treelikelihood for tree & data
	 * made public for testing purposes  
	 * */
    //@Override
	public double approxLikelihood(Node root, Double [] coalescenceRate, double u, double v) {
    	// do the real work here
    	double mu[][] = calcMu(root, u, v, coalescenceRate);
    	    	
    	double approxL = 0;
    	int nrOfTaxa = distance.length;
		for (int i = 0; i < nrOfTaxa; i++) {
			for (int j = i+ 1; j < nrOfTaxa; j++) {
				approxL += -0.5 * (distance[i][j] - mu[i][j]) * (distance[i][j] - mu[i][j]) / var[i][j];   
			}
		}

		return K  + approxL;
	}
	
	/** calc approximate distances between taxa **/
    public double[][] calcMu(Node root, double u, double v, Double[] coalescenceRate) {
		// calculate estimates of distance between taxa based on the 
		// tree and other parameters
		
		// 1. calc moment generating function
		double [] M = new double[tree.getNodeCount()];
		calcMomentGeneratingFunction(M, root, u, v, coalescenceRate);
		
		// 2. calc approx distances, store result in u
		double [][] mu = new double[var.length][var.length];
		double probVariableSites = treelikelihood.getProbVariableSites();
		calcApproxDistance(mu, M, root, u, v, probVariableSites);		
		return mu;
	}

	List<Node> calcApproxDistance(double[][] mu, double[] M, final Node node,
			final double u, final double v, final double probVariableSites) {
		int x = node.getNr();
		double t = node.getHeight();
		double pi0 = v/(u+v);
		double pi1 = 1.0 - pi0;
		
		
		if (node.isLeaf()) {
			// nx = nr of lineages for node x, does not matter whether they are missing
			int nx = data.getStateCounts().get(x);
//			if (useMatLabFormulae) {
//				mu[x][x] = 2.0 * pi0 * pi1 * (1.0 - M[x]) / probVariableSites;
//			} else {
				mu[x][x] = 2.0 * pi0 * pi1 * (1.0 - M[x]) * (1.0 - 1.0/nx)  / probVariableSites;
//			}
			List<Node> list = new ArrayList<Node>();
			list.add(node);
			return list;
		} else {
			List<Node> left = calcApproxDistance(mu, M, node.getLeft(), pi0, pi1, probVariableSites);
			List<Node> right = calcApproxDistance(mu, M, node.getRight(), pi0, pi1, probVariableSites);
			for (Node lNode : left) {
				for (Node rNode : right) {
					int i = lNode.getNr();
					int j = rNode.getNr();
					mu[i][j] = 2.0 * pi0 * pi1 * (1.0 - Math.exp(-2.0 * (u+v) * t) * M[x]) / probVariableSites;
					mu[j][i] = mu[i][j];
				}
			}
			left.addAll(right);
			return left;
		}
	}

	/**
	 * calculate moment generating function for each node in the tree based on SNAPP parameters, 
	 * and store results in M
	 */
	public void calcMomentGeneratingFunction(double[] M, Node node, double u, double v, Double [] coalescenceRate) {
		if (node.isRoot()) {
			// x = -2(u+v)
			// Mroot = -1/(.5*theta(3)*x - 1); %MGF at root
			// Mroot = -1/(theta(3)*x/2 - 1); %MGF at root
			// Mroot = -1/(x/lambda - 1); %MGF at root
			// Mroot = 1/(1 - x/lambda); %MGF at root
			// Mroot = 1/(1 + 2(u+v)/lambda); %MGF at root
			M[node.getNr()] = 1.0 / (1.0 + 2.0 * (u + v)/coalescenceRate[node.getNr()]);
		} else {
			// M=(exp(T*(x-(2/theta(1))))-1)/(.5*theta(1)*x-1) + M[parent]*exp(T*(x-(2/theta(1))));
			// M=(exp(T*(x-lambda(1)))-1)/(x/lambda(1)-1) + M[parent]*exp(T*(x-lambda(1)));
			// M=(exp(T*(x-lambda))-1)/(x/lambda-1) + exp((x-lambda)* T) * M[parent]
			// M=(exp(tx*(-pi-lambda))-1)/(-pi/lambda-1) + exp((-pi-lambda)* tx) * M[parent]

			double lambda = coalescenceRate[node.getNr()];
			int x = node.getNr();
			double tx = node.getLength();
			double pi = 2.0 * (u + v);
			int parent = node.getParent().getNr();
			if (useMatLabFormulae) {
				M[x]= (Math.exp(tx*(-pi-lambda))-1)/(-pi/lambda-1) + Math.exp((-pi-lambda)* tx) * M[parent];
			} else {		
				//M[x] = (1.0 - Math.exp(-lambda * tx)) * lambda * 
				//		(Math.exp(lambda*tx) -Math.exp(-pi*tx)) / ((lambda + pi) * (Math.exp(lambda * tx) - 1.0));

				M[x] = lambda * (1.0 - Math.exp(-(lambda + pi) * tx)) / (lambda + pi);

				M[x] += Math.exp(-(lambda + pi) * tx) * M[node.getParent().getNr()];
			}	
		}
		if (!node.isLeaf()) {
			calcMomentGeneratingFunction(M, node.getLeft(), u, v, coalescenceRate);
			calcMomentGeneratingFunction(M, node.getRight(), u, v, coalescenceRate);
		}
	}

}