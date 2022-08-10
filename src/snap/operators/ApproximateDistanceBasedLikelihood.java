package snap.operators;

import java.util.ArrayList;
import java.util.List;

import snap.likelihood.SnAPTreeLikelihood;
import snap.likelihood.SnapSubstitutionModel;
import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.distance.Distance;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;

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
        if (!treelikelihood.m_usenNonPolymorphic.get() && ascSiteCount==null) {
        	throw new RuntimeException("DelayedAcceptanceOperator: distance based approximation not available unless nonPolymorphic sites included or estimated");
        }
        if (treelikelihood.hasDominantMarkers.get()) {
        	throw new RuntimeException("DelayedAcceptanceOperator: distance-based delayed acceptance not available for dominant markers (if you want it, contact David Bryant)");
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
				
//              Here we assume that the sampled ascertained sites are *not* taken into account.
//				These come into play when we compute the likelihood.
//				if (ascSiteCount != null) {
//					Kxy += ascSiteCount.getValue(0) + ascSiteCount.getValue(1);
//				}
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
					
//					Here we assume that the sampled ascertained sites are *not* taken into account. 
//					if (ascSiteCount != null) {
//						Kxy += ascSiteCount.getValue(0) + ascSiteCount.getValue(1);
//					}
				}
				v = v / (Kxy * Kxy);
				//if (true) return 0.00005;
				return v;			
				}
		};
		
		((Distance.Base)d).setPatterns(data);
		((Distance.Base)v).setPatterns(data);
		
		int nrOfTaxa = data.getTaxonCount();
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
				//TODO: Account for ascertained sites here.
				approxL += -0.5 * (distance[i][j] - mu[i][j]) * (distance[i][j] - mu[i][j]) / var[i][j];   
			}
		}

		return K  + approxL;
	}
	
	/** 
	 * calc expected genetic distances between species, given population tree and parameters 
	 * **/
    public double[][] calcMu(Node root, double u, double v, Double[] coalescenceRate) {
		// calculate estimates of distance between taxa based on the 
		// tree and other parameters
		
		// 1. eval moment generating functions
		double [] M = new double[tree.getNodeCount()];
		calcMomentGeneratingFunction(M, root, -2*(u+v), coalescenceRate);
		
		// 2. calc approx distances, store result in u
		double [][] mu = new double[var.length][var.length];
		//double probVariableSites = treelikelihood.getProbVariableSites();
		calcApproxDistance(mu, M, root, u, v);	
		return mu;
	}

	List<Node> calcApproxDistance(double[][] mu, double[] M, final Node node,
			final double u, final double v) {
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
				mu[x][x] = 2.0 * pi0 * pi1 * (1.0 - M[x]) * (1.0 - 1.0/nx);  //Multiplier is prob of not getting same individual.
//			}
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

	/**
	 * calculate moment generating function for each node in the tree based on SNAPP parameters, 
	 * evaluate the mgf at x and store the result in M.
	 */
	public void calcMomentGeneratingFunction(double[] M, Node node, double x, Double [] coalescenceRate) {
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

}
