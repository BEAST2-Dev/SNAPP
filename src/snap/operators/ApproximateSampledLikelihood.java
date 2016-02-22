package snap.operators;


import snap.likelihood.SnAPTreeLikelihood;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;

@Description("Approximate Likelihood based on subsample of lineages")
public class ApproximateSampledLikelihood extends BEASTObject implements ApproximateLikelihood {
	public Input<Distribution> priorInput = new Input<Distribution>("prior", "prior used when likelihood is approximated", Validate.REQUIRED);
	public Input<SnAPTreeLikelihood> treeLikelihoodInput = new Input<SnAPTreeLikelihood>("treelikelihood", "SNAPP tree likelihood for the tree", Validate.REQUIRED);
	
	Distribution prior = null;
	SnAPTreeLikelihood treelikelihood;
	
	@Override
	public void initAndValidate() {
		treelikelihood = treeLikelihoodInput.get();
    	prior = priorInput.get();
	}

	@Override
	public double approxPosterior(Node root, Double[] coalescentRate, double u, double v) throws Exception {
		double logP = prior.calculateLogP();
		logP  += treelikelihood.calculateLogP();//approxLikelihood(root, coalescentRate, u, v);
		return logP;
	}

//	@Override
//	public double approxLikelihood(Node root, Double[] coalescentRate, double u, double v) {
//		return treelikelihood.getCurrentLogP();
//	}

	@Override
	public double getPriorValue() {
		return prior.getCurrentLogP();
	}
}
