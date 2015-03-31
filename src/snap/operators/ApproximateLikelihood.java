package snap.operators;

import beast.evolution.tree.Node;

public interface ApproximateLikelihood {

	public double approxPosterior(Node root, Double [] coalescentRate, double u, double v) throws Exception;
	
	//public double approxLikelihood(Node root, Double [] coalescentRate, double u, double v);
	
	public double getPriorValue();
}
