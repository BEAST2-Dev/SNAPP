package snap.operators;

import beast.evolution.tree.Node;

public interface ApproximateLikelihood {

	public void initialise();
	public double evaluate(Node root, Double [] coalescenceRate, double u, double v) throws Exception;
	public double approxLikelihood(Node root, Double [] coalescenceRate, double u, double v);
	public double getPriorValue();
	public void setProbVariableSites(double probVariableSites);
}
