package snap.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;


@Description("Operator for merging and splitting branches of a species tree " +
		"from Yang and Ranala PNAS 2010.")
public class MergeSplitSpeciesTree extends TreeOperator {
	public Input<RealParameter> coalescentRateInput = new Input<RealParameter>("coalescenceRate", "population sizes", Validate.REQUIRED);
	public Input<Double> scaleFactorInput = new Input<Double>("scalefactor", "scale factor which determines how bold proposals for theta are when splitting", 1.0);
	
	Tree tree;
	RealParameter coalescentRate;
	double epsilon;
	
	@Override
	public void initAndValidate() throws Exception {
		tree = treeInput.get();
		coalescentRate = coalescentRateInput.get();
		epsilon = scaleFactorInput.get();
	}
	
	@Override
	public double proposal() {
		if (Randomizer.nextBoolean()) {
			return split();
		} else {
			return merge();
		}
	}

	/** return list of nodes that can be merged 
	 * those are the ones at the bottom of the tree 
	 * with branch length = 0
	 **/
	List<Node> getmergecandidates() {
		Node [] nodes = tree.getNodesAsArray();
		List<Node> mergecandidates = new ArrayList<Node>();
		for (Node node : nodes) {
			Node left = node.getLeft();
			Node right = node .getRight();
			// candidates are those nodes that 
			// - are not merged already
			// - have two children that are merged OR are leafs 
			if (!node.isRoot() && 
					(left != null && left.getLength() > 0 && 
						(left.isLeaf() || (isMerged(left.getLeft()) && isMerged(left.getRight())))) &&
					(right != null && right.getLength() > 0&& 
						(right.isLeaf() || (isMerged(right.getLeft()) && isMerged(right.getRight()))))) {
				mergecandidates.add(node);
			}
		}
		return mergecandidates;
	} // getmergecandidates

	/** return list of nodes that can be split 
	 * these are nodes that are merged already 
	 * and have a parent that is not merged
	 **/
	List<Node> getsplitcandidates() {
		Node [] nodes = tree.getNodesAsArray();
		List<Node> splitcandidates = new ArrayList<Node>();
		for (Node node : nodes) {
			Node left = node.getLeft();
			Node right = node .getRight();
			// candidates are those nodes that 
			// - are merged already
			// - have a parent that is not merged 
			if (!node.isRoot() && 
					(left != null && isMerged(left)) &&
					(right != null && isMerged(right))) {
				splitcandidates.add(node);
			}
		}
		return splitcandidates;
	} // getsplitcandidates

	
	/** merge move from Yang and Ranala PNAS 2010 **/
	private double merge() {
		tree.startEditing(this);
		List<Node> mergecandidates = getmergecandidates();
		
		// check there are candidates that can be merged
		if (mergecandidates.size() == 0) {
			return Double.NEGATIVE_INFINITY;
		}
		
		// Randomly pick one of the candidates
		Node node = mergecandidates.get(Randomizer.nextInt(mergecandidates.size()));
		
		// merge with nodes below
		double u1 = node.getHeight();
		node.setHeight(node.getLeft().getHeight());
		
		// calc Hastings ratio
		List<Node> splitcandidates = getsplitcandidates();
		
		int x = mergecandidates.size(); // number of nodes feasible for joining in  S
		int y = splitcandidates.size(); // number of nodes feasible for splitting in S*
		double thetaj = coalescentRate.getArrayValue(node.getLeft().getNr());
		double thetak = coalescentRate.getArrayValue(node.getRight().getNr());
		double tau = node.getParent().getHeight();
		double logHR = Math.log(x)- Math.log(y) + 2.0*Math.log(3.0*u1)-3.0*Math.log(tau) 
				- Math.log(epsilon * thetaj) - Math.log(epsilon * thetak);
		return logHR;
	} // merge

	private boolean isMerged(Node node) {
		// check that length == 0
		// take in account numerical instability
		return node.getLength() < 1e-100;
	}

	/** split move from Yang and Ranala PNAS 2010 **/
	private double split() {
		tree.startEditing(this);		
		List<Node> splitcandidates = getsplitcandidates();
		
		// check there are candidates that can be split
		if (splitcandidates.size() == 0) {
			return Double.NEGATIVE_INFINITY;
		}
		
		// Randomly pick one of the candidates
		Node node = splitcandidates.get(Randomizer.nextInt(splitcandidates.size()));
		
		// split by chaning height of node
		double upperLimit = node.getParent().getHeight();
		double height = getParabolicDistributionSample(upperLimit); 
		node.setHeight(height);

		Node left = node.getLeft();
		Node right = node.getRight();
		double thetai = coalescentRate.getValue(node.getNr());
		double thetaj = thetai * Math.exp(epsilon * (Randomizer.nextDouble() - 0.5));
		double thetak = thetai * Math.exp(epsilon * (Randomizer.nextDouble() - 0.5));
		coalescentRate.setValue(left.getNr(), thetaj);
		coalescentRate.setValue(right.getNr(), thetak);
		
		List<Node> mergecandidates = getmergecandidates();
		int x = splitcandidates.size(); // number of nodes feasible for splitting in  S
		int y = mergecandidates.size(); // number of nodes feasible for joining in S*
		double logHR = Math.log(x) - Math.log(y) + 3 * Math.log(upperLimit) - Math.log(2 * height * height) 
				+ Math.log(epsilon * thetaj) + Math.log(epsilon * thetak);
		return logHR;
	} // split

	private double getParabolicDistributionSample(double upperLimit) {
		double u1 = Randomizer.nextDouble();
		double sample =  upperLimit * Math.pow(u1, 1.0/3.0); 
		return sample;
	}
	
	public static void main(String[] args) {
		// test parabolic distribution sampling
		MergeSplitSpeciesTree t = new MergeSplitSpeciesTree();
		List<Double> list = new ArrayList<Double>();
		int LIMIT = 10000;
		int LIMIT100 = LIMIT / 100;
		for (int i = 0; i < LIMIT; i++) {
			list.add(t.getParabolicDistributionSample(1.0));
		}
		java.util.Collections.sort(list);
		System.out.println(list);
		double mse = 0;
		for (int i = 0; i < 100; i++) {
			double sum = 0;
			for (int j = 0; j < LIMIT100; j++) {
				sum += list.get(i*LIMIT100+j);
			}
			double diff = (sum/(LIMIT100)- Math.pow((i+0.5)/100.0, 1.0/3.0));
			mse += diff * diff;
		}
		System.out.println("Mean square error: " + mse);
	}

}


// TODO: unit test for getParabolicDistributionSample
// TODO: unit test for split and merge moves
// TODO: integration test for split/merge moves -- sampling from prior
// TODO: TreeLogAnalyser for clusterings -- this is probably a bit more work than I anticipated. The code is implemented in DensiTree, but not in BEAST yet.
// TODO: update prior to take in account gammas on zero length branches (those need not taken in account in the prior)
// TODO: update operators to take in account gammas on zero length branches (those need not be sampled) NodeBudger should not select nodes with children that have zero length branches, etc..
// TODO: update treelikelihood to replace convolution on zero branch lengths with something more efficient
// TODO: unit tests for each of the above updates
