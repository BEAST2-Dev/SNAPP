package beast.evolution.tree;

import java.util.List;

import beast.core.Description;
import beast.evolution.speciation.SpeciesTreeDistribution;

@Description("Yule prior that has a 1/X prior on the birth rate, which is analytically integrated out")
public class YulePriorOneOnXBirthRatePrior extends SpeciesTreeDistribution {

	   @Override
	    public void initAndValidate() throws Exception {
	        super.initAndValidate();
	        
	        // make sure that all tips are at the same height,
	        // otherwise this Yule Model is not appropriate
	        TreeInterface tree = treeInput.get();
	        if (tree == null) {
	            tree = treeIntervalsInput.get().treeInput.get();
	        }
	        List<Node> leafs = tree.getExternalNodes();
	        double height = leafs.get(0).getHeight();
	        for (Node leaf : leafs) {
	            if (Math.abs(leaf.getHeight() - height) > 1e-8) {
	                System.err.println("WARNING: Yule Prior cannot handle dated tips. Use for example a coalescent prior instead.");
	                break;
	            }
	        }
	    }

	    /**
	    * The Yule prior (as we use it) is proportional to

	    \lambda^n exp(-\lambda H) 

	    where n is the number of species and H is the sum of heights. Throwing on the 1/lambda prior we can integrate analytically:

	    \int_0^\infty \lambda^n exp(-\lambda H) 1/\lambda d \lambda

	    = ...

	    = (n-1)! / H^n
	    **/
	    @Override
	    public double calculateTreeLogLikelihood(final TreeInterface tree) {
	    	
	        final int n = tree.getLeafNodeCount();

	        double H = 0;
	        final Node[] nodes = tree.getNodesAsArray();
	        for (int i = n; i < nodes.length; i++) {
	            assert (!nodes[i].isLeaf());
		        final double height = nodes[i].getHeight();
		        H += height;
	        }

	        logP = Math.log(n-1) - n * Math.log(H);
	        return logP;
	    }

	    @Override
	    public boolean canHandleTipDates() {
	        return false;
	    }
}
