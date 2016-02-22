package snap.likelihood;

import beast.core.Citation;
import beast.core.Description;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

@Description("Implementation of Felsentstein's theshold model")
@Citation("")
public class ThresholdTreeLikelihood extends GenericTreeLikelihood {

	
	Tree tree;
	// rate model, specifies rate for each branch, constant 1.0 by default
	BranchRateModel branchRateModel;
	Alignment data;
	// site model, like gamma heterogeneity model
	SiteModel siteModel;
	
	// contain likelihoods for each pattern
	double [] patternLogLikelihoods;
	
    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

	
	@Override
	public void initAndValidate() {
		tree = (Tree) treeInput.get();
		data = dataInput.get();
		siteModel = (SiteModel) siteModelInput.get();
        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        
        patternLogLikelihoods = new double[data.getPatternCount()];
        // TODO: initialise data structures that sets up data for leaves of the tree
        
        hasDirt = Tree.IS_FILTHY;
	}
	

	/**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLogP() {
    	traverse(tree.getRoot());
    	
    	logP = 0;
        for (int i = 0; i < data.getPatternCount(); i++) {
            logP += patternLogLikelihoods[i] * data.getPatternWeight(i);
        }
        return logP;
    }
	
    
    int traverse(final Node node) {
        int update = (node.isDirty() | hasDirt);

        final int iNode = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[iNode])) {
        	// set up site model for this branch
            m_branchLengths[iNode] = branchTime;
            final Node parent = node.getParent();
            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
                // TODO set up site model for this category for this branch
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                update |= (update1 | update2);

                if (siteModel.integrateAcrossCategories()) {
                    calculatePartials(childNum1, childNum2, iNode);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods

                    final double[] proportions = siteModel.getCategoryProportions(node);
                    integratePartials(node.getNr(), proportions);
                }

            }
        }
        return update;
    }


	// for each pattern
    // calculate patternLogLikelihoods, integrating over categories (if any) 
	private void integratePartials(int nr, double[] proportions) {
		// TODO 
		
	}

	// for each pattern
	// calculate partials for node iNode given the partials/states of children childNum1 and childNum2
	private void calculatePartials(int childNum1, int childNum2, int iNode) {
		// TODO  
	}
	
}
