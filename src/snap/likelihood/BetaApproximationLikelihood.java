package snap.likelihood;


import snap.Data;
import beast.core.Citation;
import beast.core.Description;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.Binomial;
import beast.math.GammaFunction;

import org.apache.commons.math.distribution.BetaDistribution;

@Description("Implementation of Siren et al. Beta SNP approximation using Hiscott et al. integrator")
@Citation("")
public class BetaApproximationLikelihood extends GenericTreeLikelihood {

	
	Tree tree;
	Double rootTheta; //Parameter for root distribution. 
	Data data;	
	
	int nPoints;    //Size of mesh used in integration
	int thisPattern; 
	protected double[][] partialIntegral;  //Partial likelihoods. First index = node, second index = point number.
	protected int[] lineageCounts; //Lineage counts for each leaf taxa, for current pattern.
	protected int[] redCounts; //Red allele counts for each leaf, for current pattern.
	protected double[] patternLogLikelihoods;
	
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
	public void initAndValidate() throws Exception {
		tree = (Tree) treeInput.get();
		data = dataInput.get();
		patternLogLikelihoods = new double[data.getPatternCount()];
		
        nPoints = 20; //TODO: read this from the xml
        int nNodes = tree.getNodeCount();
        partialIntegral = new double[nPoints][nNodes];
        rootTheta = 0.5; //We should get this from the PARAMETERS.
        
        hasDirt = Tree.IS_FILTHY;
	}
	

	/**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLogP() {
    	
    	
    	double logL=0.0;
    	
        for (int thisPattern = 0; thisPattern < data.getPatternCount(); thisPattern++) {
        	redCounts = data.getPattern(thisPattern);
        	lineageCounts = data.getPatternLineagCounts(thisPattern);
        	traverse(tree.getRoot());    //Compute likelihood and puts value in patternLogLikelihood[thisPattern]
        	
            logL += patternLogLikelihoods[thisPattern] * data.getPatternWeight(thisPattern);   
        }
        return logL;
    }
	
   
    
    /**
     * Evaluate a Simpson's type rule on the interval x[0]...x[n-1], where f[i] = f(x[i])
     * @param x  array of x values. Assumed to be in strictly increasing order, but need not be regular.
     * @param f  array of f values with same length as x. 
     * @return  Simpson's estimate of integral of f between x[0] and x[n-1]. 
     * If the x[i] are regular and n is even then this returns the standard Simpson's estimate. Otherwise, it 
     * divides the intervals into consecutive pairs. For each pair  [a,c][c,b] it determines the integral
     * which is exact for all quadratics (as in Simpson's).
     * If n is odd then the last interval is evaluated using the the same rule but to evalute the
     * last integral using the final three points.
     *     
     */
    double simpsons(double[] x, double[] f) {
    	int n = x.length;
    	double total = 0.0;
    	for (int i = 0;i<n-2;i=i+2) {
    		double a,b,c,wa,wb,wc;
    		a = x[i]; c = x[i+2]; b = x[i+1]; 
    		wa = (b-a)*(2.0*a + b - 3.0*c)/(c-a)/6.0;
    		wb = (b-a)*(a + 2.0*b - 3.0*c)/(c-b)/6.0;
    		wc = (b-a)*(a-b)*(a-b)/(b-c)/(c-a)/6.0;
    		total = total + f[i]*wa + f[i+1]*wc+f[i+2]*wb;
    	}
    	if (n%2!=0){
    		double a,b,c,wa,wb,wc;
    		a = x[n-2]; c = x[n-3]; b = x[n-1]; 
    		wa = (b-a)*(2.0*a + b - 3.0*c)/(c-a)/6.0;
    		wb = (b-a)*(a + 2.0*b - 3.0*c)/(c-b)/6.0;
    		wc = (b-a)*(a-b)*(a-b)/(b-c)/(c-a)/6.0;
    		total = total + f[n-2]*wa + f[n-3]*wc+f[n-1]*wb;
    	}
    	return total;
    }
    
    
    
    int traverse(final Node node) {
        int update = (node.isDirty() | hasDirt| Tree.IS_DIRTY);

        final int iNode = node.getNr();
       

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

                calculatePartials(childNum1, childNum2, node);
                
                if (node.isRoot()) {     
                	double l = 0.0;
                	for (int s=0;s<nPoints;s++)
                		l += partialIntegral[iNode][s];
                	patternLogLikelihoods[thisPattern] = Math.log(l);              	
                }

            }
        } else {
        	//Get the data for this node, and compute the partial likelihoods.
        	calculateLeafPartials(node);
        }
        return update;
    }


	

	// for each pattern
	// calculate partials for node iNode given the partials/states of children childNum1 and childNum2
	private void calculatePartials(int childNum1, int childNum2, Node node) {
		double[] xvals = new double[nPoints];
		double[] fvals = new double[nPoints];
		int iNode = node.getNr();
		
		for (int j=0;j<nPoints;j++)
			xvals[j] = (1.0/(nPoints-1))*j; //TODO: Use GaussQ points
		
		for (int s=0;s<nPoints;s++) {
			if (node.isRoot()) {
				for (int j=0;j<nPoints;j++) {					
					fvals[j] = rootProb(xvals[j]) * partialIntegral[childNum1][j]*partialIntegral[childNum2][j];
				}
			} else {		
				for (int j=0;j<nPoints;j++) {			
					fvals[j] = transProb(xvals[s],xvals[j],node.getLength()) * partialIntegral[childNum1][j]*partialIntegral[childNum2][j];
				}
			}
			partialIntegral[iNode][s] = simpsons(xvals,fvals); 
		}			
	}
	
	private void calculateLeafPartials(Node node) {
		double[] xvals = new double[nPoints];
		double[] fvals = new double[nPoints];
		double[] lvals = new double[nPoints]; //Stores probability of data given leaf value.
				
		int n = lineageCounts[node.getNr()];
		int r = redCounts[node.getNr()];
		
		
		for (int j=0;j<nPoints;j++)
			xvals[j] = (1.0/(nPoints-1))*j; //TODO: Use GaussQ points
		
		double logBinom = Binomial.logChoose(n, r);
		for (int j=0;j<nPoints;j++) {
			lvals[j] = Math.exp(logBinom + r*xvals[j] + (n-r)*(1.0-xvals[j]));
		}
		for (int s = 0;s<nPoints;s++) {
			for (int j=0;j<nPoints;j++) {
				fvals[j] = transProb(xvals[s],xvals[j],node.getLength())*lvals[j];
			}
		}
	}
	
	private double betaDensity(double x, double alpha, double beta) {
		double logP = GammaFunction.lnGamma(alpha+beta) - GammaFunction.lnGamma(alpha) - GammaFunction.lnGamma(beta);
		logP += (alpha - 1)*Math.log(x) + (beta-1)*Math.log(1-x);
		return Math.exp(logP);
	}
		
	private double transProb(double x,double y, double t) {
		//From Siren et al 2011. Compute density of y, for a Beta with parameters x*(1-t)/t, (1-x)(1-t)/t.
		if (t==0)
			return (x==y)?1.0:0.0;		
		return betaDensity(y,x*(1-t)/t,(1-x)*(1-t)/t);
	}
	
	private double rootProb(double x) {
		return betaDensity(x,rootTheta,rootTheta);
	}
	
	
}
