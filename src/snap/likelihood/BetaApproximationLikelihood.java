package snap.likelihood;



import java.util.Arrays;

import snap.Data;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.Binomial;
import beast.math.GammaFunction;


@Description("Implementation of Siren et al. Beta SNP approximation using Hiscott et al. integrator")
@Citation("")
public class BetaApproximationLikelihood extends GenericTreeLikelihood {
	public Input<RealParameter> coalescenceRateInput = new Input<RealParameter>("coalescenceRate", "population size parameter with one value for each node in the tree");

	public BetaApproximationLikelihood() {
		siteModelInput.setRule(Validate.OPTIONAL);
	}

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

	RealParameter coalescenceRate;


	@Override
	public void initAndValidate() throws Exception {
		tree = (Tree) treeInput.get();
		data = (Data) dataInput.get();
		patternLogLikelihoods = new double[data.getPatternCount()];
		coalescenceRate = coalescenceRateInput.get();

		nPoints = 50; //We need to add a new parameter here. Note that nPoints higher -> more accurate, so we could make
						//approximations by keeping npoints small.
		int nNodes = tree.getNodeCount();
		partialIntegral = new double[nPoints][nNodes];

		testSolve();
		
		hasDirt = Tree.IS_FILTHY;
	}


	/**
	 * Calculate the log likelihood of the current state.
	 *
	 * @return the log likelihood.
	 */
	@Override
	public double calculateLogP() {
		rootTheta = coalescenceRate.getValue(tree.getRoot().getNr()); //We should get this from the PARAMETERS.

		System.err.println("Do we get here?");
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
	 * Solves Ax = b where A is square and assumed non-singular.
	 * Uses LU decomposition: Algorithm 3.4.1 of Golub and Van Loan edition 4.
	 * @param A
	 * @param b
	 * @return x
	 */

	private double[] solveMatrixEqn(double[][] A, double[] b) {  	
		int n = b.length;
		int[] piv = new int[n];

		for (int k=1;k<n;k++) {
			int mu = k;
			double rowMax = Math.abs(A[mu-1][k-1]);
			for (int i=k+1;i<=n;i++) {
				double Aik = Math.abs(A[i-1][k-1]);
				if (Aik > rowMax) {
					mu = i;
					rowMax = Aik;
				}
			}
			piv[k-1] = mu;	
			for (int i=1;i<=n;i++) {
				double x = A[k-1][i-1];
				A[k-1][i-1] = A[mu-1][i-1];
				A[mu-1][i-1] = x;
			}
			if (A[k-1][k-1] != 0.0) {
				for (int rho=k+1;rho<=n;rho++) {
					A[rho-1][k-1] = A[rho-1][k-1]/A[k-1][k-1];
					for(int rho2 = k+1;rho2<=n;rho2++) {
						A[rho-1][rho2-1] = A[rho-1][rho2-1] - A[rho-1][k-1]*A[k-1][rho2-1];
					}
				}
			}
		}

		//Overwrite b by Pb
		for (int k=1;k<n;k++) {
			double x = b[k-1];
			b[k-1] = b[piv[k-1]-1];
			b[piv[k-1]-1] = x;
		}

		double[] x = new double[n];
		double[] z = new double[n];
		//solve Lz = Pb, (note - could put x and z in same memory)
		for (int i=1;i<=n;i++) {
			z[i-1] = b[i-1];
			for (int j=1;j<i;j++)
				z[i-1] = z[i-1] - A[i-1][j-1] * z[j-1];
		}
			
		//solve Ux = z
		for (int i=n;i>=1;i--) {
			x[i-1] = z[i-1];
			for (int j=i+1;j<=n;j++)
				x[i-1] -= A[i-1][j-1]*x[j-1];
			x[i-1] = x[i-1]/A[i-1][i-1];	
		}
		
		return x;

	}


	private void testSolve() {
		double[][]A = {{3,17,10},{2,4,-2},{6,18,-12}};
		double[] b = {1,3,4};
		double[] x = solveMatrixEqn(A,b);
		for (int i=0;i<3;i++)
			System.out.println(x[i]);
	}

	/**
	 * Evaluates \int x^(alpha-1) f(x)  dx from X[0] ... X[n] by interpolating
	 * an nth degree polynomial.
	 * @param X
	 * @param fx
	 * @return
	 */
	private double alphaKernelQuadrature(double[] X, double[] fx, double alpha) {
		int N = X.length-1;
		double a = X[0], b = X[N];

		//z[k] =  \int_a^b x^(alpha-1) x^k dx 
		double[] z = new double[N+1];
		for (int k=0;k<N+1;k++)
			z[k] = (Math.pow(b, alpha+k) - Math.pow(a, alpha+k))/(alpha+k);

		double[][] M = new double[N+1][N+1];
		for (int k=0;k<N+1;k++)
			for(int i=0;i<N+1;i++)
				M[k][i] = Math.pow(X[i],k);

		double[] w = solveMatrixEqn(M,z);
		double total = 0.0;
		for (int k = 0;k<N+1;k++)
			total += fx[k]*w[k];
		return total;
	}


	/**
	 * Evaluates \int (1-x)^(beta-1) f(x)  dx from X[0] ... X[n] by interpolating
	 * an nth degree polynomial.
	 * @param X
	 * @param fx
	 * @return
	 */
	private double betaKernelQuadrature(double[] X, double[] fx, double beta) {
		int N = X.length-1;
		double a = X[0], b = X[N];

		//z[k] =  \int_a^b x^(alpha-1) x^k dx 
		double[] z = new double[N+1];
		for (int k=0;k<N+1;k++)
			z[k] = (Math.pow((1-a), beta+k) - Math.pow(1-b, beta+k))/(beta+k);

		double[][] M = new double[N+1][N+1];
		for (int k=0;k<N+1;k++)
			for(int i=0;i<N+1;i++)
				M[k][i] = Math.pow(1.0-X[i],k);

		double[] w = solveMatrixEqn(M,z);
		double total = 0.0;
		for (int k = 0;k<N+1;k++)
			total += fx[k]*w[k];
		return total;
	}

	/**
	 * Given points X = [X0,...,XN], where X0=0, XN=1, and function evaluations F[i] = f(X[i])
	 * estimates integral
	 *    1/Beta(alpha,beta) \int x^(alpha-1) (1-x)^(beta-1) f(x) dx
	 * Uses interpolatory quadrature, on sets of 8 points (or thereabouts).
	 * For the first part of the integral, uses quadrature rule of the form
	 * int x^(alpha-1) f(x) =  \sum w_i f_i
	 * and in the second part,
	 * int (1-x)^(beta-1) f(x) = \sum w_i f_i
	 * Thereby avoiding the potential singularities with alpha<1 or beta<1.
	 * 
	 * @param X
	 * @param F
	 * @param alpha
	 * @param beta
	 * @return
	 */
	private double betaDensityQuadrature(double[] X, double[] F, double alpha, double beta) {
		//X = [X0,X1,...,XN], X0 = 0, XN = 1
		double Bab = Math.exp(GammaFunction.lnGamma(alpha) + GammaFunction.lnGamma(beta) - GammaFunction.lnGamma(alpha+beta));
		int N = F.length-1;
		int[] Ns;
		int K = 7;

		//Divide points into blocks of size K=7.
		if (N<=K) {
			Ns = new int[2];
			Ns[0] = 0; Ns[1] = N;
		} else {
			int m = N/K;
			int r = N%K;
			if (r<=0) {
				Ns = new int[m+1];
				for(int i=0;i<m;i++)
					Ns[i]=K*i;
				Ns[m] = N;
			} else {
				Ns = new int[m+2];
				for(int i=0;i<m/2;i++)
					Ns[i]=K*i;
				Ns[m/2]=K*(m/2)+r;
				for(int i=m/2+1;i<=m;i++)
					Ns[i]=K*i+r;
			}
		}

		double total = 0.0;

		int m = Ns.length-1;
		if (m/2>0) {
			for (int i=0;i<m/2;i++) {
				double[] Xlocal = new double[Ns[i+1]-Ns[i]+2];
				double[] Flocal = new double[Ns[i+1]-Ns[i]+2];
				for (int j=Ns[i];j<=Ns[i+1];j++) {
					double x = X[j];
					double f = F[j];
					double newf = f*Math.pow(1-x, beta-1)/Bab;
					Xlocal[j-Ns[i]] = x;
					Flocal[j-Ns[i]] = newf;
					total = total+alphaKernelQuadrature(Xlocal,Flocal,alpha);
				}
			}  		
		}
		for (int i=m/2;i<=m-1;i++) {
			double[] Xlocal = new double[Ns[i+1]-Ns[i]+2];
			double[] Flocal = new double[Ns[i+1]-Ns[i]+2];
			for (int j=Ns[i];j<=Ns[i+1];j++) {
				double x = X[j];
				double f = F[j];
				double newf = f*Math.pow(x, alpha-1)/Bab;
				Xlocal[j-Ns[i]] = x;
				Flocal[j-Ns[i]] = newf;
				total = total+betaKernelQuadrature(Xlocal,Flocal,beta);
			}
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
					fvals[j] = partialIntegral[childNum1][j]*partialIntegral[childNum2][j];
				}
			} else {		
				for (int j=0;j<nPoints;j++) {
					fvals[j] = partialIntegral[childNum1][j]*partialIntegral[childNum2][j];
				}
			}

			double alpha,beta;
			if (node.isRoot()) {
				alpha = rootTheta/(1-2.0*rootTheta);
				beta = alpha;
			} else {
				double coalRate=coalescenceRate.getValue(node.getNr());
				double delta = 1.0 - Math.exp(-2.0*coalRate * node.getLength());
				alpha = xvals[s]*(1-delta)/delta;
				beta = (1-xvals[s])*(1-delta)/delta;
			}

			partialIntegral[iNode][s] = betaDensityQuadrature(xvals,fvals,alpha,beta); 
		}			
	}

	private void calculateLeafPartials(Node node) {
		double[] xvals = new double[nPoints];
		double[] fvals = new double[nPoints];

		int n = lineageCounts[node.getNr()];
		int r = redCounts[node.getNr()];

		for (int j=0;j<nPoints;j++)
			xvals[j] = (1.0/(nPoints-1))*j; //TODO: Use GaussQ points

		double logBinom = Binomial.logChoose(n, r);
		for (int j=0;j<nPoints;j++) {
			fvals[j] = Math.exp(logBinom + r*xvals[j] + (n-r)*(1.0-xvals[j]));
		}
		for (int s = 0;s<nPoints;s++) {
			double coalRate=coalescenceRate.getValue(node.getNr());
			double delta = 1.0 - Math.exp(-2.0*coalRate * node.getLength());
			double alpha = xvals[s]*(1-delta)/delta;
			double beta = (1-xvals[s])*(1-delta)/delta;
			partialIntegral[node.getNr()][s] = betaDensityQuadrature(xvals,fvals,alpha,beta);
		}
	}




}
