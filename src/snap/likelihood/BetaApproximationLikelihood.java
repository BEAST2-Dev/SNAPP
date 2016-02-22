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

	//TODO: Do the analysis to bound the error introduced by these.
	private static final double EPSILON = 1e-10;
	private static final double M_CUTOFF = 100;
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
	public void initAndValidate() {
		tree = (Tree) treeInput.get();
		data = (Data) dataInput.get();
		patternLogLikelihoods = new double[data.getPatternCount()];
		coalescenceRate = coalescenceRateInput.get();


		nPoints = 50; //We need to add a new parameter here. Note that nPoints higher -> more accurate, so we could make
		//approximations by keeping npoints small.
		int nNodes = tree.getNodeCount();
		partialIntegral = new double[nNodes][nPoints];

		//testSolve();
		
		double val = testBetaDensityQuadrature(nPoints,295, 203);

		
		hasDirt = Tree.IS_FILTHY;
	}


	/**
	 * Calculate the log likelihood of the current state.
	 *
	 * @return the log likelihood.
	 */
	@Override
	public double calculateLogP() {
		rootTheta = 2.0/coalescenceRate.getValue(tree.getRoot().getNr()); 

		//ToDO: the formula we really want is 
		//rootTheta = 2.0*mutationRate / coalescenceRate.getValue(tree.getRoot().getNr()); 
		//where
		//mutationRate = 2.0*u*v/(u+v)

		//System.err.println("Do we get here?");
		double logL=0.0;
		double sumL = 0.0;     //For debugging only. Sum of probabilities
		
		double[] patternLikelihoods = new double[data.getPatternCount()];
		
		for (thisPattern = 0; thisPattern < data.getPatternCount(); thisPattern++) {
			redCounts = data.getPattern(thisPattern);
			lineageCounts = data.getPatternLineagCounts(thisPattern);
			traverse(tree.getRoot());    //Compute likelihood and puts value in patternLogLikelihood[thisPattern]

			logL += patternLogLikelihoods[thisPattern] * data.getPatternWeight(thisPattern);   
			sumL += Math.exp(patternLogLikelihoods[thisPattern] * data.getPatternWeight(thisPattern));
			patternLikelihoods[thisPattern] = Math.exp(patternLogLikelihoods[thisPattern]);
		}
		System.err.println("Sum of probabilities = "+sumL);
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
		//X = [X0,X2,...,X{N-1}], X0 = 0, XN = 1


		int N = F.length;
		int[] Ns;
		int K = 7;


		double Bab = Math.exp(GammaFunction.lnGamma(alpha) + GammaFunction.lnGamma(beta) - GammaFunction.lnGamma(alpha+beta));


		//Divide points into blocks of size approx K=7 points.  
		if (N<=K) {
			Ns = new int[2];
			Ns[0] = 0; Ns[1] = N-1;
		} else {
			int m = (N-1)/(K-1);
			int r = (N-1)%(K-1);
			if (r<=1) {
				Ns = new int[m+1];
				for(int i=0;i<m;i++)
					Ns[i] = i*(K-1);
				Ns[m] = N-1;
			} else {
				Ns = new int[m+1];
				for(int i=0;i<=m/2;i++) 
					Ns[i]=(K-1)*i;
				Ns[m/2+1] = Ns[m/2]+r;
				for(int i=m/2+2;i<=m;i++)
					Ns[i]=Ns[i-1]+(K-1);					
			} 
		}

		double total = 0.0;
		int m = Ns.length-1;
		if (m/2>0) {
			for (int i=0;i<m/2;i++) {
				double[] Xlocal = new double[Ns[i+1]-Ns[i]+1];
				double[] Flocal = new double[Ns[i+1]-Ns[i]+1];
				for (int j=Ns[i];j<=Ns[i+1];j++) {
					double x = X[j];
					double f = F[j];
					double newf = f*Math.pow(1-x, beta-1)/Bab;
					Xlocal[j-Ns[i]] = x;
					Flocal[j-Ns[i]] = newf;
				}
				double tot = alphaKernelQuadrature(Xlocal,Flocal,alpha);
				total = total+tot;
			}  		
		}
		for (int i=m/2;i<=m-1;i++) {
			double[] Xlocal = new double[Ns[i+1]-Ns[i]+1];
			double[] Flocal = new double[Ns[i+1]-Ns[i]+1];
			for (int j=Ns[i];j<=Ns[i+1];j++) {
				double x = X[j];
				double f = F[j];
				double newf = f*Math.pow(x, alpha-1)/Bab;
				Xlocal[j-Ns[i]] = x;
				Flocal[j-Ns[i]] = newf;
			}
			double tot = betaKernelQuadrature(Xlocal,Flocal,beta);
			total = total+tot;
		}

		return total;
	}

	private static double erf(double z) {  //TODO get a more accurate erf
		double t = 1.0 / (1.0 + 0.5 * Math.abs(z));

		// use Horner's method
		double ans = 1 - t * Math.exp( -z*z   -   1.26551223 +
				t * ( 1.00002368 +
						t * ( 0.37409196 + 
								t * ( 0.09678418 + 
										t * (-0.18628806 + 
												t * ( 0.27886807 + 
														t * (-1.13520398 + 
																t * ( 1.48851587 + 
																		t * (-0.82215223 + 
																				t * ( 0.17087277))))))))));
		if (z >= 0) return  ans;
		else        return -ans;
    }
	

	/**
	 * Given points X = [X0,...,XN], where X0=0, XN=1, and function evaluations F[i] = f(X[i])
	 * estimates integral
	 *     \int K(x;mu,v) f(x) dx
	 *     where K(x;mu,v) is a (peaked) Gaussian Kernel
	 *
	 * Uses interpolatory quadrature.
	 * @param X
	 * @param F
	 * @param alpha
	 * @param beta
	 * @return
	 */
	private double gaussianDensityQuadrature(double[] X, double[] F, double mu, double v) {
	
		/**
		 * Solves:
		 * 
		 * 
		 */
	
		int N = F.length;
		double total = 0.0;
		double root2 = Math.sqrt(2);
		double rootv = Math.sqrt(v);
		for (int i=0;i<N-1;i++) {
			double z1 = 0.5*erf(0.5*root2*(X[i+1]-mu)/rootv) - 0.5*erf(0.5*root2*(X[i]-mu)/rootv);
			double z2 = 0.5*Math.exp(0.5*v-mu)*(erf(0.5*root2*(X[i+1]-mu+v)/rootv) - erf(0.5*root2*(X[i]-mu+v)/rootv) );
			double diff = Math.exp(-X[i+1])-Math.exp(-X[i]);
			double w1 = (Math.exp(-X[i+1])*z1 - z2)/diff;
			double w2 = (-Math.exp(-X[i])*z1 + z2)/diff;
			total = total + w1*F[i] + w2*F[i+1];
		}
		return total;	
	}
		
		/**
		 * Implementation of Simpsons rules on [0,1].
		 * @param F
		 * @return
		 */
	private double simpsons(double[] F) {
		int N = F.length;
		double total = 0.0;
		double h = 1.0/(N-1.0);
		for (int i=1;i<N-1;i+=2) {
			total = total + h/3*F[i-1] + 4*h/3*F[i] + h/3*F[i+1];
		}
		if (N%2==0) {
			//Odd number of intervals. Interpolate through F[N-3]...F{N-1] and 
			//integrate last interval.
			total = total - h/12 * F[N-3] + 2.0*h/3 * F[N-2] + 5.0/12 * F[N-1];
		}
		return total;
	}
	
	private double testBetaDensityQuadrature(int N,double alpha, double beta) {
		double[] F = new double[N];
		double[] x = new double[N];
		
		for(int i=0;i<N;i++) {
			F[i] = 1.0;
			x[i] = 1.0/(N-1.0) * i;
		}
		return betaDensityQuadrature(x,F,alpha,beta);
		
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

				if (node.isRoot())      
					patternLogLikelihoods[thisPattern] = Math.log(partialIntegral[iNode][0]);              	
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
			xvals[j] = (1.0/(nPoints-1))*j; 

		for (int j=0;j<nPoints;j++) 				
			fvals[j] = partialIntegral[childNum1][j]*partialIntegral[childNum2][j];

		if (node.isRoot()) {
			double alpha = rootTheta/(1-2.0*rootTheta);
			double beta = alpha;
			if (alpha*alpha<2*EPSILON) 
				partialIntegral[iNode][0] = 0.5*fvals[0]+0.5*fvals[nPoints-1];
			else
				partialIntegral[iNode][0] = betaDensityQuadrature(xvals,fvals,alpha,beta); 
		}
		else {
			double coalRate=coalescenceRate.getValue(node.getNr());
			double delta = 1.0 - Math.exp(-2.0*coalRate * node.getLength());
			double m = (1.0 - delta)/delta;
			for (int s=0;s<nPoints;s++) {
				if (delta<EPSILON) {
					partialIntegral[iNode][s] = fvals[s]; //Zero branch length - copy probability.
				} else if (xvals[s]<EPSILON || xvals[s]>1.0-EPSILON || delta>1.0-EPSILON) {   //Saturation/Fixation with expected frequency xvals[s]. 
					partialIntegral[iNode][s] = xvals[s]*fvals[nPoints-1] + (1.0-xvals[s])*fvals[0]; 
				} else if (m>M_CUTOFF) {
					double mu = xvals[s];
					double v = xvals[s]*(1-xvals[s])/(m+1);
					partialIntegral[iNode][s] = gaussianDensityQuadrature(xvals,fvals,mu,v);					
				} else {
					double alpha = xvals[s]*(1-delta)/delta;
					double beta = (1-xvals[s])*(1-delta)/delta;
					partialIntegral[iNode][s] = betaDensityQuadrature(xvals,fvals,alpha,beta); 
				}
			}			
		}
	}
	
	
	private void calculateLeafPartials(Node node) {
		double[] xvals = new double[nPoints];
		double[] fvals = new double[nPoints];

		int n = lineageCounts[node.getNr()];
		int r = redCounts[node.getNr()];

		for (int j=0;j<nPoints;j++)
			xvals[j] = (1.0/(nPoints-1))*j; 

		double logBinom = Binomial.logChoose(n, r);
		for (int j=0;j<nPoints;j++) {
			fvals[j] = Math.exp(logBinom + r*Math.log(xvals[j]) + (n-r)*Math.log(1.0-xvals[j]));
		}
		if (r!=0)
			fvals[0]=0.0;  //Probability of observing > 0 is 0 when prob = 0
		else
			fvals[0]=1.0;  //Probability of observing  0 is 1 when prob = 0
		if (r!=n)
			fvals[nPoints-1]=0.0; //Probability of observing < n is 0 when prob = 1
		else
			fvals[nPoints-1]=1;   //Probability of observing  n is 1 when prob = 1

		double coalRate=coalescenceRate.getValue(node.getNr());
		double delta = 1.0 - Math.exp(-2.0*coalRate * node.getLength());
		double m = (1.0 - delta)/delta;
		
		for (int s = 0;s<nPoints;s++) {			
			if (delta<EPSILON) {
				partialIntegral[node.getNr()][s] = fvals[s]; //Zero branch length - copy probability.
			} else if (xvals[s]<EPSILON || xvals[s]>1.0-EPSILON || delta>1.0-EPSILON) {   //Saturation/Fixation with expected frequency xvals[s]. 
				partialIntegral[node.getNr()][s] = xvals[s]*fvals[nPoints-1] + (1.0-xvals[s])*fvals[0]; 
			} else if (m>M_CUTOFF) {
				double mu = xvals[s];
				double v = xvals[s]*(1-xvals[s])/(m+1);
				partialIntegral[node.getNr()][s] = gaussianDensityQuadrature(xvals,fvals,mu,v);					
			} else {
				double alpha = xvals[s]*(1-delta)/delta;
				double beta = (1-xvals[s])*(1-delta)/delta;
				partialIntegral[node.getNr()][s] = betaDensityQuadrature(xvals,fvals,alpha,beta); 
			}
		}
	}




}
