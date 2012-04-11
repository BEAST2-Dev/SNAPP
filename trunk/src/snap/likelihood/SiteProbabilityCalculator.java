
/*
 * File SiteProbabilityCalculator.java
 *
 * Copyright (C) 2010 Remco Bouckaert, David Bryant remco@cs.auckland.ac.nz
 *
 * This file is part of SnAP.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * SnAP is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  SnAP is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SnAP; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */
package snap.likelihood;

import java.util.Arrays;

import snap.FMatrix;
import snap.NodeData;
import snap.matrix.QMatrix;

public class SiteProbabilityCalculator {

    void convolution2Dfft( double[][]A,  double[][]B, double[][] result) {

        //int ma,mb,na,nb;
        //getSize(A,ma,na);
        int na = A.length - 1; //A.getNrOfRows();
        int ma = na;//A.getNrOfCols();
        //getSize(B,mb,nb);
        int nb = B.length - 1; // B.getNrOfRows();
        int mb = nb;//B.getNrOfCols();

        //Determine smallest FFT array that is larger than A and B.
        int mr = ma+mb-1;
        int nr = na+nb-1;

        //result.resize(mr+1, nr+1);
        //double[][]result = new double[mr+1][nr+1];

        //typedef vector<double>::const_iterator DBL_PTR;

        for(int i3=1;i3<=mr;i3++) {
            for(int j3=1;j3<=nr;j3++) {
                double sum=0.0;
                //Sum over all i1,i2 such that i1+i2 = i3+1.
                //and over all j1,j2 such that j1+j2 = j3+1.
                //Note: 1<=i1<=ma  and 1<=i2<=mb
                int min_i = (i3>mb) ? (i3+1) - mb : 1;
                int max_i = (i3<ma) ? i3 : ma;

                int min_j = (j3>nb) ? (j3+1) - nb : 1;
                int max_j = (j3<na) ? j3 : na;

                for(int i1 = min_i;i1<=max_i;i1++) {
                    int i2 = (i3 + 1) - i1;
                    for(int j1 = min_j;j1<=max_j;j1++) {
                        int j2 = (j3+1) - j1;
                        sum+=A[i1][j1]*B[i2][j2];
                    }
                }
                result[i3][j3] = sum;
            }
        }
    } // convolution2Dfft

    /**
     Determines a non-zero right e-vector for the matrix Q, defined by u,v,coalescenceRate and N.
     The e-vector is normalised so that the entries for each n sum to 1.0

     //TODO: incorporate into code for abstract matrix
     */
    double [][]  findRootProbabilities(int N, double u, double v, double coalescenceRate, boolean dprint) throws Exception {
        double [][] x;
        QMatrix Qt = new QMatrix(N,u,v,coalescenceRate);
        double [] xcol;
        xcol = Qt.findOrthogonalVector(dprint);
        if (dprint) {
            System.out.println("xcol = " +Arrays.toString(xcol));
        }

        int index = 1;
        x = new double[N+1][];
        for(int n=1;n<=N;n++) {
            x[n] = new double[n+1];
            double rowsum = 0.0;
            for(int r=0;r<=n;r++) {
                double xcol_index = Math.max(xcol[index], 0.0);
                rowsum += xcol_index;
                x[n][r] = xcol_index;
                index++;
            }
            for(int r=0;r<=n;r++)
                x[n][r] = x[n][r] / rowsum;
        }
        return x;
    } // findRootProbabilities

    double doRootLikelihood(NodeData rootData, double u, double v, Double [] coalescenceRate, boolean dprint) throws Exception
    {
        //int N=rootData.m_n;
        int N=rootData.getSize();
		
		//System.err.println("Coalescent Rate at Root = "+coalescenceRate[rootData.getNr()]);
		
		
        double[][] conditional = findRootProbabilities(N, u, v, coalescenceRate[rootData.getNr()], dprint);

		//System.err.println("Root theta = "+rootData.gamma());
		

		if (dprint) {
			for(int n=1;n<=N;n++) {
//				if (rootData.m_Nb[n]<0.0)
//					System.out.println("Nb["+n+"] = "+rootData.m_Nb[n]);
				for(int r=0;r<=n;r++) {
					if (conditional[n][r]<0.0)
						System.out.println("conditional["+n+", "+r+"] = "+conditional[n][r]);;
					if (rootData.getFb().get(n,r)<0.0)
						System.out.println("Fb["+n+", "+r+"] = "+rootData.getFb().get(n,r));;
                }
                //cout+endl;
            }
        }

        double sum = 0.0;
        for(int n=1;n<=N;n++) {
            for(int r=0;r<=n;r++) {
                double term =  conditional[n][r] * rootData.getFb().get(n,r);
                sum += term;
                if (sum<0.0)
                    System.out.println("Numerical problems");
            }
        }
        return sum;
    } // doRootLikelihood


    static
    FCache m_cache;

    public void clearCache(int nNodeNrMax, int nRedsMax) {
        //m_cache = new FCache(nNodeNrMax, nRedsMax + 1);
        m_cache = new FCacheT(nNodeNrMax, nRedsMax + 1);
    }

    void doCachedInternalLikelihood(NodeData u1, NodeData u2, NodeData parent) {
        FCache.CacheObject o = m_cache.getBottomOfBrancheF(u1.getCacheIDT(), u2.getCacheIDT(), u1, u2, parent, this);
        parent.setCacheIDB(o.m_nCacheID);
        //parent.assignFb(o.getF());
        parent.initFb(o.getF());
    }

    void doCachedLeafLikelihood(NodeData node, int numReds, int totalCount, boolean bHasDominantMarkers) {
        FCache.CacheObject o = m_cache.getLeafF(node, numReds, totalCount, bHasDominantMarkers, this);
        node.setCacheIDB(o.m_nCacheID);
        //node.assignFb(o.getF());
		
		
		
		
		
        node.initFb(o.getF());
		
		//System.err.print(">>>>Leaf = "+node.getNr());
		//System.err.print(" ["+numReds+","+node.m_n+"::");
		//for(int k=0;k<=node.m_n;k++)
		//	System.err.print(""+node.getFb().get(node.m_n,k)+", ");
		//System.err.println(); 
		
    }

    void doCachedTopOfBranchLikelihood(NodeData node, double u, double v, Double [] coalescenceRate) throws Exception {
        FCache.CacheObject o = m_cache.getTopOfBrancheF(node.getCacheIDB(), node, u, v, coalescenceRate, this);
        node.setCacheIDT(o.m_nCacheID);
        //node.assignFt(o.getF());
        node.initFt(o.getF());
    }

    double [] getFt(int N, NodeData u) {
        double [] uFt = u.getFt().asVectorCopy();
        for(int n=1; n<=N; n++) {
            double b_nr = 1.0;
            for(int r=0;r<=n;r++) {
                uFt[n*(n+1)/2-1+r] *= b_nr;
                b_nr *= ((double)n - r)/(r+1);
            }
        }
        return uFt;
    }

    double [] getConvolution(int N1, int N2, double [] u1Ft, double [] u2Ft) {
        int N = N1 + N2;
        double [] parentFb = new double [(N+1)*(N+2)/2-1];
        for(int n1=1;n1<=N1;n1++) {
            for(int r1=0;r1<=n1;r1++) {
                double f11  =  u1Ft[n1*(n1+1)/2-1+r1];//u1.getFt().get(n1,r1);
                for(int n2=1;n2<=N2;n2++) {
                    for(int r2=0;r2<=n2;r2++) {
                        //parent.F[n1+n2][r1+r2] += u1.F[n1][r1] * u2.F[n2][r2];
                        //parent.getFb().add(n1+n2,r1+r2, f11 * u2.getFt().get(n2,r2));
                        //FAsVector[(n1+n2)*(n1+n2+1)/2-1+(r1+r2)] += f11 * u2.getFt().get(n2,r2);
                        parentFb[(n1+n2)*(n1+n2+1)/2-1+(r1+r2)] += f11 * u2Ft[n2*(n2+1)/2-1+r2];
                    }
                }
            }
        }
        //}
        return parentFb;
    }

    void doInternalLikelihood(NodeData u1, NodeData u2, NodeData parent, boolean dprint) {
        /**
         Let Y and Z be the F tables for u1 and u2.

         First compute Y <- Y .* B ; Z <- Z.* B

         where B is the the table of binomial values.

         we then have   F .* B  = Y conv Z

         where conv means a 2D convolution.
         **/
        // binom(n,r+1) = n! / (r+1)! / (n-r-1)! = binom(n,r) * (n-r)/(r+1)

        //Construct table of binomials

//        int N1 = u1.m_n;
//        int N2 = u2.m_n;

        int N1 = u1.getSize();
        int N2 = u2.getSize();


//		double [] u1Ft = getFt(N1, u1);
//		double [] u2Ft = getFt(N2, u2); 
//		double [] parentFb = getConvolution(N1, N2, u1Ft, u2Ft);


        double [] u1Ft = u1.getFt().asVectorCopy();
        for(int n=1; n<=N1; n++) {
            double b_nr = 1.0;
            for(int r=0;r<=n;r++) {
                //u2.getFt().mul(n,r, b_nr);
                u1Ft[n*(n+1)/2-1+r] *= b_nr;
                b_nr *= ((double)n - r)/(r+1);
            }
        }

        double [] u2Ft = u2.getFt().asVectorCopy();
        for(int n=1; n<=N2; n++) {
            double b_nr = 1.0;
            for(int r=0;r<=n;r++) {
                //u2.getFt().mul(n,r, b_nr);
                u2Ft[n*(n+1)/2-1+r] *= b_nr;
                b_nr *= ((double)n - r)/(r+1);
            }
        }


        //if (N1+N2>60) {
        //	//When problem is large, use fast fourier transform
        //	convolution2Dfft(u1.F, u2.F, parent.F);
        //} else {
        //otherwise, compute convolution directly.

        //parent.resizeF(N1+N2);
        int N = N1 + N2;
        double [] parentFb = new double [(N+1)*(N+2)/2-1];
        for(int n1=1;n1<=N1;n1++) {
            for(int r1=0;r1<=n1;r1++) {
                double f11  =  u1Ft[n1*(n1+1)/2-1+r1];//u1.getFt().get(n1,r1);
                for(int n2=1;n2<=N2;n2++) {
                    for(int r2=0;r2<=n2;r2++) {
                        //parent.F[n1+n2][r1+r2] += u1.F[n1][r1] * u2.F[n2][r2];
                        //parent.getFb().add(n1+n2,r1+r2, f11 * u2.getFt().get(n2,r2));
                        //FAsVector[(n1+n2)*(n1+n2+1)/2-1+(r1+r2)] += f11 * u2.getFt().get(n2,r2);
                        parentFb[(n1+n2)*(n1+n2+1)/2-1+(r1+r2)] += f11 * u2Ft[n2*(n2+1)/2-1+r2];
                    }
                }
            }
        }
        //}






//		parent.getF().setZero();
		for(int n=1;n<=N1+N2;n++) {
//			if (false && parent.m_Nb[n]==0) {
//				//parent.getFb().setZero(n);
//				Arrays.fill(parentFb, n*(n+1)/2-1, (n+1)*(n+2)/2-1, 0.0);
////				//Likelihood is zero if there can't be n lineages (e.g. n=1 when there are two children)
//			} else {
				double b_nr = 1.0;
				for(int r=0;r<=n;r++) {
					//double Fnr = parent.getFb().get(n,r);
					//Fnr /=  b_nr;
					//Fnr = Math.max(Fnr,0.0); //TODO: Fix this dodgy fix!!!!!!!
					//parent.getFb().set(n,r, Fnr);
					double Fnr = parentFb[n*(n+1)/2-1+r];
					Fnr /=  b_nr;
					Fnr = Math.max(Fnr,0.0); //TODO: Fix this dodgy fix!!!!!!!
					parentFb[n*(n+1)/2-1+r] = Fnr;
					b_nr *= ((double)n - r)/(r+1);
				}
			}
//		}
		parent.initFb(N1+N2, parentFb);
	}





    void doInternalLikelihood(NodeData u1, NodeData parent, boolean dprint)
    {
        /**
         This node has a single child. We simply copy the likelihood table from the
         top of the child branch.
         **/

        //int N = u1.m_n;
        parent.assignFt(u1.getFb());
        /*
          parent.resizeF(N);
          for(int n=1;n<=N;n++) {
              System.arraycopy(u1.F[n], 0, parent.F[n],0, u1.F[n].length);
          }
          */
    }


    /**
     Computes likelihood at a leaf. That is, one for the correct number of lineages and zero otherwise.
     **/
    void doLeafLikelihood(NodeData node, int nReds, int nTotal, boolean bHasDominantMarkers, boolean dprint)
    {
        
		if (bHasDominantMarkers && nReds>0) {
			
			//Need to account for the fact that the markers are dominant.
			
//			node.initFb(node.m_n, nReds);
//			node.getFb().set(node.m_n, nReds, 0);
//			int n = node.m_n/2; //Sample size in # individuals (rather than number of gametes)

			node.initFb(node.getSize(), nReds);
			node.getFb().set(node.getSize(), nReds, 0);
			int n = node.getSize()/2; //Sample size in # individuals (rather than number of gametes)
			//Compute p(r,k,n), which is the probability of r individuals having at least one copy of the 1 allele, 
			// given that k of their gameters carry the allele.
			double p_r_k_n = 1.0; //p(0,0,n)=1
			for (int r = 1;r<=nReds;r++)
				p_r_k_n = (p_r_k_n*2.0*(n-r+1.0))/(2.0*n-r+1.0);
			//Now p_r_k_n = p(r,r,n)
			
			for(int k=nReds;k<=2*nReds;k++) {
				if (k>nReds)
					p_r_k_n = (p_r_k_n * (2.0*nReds-k+1)*k) / (2.0*(k-nReds)*(2.0*n-k+1.0));
//				node.getFb().set(node.m_n,k,p_r_k_n);
				node.getFb().set(node.getSize(), k, p_r_k_n);
			}
			
		}
		else
//			node.initFb(node.m_n, nReds);
			node.initFb(node.getSize(), nReds);
		
		//System.err.print("Leaf = "+node.getNr());
		//System.err.print(" ["+nReds+","+node.m_n+"::");
		//for(int k=0;k<=node.m_n;k++)
		//	System.err.print(""+node.getFb().get(node.m_n,k)+", ");
		//System.err.println();
		
		
        //node.resizeF(node.n);
        //node.getFb().set(node.n,numReds,1.0);
    } // doLeafLikelihood

    /**

     Updates the likelihood at the top of a branch, assuming that the likelihood at the
     bottom of the branch as already been computed.
     * @throws Exception

     **/
    void doTopOfBranchLikelihood(NodeData node, double u, double v, Double [] coalescenceRate, boolean dprint) throws Exception {

        //int N = node.m_n;

        if (dprint) {
            System.err.print("BEFORE\t");
            System.err.println(node.getFt().toString());
        }

		//System.err.println("node.gamma = "+coalescenceRate[node.getNr()]+"\tnot t = "+node.t() + "\t2/node.gamma = "+(2.0/coalescenceRate[node.getNr()]));
		
        //FMatrix tmp = MatrixExponentiator.expQTtx(node.m_n, u, v, coalescenceRate[node.getNr()], node.t(), node.getFb());
        FMatrix tmp = MatrixExponentiator.expQTtx(node.getSize(), u, v, coalescenceRate[node.getNr()], node.t(), node.getFb());
        //TODO: What is the effect of the tolerance?

        node.initFt(tmp);
        //for(int i = 1;i <= node.n; i++) {
        //  System.arraycopy(tmp[i],0,node.F[i],0,tmp[i].length);
        //}


        if (dprint) {
            System.err.print("MIDWAY\t");
            System.err.println(node.getFt().toString());
        }
    } // doTopOfBranchLikelihood



    /**
     Compute the site likelihood of the tree.
     If updateAll is set to true, then all likelihoods are updated. Otherwise, only those leaf
     nodes with mustUpdate = true, or internal ancestors of those leaf nodes, are updated.
     @param tree  The species tree
     @param redCount  The number of red alleles in the sample: redCount[i] is the number for species i.
      * @throws Exception
     @bool updateAll Update the partial likelihoods for all nodes.
     */
    void computeSiteLikelihood2(NodeData tree, double u, double v, Double [] coalescenceRate, int [] redCount, int [] totalCount, boolean bHasDominantMarkers, boolean dprint/*=false*/) throws Exception {
        //Post-order traversal
        if (tree.isLeaf()) {
            doLeafLikelihood(tree, redCount[tree.getNr()], totalCount[tree.getNr()], bHasDominantMarkers, dprint);
        } else if (tree.getNrOfChildren() == 1) {
            NodeData p = tree.getChild(0);
            computeSiteLikelihood2(p, u, v, coalescenceRate, redCount, totalCount, bHasDominantMarkers, dprint);
            doTopOfBranchLikelihood(p, u,v, coalescenceRate, dprint);
            doInternalLikelihood(p, tree, dprint);
        } else { // assume two children
            NodeData leftChild = tree.getChild(0);
            NodeData rightChild = tree.getChild(1);
            computeSiteLikelihood2(leftChild, u, v, coalescenceRate, redCount, totalCount, bHasDominantMarkers, dprint);
            computeSiteLikelihood2(rightChild, u, v, coalescenceRate, redCount, totalCount, bHasDominantMarkers, dprint);
            doTopOfBranchLikelihood(leftChild, u,v,coalescenceRate, dprint);
            doTopOfBranchLikelihood(rightChild, u,v,coalescenceRate, dprint);
            doInternalLikelihood(leftChild, rightChild, tree, dprint);
        }
    } // computeSiteLikelihood2

    void computeCachedSiteLikelihood2(NodeData tree, double u, double v, Double [] coalescenceRate, int [] redCount, int [] totalCount, boolean bHasDominantMarkers, boolean dprint/*=false*/) throws Exception {
        //Post-order traversal
        if (tree.isLeaf()) {
            doCachedLeafLikelihood(tree, redCount[tree.getNr()], totalCount[tree.getNr()], bHasDominantMarkers);
        } else if (tree.getNrOfChildren() == 1) {
            NodeData p = tree.getChild(0);
            computeCachedSiteLikelihood2(p, u, v, coalescenceRate, redCount, totalCount, bHasDominantMarkers, dprint);
            doCachedTopOfBranchLikelihood(p, u,v, coalescenceRate);
            doInternalLikelihood(p, tree, false);
        } else { // assume two children
            NodeData leftChild = tree.getChild(0);
            NodeData rightChild = tree.getChild(1);
            computeCachedSiteLikelihood2(leftChild, u, v, coalescenceRate, redCount, totalCount, bHasDominantMarkers, dprint);
            computeCachedSiteLikelihood2(rightChild, u, v, coalescenceRate, redCount, totalCount, bHasDominantMarkers, dprint);
            doCachedTopOfBranchLikelihood(leftChild, u,v, coalescenceRate);
            doCachedTopOfBranchLikelihood(rightChild, u,v, coalescenceRate);
            doCachedInternalLikelihood(leftChild, rightChild, tree);
        }
    } // computeSiteLikelihood2

    public double computeSiteLikelihood(NodeData tree, double u, double v, Double [] coalescenceRate, int [] redCount, int [] totalCount, boolean bMutationOnlyAtRoot, boolean bHasDominantMarkers, boolean useCache, boolean dprint/*=false*/) throws Exception {
  
		//For some models, we want to allow mutation at the root but not along the branches. 
		double branch_u = u;
		double branch_v = v;
		if (bMutationOnlyAtRoot)
			branch_u = branch_v = 0.0;
				
		if (useCache) {
            computeCachedSiteLikelihood2(tree, branch_u, branch_v, coalescenceRate, redCount, totalCount, bHasDominantMarkers, dprint/*=false*/);
        } else {
            computeSiteLikelihood2(tree, branch_u, branch_v, coalescenceRate, redCount, totalCount, bHasDominantMarkers, dprint/*=false*/);
        }

        if (dprint)
            System.err.println(tree.getFb().toString());

        return doRootLikelihood(tree, u,v, coalescenceRate, dprint);
    } // computeSiteLikelihood


} // class SiteProbabilityCalculator
