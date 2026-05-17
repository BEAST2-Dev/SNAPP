
/* *
 * Title:        <p>
 * Description:  <p>
 * Copyright:    Copyright (c) <p>
 * Company:      <p>
 * @author
 * @version 1.0
 */


package snap.distribution;

/**
 * Extends the class {@link DiscreteDistributionInt} for the
 * <EM>Poisson</EM> distribution with mean 
 * <SPAN CLASS="MATH"><I>&#955;</I>&nbsp;&gt;=&nbsp; 0</SPAN>.
 * The mass function is
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay">
 * <I>p</I>(<I>x</I>) = <I>e</I><SUP>-<I>&#955;</I></SUP><I>&#955;</I><SUP>x</SUP>/(<I>x</I>!),&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;for <I>x</I> = 0, 1,...
 * </DIV><P></P>
 * and the distribution function is
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay">
 * <I>F</I>(<I>x</I>) = <I>e</I><SUP>-<I>&#955;</I></SUP>&sum;<SUB>j=0</SUB><SUP>x</SUP> &nbsp;<I>&#955;</I><SUP>j</SUP>/(<I>j</I>!),&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;for <I>x</I> = 0, 1,....
 * </DIV><P></P>
 * If one has to compute <SPAN CLASS="MATH"><I>p</I>(<I>x</I>)</SPAN> and/or <SPAN CLASS="MATH"><I>F</I>(<I>x</I>)</SPAN> for several values of <SPAN CLASS="MATH"><I>x</I></SPAN>
 * with the same <SPAN CLASS="MATH"><I>&#955;</I></SPAN>, where <SPAN CLASS="MATH"><I>&#955;</I></SPAN> is not too large, then it is more
 * efficient to instantiate an object and use the non-static methods, since
 * the functions will then be computed once and kept in arrays.
 * 
 * <P>
 * For the static methods that compute <SPAN CLASS="MATH"><I>F</I>(<I>x</I>)</SPAN> and 
 * <SPAN CLASS="MATH">bar(F)(<I>x</I>)</SPAN>,
 * we exploit the relationship
 * 
 * <SPAN CLASS="MATH"><I>F</I>(<I>x</I>) = 1 - <I>G</I><SUB>x+1</SUB>(<I>&#955;</I>)</SPAN>, where <SPAN CLASS="MATH"><I>G</I><SUB>x+1</SUB></SPAN> is the <SPAN  CLASS="textit">gamma</SPAN>
 * distribution function with parameters 
 * <SPAN CLASS="MATH">(<I>&#945;</I>, <I>&#955;</I>) = (<I>x</I> + 1, 1)</SPAN>.
 * 
 */
public class PoissonDist extends DiscreteDistributionInt {

   private double lambda;



   public static double MAXLAMBDA = 100000;


   /**
    * Creates an object that contains
    *    the probability and distribution functions, for the Poisson
    *    distribution with parameter <TT>lambda</TT>, which are
    *    computed and stored in dynamic arrays inside that object.
    * 
    */
   public PoissonDist (double lambda) {
      setLambda (lambda);
   }


   public double prob (int x) {
      if (x < 0)
         return 0.0;
      if (pdf == null)
         return prob (lambda, x);
      if (x > xmax || x < xmin)
         return prob (lambda, x);
      return pdf[x - xmin];
   }

   public double cdf (int x) {
      double Sum = 0.0;
      int j;

      if (x < 0)
         return 0.0;
      if (lambda == 0.0)
         return 1.0;

      /* For large lambda, we use the Chi2 distribution according to the exact
         relation, with 2x + 2 degrees of freedom

         cdf (lambda, x) = 1 - chiSquare (2x + 2, 2*lambda)

         which equals also 1 - gamma (x + 1, lambda) */
      if (cdf == null)
         return GammaDist.barF (x + 1.0, 15, lambda);

      if (x >= xmax)
         return 1.0;

      if (x < xmin) {
         // Sum a few terms to get a few decimals far in the lower tail. One
         // could also call GammaDist.barF instead.
         final int RMAX = 20;
	 int i;
	 double term = prob(lambda, x);
	 Sum = term;
	 i = x;
	 while (i > 0 && i >= x - RMAX) {
	    term = term * i / lambda;
	    i--;
	    Sum += term;
	 }
	 return Sum;
      }

      if (x <= xmed)
         return cdf[x - xmin];
      else
         // We keep the complementary distribution in the upper part of cdf
         return 1.0 - cdf[x + 1 - xmin];
   }


   public double barF (int x) {
      /*
       * poisson (lambda, x) = 1 - cdf (lambda, x - 1)
       */

      if (x <= 0)
         return 1.0;

      /* For large lambda,  we use the Chi2 distribution according to the exact
         relation, with 2x + 2 degrees of freedom

         cdf (lambda, x) = 1 - ChiSquare.cdf (2x + 2, 2*lambda)
         cdf (lambda, x) = 1 - GammaDist.cdf (x + 1, lambda)
       */

      if (cdf == null)
         return GammaDist.cdf ((double)x, 15, lambda);

      if (x > xmax)
//         return GammaDist.cdf ((double)x, 15, lambda);
         return PoissonDist.barF(lambda, x);
      if (x <= xmin)
         return 1.0;
      if (x > xmed)
         // We keep the complementary distribution in the upper part of cdf
         return cdf[x - xmin];
      else
         return 1.0 - cdf[x - 1 - xmin];
   }

   public int inverseFInt (double u) {
      if (cdf == null)
         return inverseF (lambda, u);
      else
         return super.inverseFInt (u);
   }

   public double getMean() {
      return PoissonDist.getMean (lambda);
   }

   public double getVariance() {
      return PoissonDist.getVariance (lambda);
   }

   public double getStandardDeviation() {
      return PoissonDist.getStandardDeviation (lambda);
   }

   /**
    * Computes and returns the Poisson probability
    *   <SPAN CLASS="MATH"><I>p</I>(<I>x</I>)</SPAN> for <SPAN CLASS="MATH"><I>&#955;</I> =</SPAN> <TT>lambda</TT>..
    * 
    */
   public static double prob (double lambda, int x) {
      final double LAMBDALIM = 20.0;
      double y;
      double Res;

      if (x < 0)
         return 0.0;

      if (lambda >= 100.0) {
         if ((double) x >= 10.0*lambda)
            return 0.0;
      } else if (lambda >= 3.0) {
         if ((double) x >= 100.0*lambda)
            return 0.0;
      } else {
         if ((double) x >= 200.0*Math.max(1.0, lambda))
            return 0.0;
      }   

      if (lambda < LAMBDALIM && x <= 100)
         Res = Math.exp (-lambda)*Math.pow (lambda, x)/Num.factorial (x);
      else {
         y = x*Math.log (lambda) - Num.lnGamma (x + 1.0) - lambda;
         Res = Math.exp (y);
      }

      return Res;
   }


   /**
    * Computes and returns the value of the Poisson
    *   distribution function <SPAN CLASS="MATH"><I>F</I>(<I>x</I>)</SPAN> for <SPAN CLASS="MATH"><I>&#955;</I> =</SPAN> <TT>lambda</TT>.
    * 
    */
   public static double cdf (double lambda, int x) {
   /*
    * On our machine, computing a value using gamma is faster than the
    * naive computation for lambdalim > 200.0, slower for lambdalim < 200.0
    */
      final double LAMBDALIM = 200.0;
      int i;
      double term, sum;

      if (lambda < 0.0)
        throw new IllegalArgumentException ("lambda < 0");
      if (lambda == 0.0)
         return 1.0;
      if (x < 0)
         return 0.0;

      if (lambda >= 100.0) {
         if ((double) x >= 10.0*lambda)
            return 1.0;
      } else {
         if ((double) x >= 100.0*Math.max(1.0, lambda))
            return 1.0;
      }   

      /* If lambda > LAMBDALIM, use the Chi2 distribution according to the
         exact relation, with 2x + 2 degrees of freedom

         poisson (lambda, x) = 1 - chiSquare (2x + 2, 2*lambda)

         which also equals 1 - gamma (x + 1, lambda) */
      if (lambda > LAMBDALIM)
         return GammaDist.barF (x + 1.0, 15, lambda);
         
      if(x >= lambda)
         return 1 - PoissonDist.barF(lambda, x+1);

      // Naive computation: sum all prob. from i = x
      sum = term = PoissonDist.prob(lambda, x);
      i = x;
      while(term > EPSILON && i > 0) {
         term *= i/lambda;
         i--;
      }
      sum = term;
      for(int j = i+1; j <= x; j++) {
         term *= lambda/j;
         sum += term;
      }
      return sum;
   }


   /**
    * Computes and returns the value of the complementary Poisson
    *    distribution function, for <SPAN CLASS="MATH"><I>&#955;</I> =</SPAN> <TT>lambda</TT>.
    * <SPAN  CLASS="textit">WARNING:</SPAN> The complementary distribution function is defined as 
    *     
    * <SPAN CLASS="MATH">bar(F)(<I>x</I>) = <I>P</I>[<I>X</I>&nbsp;&gt;=&nbsp;<I>x</I>]</SPAN>.
    * 
    */
   public static double barF (double lambda, int x) {
      final double LAMBDALIM = 200.0;
      int i;
      double term, sum;

      if (lambda < 0)
         throw new IllegalArgumentException ("lambda < 0");
      if (x <= 0)
         return 1.0;

      if (lambda >= 100.0) {
         if ((double) x >= 10.0*lambda)
            return 0.0;
      } else {
         if ((double) x >= 100 + 100.0*Math.max(1.0, lambda))
            return 0.0;
      }   

      /* If lambda > LAMBDALIM, we use the Chi2 distribution according to the
         exact relation, with 2x + 2 degrees of freedom

         cdf (lambda, x) = 1 - ChiSquare.cdf (2x + 2, 2*lambda)

         which also equals   1 - GammaDist.cdf (x + 1, lambda) */

      if (lambda > LAMBDALIM)
         return GammaDist.cdf ((double)x, 15, lambda);

      if (x <= lambda)
         return 1.0 - PoissonDist.cdf(lambda, x - 1);

      // Naive computation: sum all prob. from i = x to i = oo
      final int IMAX = 20;
      // Sum at least IMAX prob. terms from i = s to i = oo
      sum = term = PoissonDist.prob(lambda, x);
      i = x + 1;
      while (term > EPSILON || i <= x + IMAX) {
         term *= lambda/i;
         sum += term;
         i++;
      }
      return sum;
   } 


   /**
    * Performs a linear search to get the inverse function without
    *    precomputed tables.
    * 
    */
   public static int inverseF (double lambda, double u) {
      final double LAMBDALIM = 700.0;
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u is not in range [0,1]");
      if (lambda < 0.0)
        throw new IllegalArgumentException ("lambda < 0");
      if (u >= 1.0)
         return Integer.MAX_VALUE;
      if (u <= 0.0)
         return 0;
      int i = 0;
      
      if (lambda < LAMBDALIM) { // Version la plus rapide
         double sumprev = -1.0;
         double sum, term;
         u *= Math.exp (lambda);
         sum = term = 1.0;
         while (sum < u && sum > sumprev) {
            i++;
            term *= lambda/i;
            sumprev = sum;
            sum += term;
         }
      }
      else { // version legerement plus lente mais ne cause pas d'overflow
         double sum, term;
         term = PoissonDist.prob(lambda, (int)lambda);
         i = (int)lambda;
         if (u <= 0.5) {
            while (term > EPSILON && i > 0) {
               term *= i/lambda;
               i--;
            }
            sum = term;
            while (sum < u) {
               i++;
               term *= lambda/i;
               sum += term;
            }
         }
         else {
            sum = PoissonDist.cdf(lambda, (int)lambda);
            while (term > EPSILON && sum < u) {
               i++;
               term *= lambda/i;
               sum += term;
            }
         }
      }
      return i;
   }

   
   /**
    * Estimates the parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN> of the Poisson distribution
    *    using the maximum likelihood method, from the <SPAN CLASS="MATH"><I>n</I></SPAN> observations 
    *    <SPAN CLASS="MATH"><I>x</I>[<I>i</I>]</SPAN>, 
    * <SPAN CLASS="MATH"><I>i</I> = 0, 1,&#8230;, <I>n</I> - 1</SPAN>.
    * The maximum likelihood estimator 
    * <SPAN CLASS="MATH">hat(&lambda;)</SPAN> satisfy the equation
    *    
    * <SPAN CLASS="MATH">hat(&lambda;) = bar(x)<SUB>n</SUB></SPAN>,
    *    where  <SPAN CLASS="MATH">bar(x)<SUB>n</SUB></SPAN> is the average of 
    * <SPAN CLASS="MATH"><I>x</I>[0],&#8230;, <I>x</I>[<I>n</I> - 1]</SPAN>
    *    (see).
    * 
    * @param x the list of observations used to evaluate parameters
    * 
    *    @param n the number of observations used to evaluate parameters
    * 
    *    @return returns the parameter [
    * <SPAN CLASS="MATH">hat(&lambda;)</SPAN>]
    * 
    */
   public static double[] getMLE (int[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[];
      parameters = new double[1];
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         sum += x[i];
      }

      parameters[0] = (double) sum / (double) n;
      return parameters;
   }


   /**
    * Creates a new instance of a Poisson distribution with parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN>
    *    estimated using the maximum likelihood method based on the <SPAN CLASS="MATH"><I>n</I></SPAN> observations
    *    <SPAN CLASS="MATH"><I>x</I>[<I>i</I>]</SPAN>, 
    * <SPAN CLASS="MATH"><I>i</I> = 0, 1,&#8230;, <I>n</I> - 1</SPAN>.
    * 
    * @param x the list of observations to use to evaluate parameters
    * 
    *    @param n the number of observations to use to evaluate parameters
    * 
    * 
    */
   public static PoissonDist getInstanceFromMLE (int[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new PoissonDist (parameters[0]);
   }


   /**
    * Computes and returns the mean 
    * <SPAN CLASS="MATH"><I>E</I>[<I>X</I>] = <I>&#955;</I></SPAN> of the
    *    Poisson distribution with parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN>.
    * 
    * @return the mean of the Poisson distribution 
    * <SPAN CLASS="MATH"><I>E</I>[<I>X</I>] = <I>&#955;</I></SPAN>
    * 
    */
   public static double getMean (double lambda) {
      if (lambda < 0.0)
       throw new IllegalArgumentException ("lambda < 0");

      return lambda;
   }


   /**
    * Computes and returns the variance <SPAN CLASS="MATH">= <I>&#955;</I></SPAN>
    *    of the Poisson distribution with parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN>.
    * 
    * @return the variance of the Poisson distribution
    *     
    * <SPAN CLASS="MATH">Var[<I>X</I>] = <I>&#955;</I></SPAN>
    * 
    */
   public static double getVariance (double lambda) {
      if (lambda < 0.0)
       throw new IllegalArgumentException ("lambda < 0");

      return lambda;
   }


   /**
    * Computes and returns the standard deviation of the
    *    Poisson distribution with parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN>.
    * 
    * @return the standard deviation of the Poisson distribution
    * 
    */
   public static double getStandardDeviation (double lambda) {
      if (lambda < 0.0)
       throw new IllegalArgumentException ("lambda < 0");

      return Math.sqrt (lambda);
   }


   /**
    * Returns the <SPAN CLASS="MATH"><I>&#955;</I></SPAN> associated with this object.
    * 
    */
   public double getLambda() {
      return lambda;
   }



   /**
    * Sets the <SPAN CLASS="MATH"><I>&#955;</I></SPAN> associated with this object.
    * 
    */
   public void setLambda (double lambda) {
      supportA = 0;
      double epsilon;
      int i, mid, Nmax;
      int imin, imax;
      double sum;
      double[] P;    // Poisson probability terms
      double[] F;    // Poisson cumulative probabilities

      if (lambda < 0.0)
       throw new IllegalArgumentException ("lambda < 0");
      this.lambda = lambda;

      // For lambda > MAXLAMBDAPOISSON, we do not use pre-computed arrays
      if (lambda > MAXLAMBDA) {
         pdf = null;
         cdf = null;
         return;
      }

      // In theory, the Poisson distribution has an infinite range. But
      // for i > Nmax, probabilities should be extremely small.
      Nmax = (int)(lambda + 16*(2 + Math.sqrt (lambda)));
      P = new double[1 + Nmax];
      F = new double[1 + Nmax];

      mid = (int)lambda;
      epsilon = EPSILON * EPS_EXTRA/prob (lambda, mid);
      // For large lambda, mass will lose a few digits of precision
      // We shall normalize by explicitly summing all terms >= epsilon
      sum = P[mid] = 1.0;

      // Start from the maximum and compute terms > epsilon on each side.
      i = mid;
      while (i > 0 && P[i] > epsilon) {
         P[i - 1] = P[i]*i/lambda;
         i--;
         sum += P[i];
      }
      xmin = imin = i;

      i = mid;
      while (P[i] > epsilon) {
         P[i + 1] = P[i]*lambda/(i + 1);
         i++;
         sum += P[i];
         if (i >= Nmax - 1) {
            Nmax *= 2;
            double[] nT = new double[1 + Nmax];
            System.arraycopy (P, 0, nT, 0, P.length);
            P = nT;
            nT = new double[1 + Nmax];
            System.arraycopy (F, 0, nT, 0, F.length);
            F = nT;
         }
      }
      xmax = imax = i;

      // Renormalize the sum of probabilities to 1
      for (i = imin; i <= imax; i++)
         P[i] /= sum;

      // Compute the cumulative probabilities until F >= 0.5, and keep them in
      // the lower part of array, i.e. F[s] contains all P[i] for i <= s 
      F[imin] = P[imin];
      i = imin;
      while (i < imax && F[i] < 0.5) {
         i++;
         F[i] = P[i] + F[i - 1];
      }
      // This is the boundary between F and 1 - F in the CDF
      xmed = i;

      // Compute the cumulative probabilities of the complementary distribution
      // and keep them in the upper part of the array. i.e. F[s] contains all
      // P[i] for i >= s
      F[imax] = P[imax];
      i = imax - 1;
      do {
         F[i] = P[i] + F[i + 1];
         i--;
      } while (i > xmed);

       /* Reset imin because we lose too much precision for a few terms near
      imin when we stop adding terms < epsilon. */
      i = imin;
      while (i < xmed && F[i] < EPSILON)
         i++;
      xmin = imin = i;

      /* Same thing with imax */
      i = imax;
      while (i > xmed && F[i] < EPSILON)
         i--;
      xmax = imax = i;

      pdf = new double[imax + 1 - imin];
      cdf = new double[imax + 1 - imin];
      System.arraycopy (P, imin, pdf, 0, imax-imin+1);
      System.arraycopy (F, imin, cdf, 0, imax-imin+1);

   }


   /**
    * Return a table containing the parameter of the current distribution.
    * 
    * 
    */
   public double[] getParams () {
      double[] retour = {lambda};
      return retour;
   }


   public String toString () {
      return getClass().getSimpleName() + ": lambda = " + lambda;
   }

}
