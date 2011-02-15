
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
 * Extends the class {@link ContinuousDistribution} for
 * the <EM>gamma</EM> distribution with
 * shape parameter 
 * <SPAN CLASS="MATH"><I>&#945;</I> &gt; 0</SPAN> and scale parameter
 * 
 * <SPAN CLASS="MATH"><I>&#955;</I> &gt; 0</SPAN>.
 * The density is
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay">
 * <I>f</I> (<I>x</I>) = <I>&#955;</I><SUP><I>&#945;</I></SUP><I>x</I><SUP><I>&#945;</I>-1</SUP><I>e</I><SUP>-<I>&#955;</I>x</SUP>/<I>&#915;</I>(<I>&#945;</I>),&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;for <I>x</I> &gt; 0,
 * </DIV><P></P>
 * where <SPAN CLASS="MATH"><I>&#915;</I></SPAN> is the gamma function, defined by
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay">
 * <I>&#915;</I>(<I>&#945;</I>) = &int;<SUB>0</SUB><SUP>&#8734;</SUP><I>x</I><SUP><I>&#945;</I>-1</SUP><I>e</I><SUP>-x</SUP><I>dx</I>.
 * </DIV><P></P>
 * In particular, 
 * <SPAN CLASS="MATH"><I>&#915;</I>(<I>n</I>) = (<I>n</I> - 1)!</SPAN> when <SPAN CLASS="MATH"><I>n</I></SPAN> is a positive integer.
 * 
 */
public class GammaDist extends ContinuousDistribution {
   private double alpha;
   private double lambda;
   private double logFactor;      // Log (lambda^alpha / Gamma (alpha))
   private static final double ALIM = 1.0E5;

   private static class Function implements MathFunction {
      // For MLE
      private int n;
      private double empiricalMean;
      private double sumLn;

      public Function (int n, double empiricalMean, double sumLn) {
         this.n = n;
         this.empiricalMean = empiricalMean;
         this.sumLn = sumLn;
      }

      public double evaluate (double x) {
         if (x <= 0.0) return 1.0e200;
         return (n * Math.log (empiricalMean / x) + n * Num.digamma (x) - sumLn);
      }
   }


   private static class myFunc implements MathFunction {
      // For inverseF
      protected int d;
      protected double alp, u;

      public myFunc (double alp, int d, double u) {
         this.alp = alp;
         this.d = d;
         this.u = u;
      }

      public double evaluate (double x) {
         return u - GammaDist.cdf(alp, d, x);
      }
   }




   /**
    * Constructs a <TT>GammaDist</TT> object with parameters
    *   <SPAN CLASS="MATH"><I>&#945;</I></SPAN> = <TT>alpha</TT> and <SPAN CLASS="MATH"><I>&#955;</I> = 1</SPAN>.
    * 
    */
   public GammaDist (double alpha) {
      setParams (alpha, 1.0, decPrec);
   }


   /**
    * Constructs a <TT>GammaDist</TT> object with parameters
    *   <SPAN CLASS="MATH"><I>&#945;</I></SPAN> = <TT>alpha</TT>  and <SPAN CLASS="MATH"><I>&#955;</I></SPAN> = <TT>lambda</TT>.
    * 
    */
   public GammaDist (double alpha, double lambda) {
      setParams (alpha, lambda, decPrec);
   }


   /**
    * Constructs a <TT>GammaDist</TT> object with parameters <SPAN CLASS="MATH"><I>&#945;</I></SPAN> =
    *    <TT>alpha</TT> and  <SPAN CLASS="MATH"><I>&#955;</I></SPAN> = <TT>lambda</TT>, and approximations of
    *    roughly <TT>d</TT> decimal digits of precision when computing functions.
    * 
    */
   public GammaDist (double alpha, double lambda, int d) {
      setParams (alpha, lambda, d);
   }


   private static double mybelog (double x)
   /*
    * This is the function  1 + (1 - x*x + 2*x*log(x)) / ((1 - x)*(1 - x))
    */
   {
      if (x < 1.0e-30)
         return 0.0;
      if (x > 1.0e30)
         return 2.0*(Math.log(x) - 1.0) / x;
      if (x == 1.0)
         return 1.0;

      double t = 1.0 - x;
      if (x < 0.9 || x > 1.1) {
         double w = (t + x*Math.log(x)) / (t*t);
         return 2.0 * w;
      }

      // For x near 1, use a series expansion to avoid loss of precision.
      double term;
      final double EPS = 1.0e-12;
      double tpow = 1.0;
      double sum = 0.5;
      int j = 3;
      do {
         tpow *= t;
         term = tpow / (j * (j - 1));
         sum += term;
         j++;
      } while (Math.abs (term / sum) > EPS);
      return 2.0*sum;
   }


   public double density (double x) {
      if (x <= 0)
         return 0.0;
      double z = logFactor + (alpha - 1.0) * Math.log(x) - lambda * x;
      if (z > -XBIGM)
         return Math.exp (z);
      else
         return 0.0;
   }

   public double cdf (double x) {
      return cdf (alpha, lambda, decPrec, x);
   }

   public double barF (double x) {
      return barF (alpha, lambda, decPrec, x);
   }

   public double inverseF (double u) {
      return inverseF (alpha, decPrec, u)/lambda;
   }

   public double getMean() {
      return GammaDist.getMean (alpha, lambda);
   }

   public double getVariance() {
      return GammaDist.getVariance (alpha, lambda);
   }

   public double getStandardDeviation() {
      return GammaDist.getStandardDeviation (alpha, lambda);
   }


   /**
    * Computes the density function at <SPAN CLASS="MATH"><I>x</I></SPAN>.
    * 
    */
   public static double density (double alpha, double lambda, double x) {
      if (alpha <= 0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");

      if (x <= 0)
         return 0.0;
      double z = alpha * Math.log (lambda*x) - lambda*x - Num.lnGamma (alpha);
      if (z > -XBIGM)
         return Math.exp (z) / x;
      else
         return 0.0;
   }


   /**
    * Returns an approximation of the gamma distribution
    * function with  parameters <SPAN CLASS="MATH"><I>&#945;</I></SPAN> = <TT>alpha</TT> and
    * <SPAN CLASS="MATH"><I>&#955;</I></SPAN> = <TT>lambda</TT>.
    *  The function tries to return <SPAN CLASS="MATH"><I>d</I></SPAN> decimals digits of precision.
    * For <SPAN CLASS="MATH"><I>&#945;</I></SPAN> not too large (e.g., 
    * <SPAN CLASS="MATH"><I>&#945;</I>&nbsp;&lt;=&nbsp;1000</SPAN>),
    *  <SPAN CLASS="MATH"><I>d</I></SPAN> gives a good idea of the precision attained.
    * 
    */
   public static double cdf (double alpha, double lambda, int d, double x) {
      return cdf (alpha, d, lambda*x);
   }


   /**
    * Equivalent to <TT>cdf (alpha, 1.0, d, x)</TT>.
    * 
    */
   public static double cdf (double alpha, int d, double x) {
      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (d <= 0)
        throw new IllegalArgumentException ("d <= 0");
      if (x <= 0.0)
         return 0.0;
      if (1.0 == alpha)
         return ExponentialDist.cdf (1.0, x);

      if (alpha > 10.0) {
         if (x > alpha * 10.0)
            return 1.0;
      } else {
         if (x > XBIG)
            return 1.0;
      }

      if (alpha >= ALIM) {
         double d2 = x + 1.0/3.0 - alpha - 0.02/alpha;
         double S = alpha - 1.0/2.0;
         double w = mybelog(S/x);
         double y = d2 * Math.sqrt(w/x);
         return NormalDist.cdf01 (y);
      }

      if (x <= 1.0 || x < alpha) {
         double factor, z, rn, term;
         factor = Math.exp (alpha*Math.log (x) - x - Num.lnGamma (alpha));
         final double EPS = EPSARRAY[d];
         z = 1.0;
         term = 1.0;
         rn = alpha;
         do {
            rn += 1.0;
            term *= x/rn;
            z += term;
         } while (term >= EPS * z);
         return z*factor/alpha;

      } else
         return 1.0 - barF (alpha, d, x);
   }


   /**
    * Computes the complementary distribution function.
    * 
    */
   public static double barF (double alpha, double lambda, int d, double x) {
      return barF (alpha, d, lambda*x);
   }


   /**
    * Same as {@link #barF(double,double,int,double) barF}&nbsp;<TT>(alpha, 1.0, d, x)</TT>.
    * 
    */
   public static double barF (double alpha, int d, double x) {
      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (d <= 0)
        throw new IllegalArgumentException ("d <= 0");
      if (x <= 0.0)
         return 1.0;
      if (1.0 == alpha)
         return ExponentialDist.barF (1.0, x);

      if (alpha >= 70.0) {
         if (x >= alpha * XBIG)
            return 0.0;
      } else {
         if (x >= XBIGM)
            return 0.0;
      }

      if (alpha >= ALIM) {
         double d2 = x + 1.0/3.0 - alpha - 0.02/alpha;
         double S = alpha - 1.0/2.0;
         double w = mybelog(S/x);
         double y = d2 * Math.sqrt(w/x);
         return NormalDist.barF01 (y);
      }

      if (x <= 1.0 || x < alpha)
         return 1.0 - cdf (alpha, d, x);

      double[] V = new double[6];
      final double EPS = EPSARRAY[d];
      final double RENORM = 1.0E100;
      double R, dif;
      int i;
      double factor = Math.exp (alpha*Math.log (x) - x - Num.lnGamma (alpha));

      double A = 1.0 - alpha;
      double B = A + x + 1.0;
      double term = 0.0;
      V[0] = 1.0;
      V[1] = x;
      V[2] = x + 1.0;
      V[3] = x * B;
      double res = V[2]/V[3];

      do {
         A += 1.0;
         B += 2.0;
         term += 1.0;
         V[4] = B * V[2] - A * term * V[0];
         V[5] = B * V[3] - A * term * V[1];
         if (V[5] != 0.0) {
            R = V[4]/V[5];
            dif = Math.abs (res - R);
            if (dif <= EPS*R)
               return factor*res;
            res = R;
         }
         for (i = 0; i < 4; i++)
            V[i] = V[i + 2];
         if (Math.abs (V[4]) >= RENORM) {
            for (i = 0; i < 4; i++)
               V[i] /= RENORM;
         }
      } while (true);
   }


   /**
    * Computes the inverse distribution function.
    * 
    */
   public static double inverseF (double alpha, double lambda, int d,
                                  double u) {
      return inverseF (alpha, d, u)/lambda;
   }


   /**
    * Same as
    * {@link #inverseF(double,double,int,double) inverseF}&nbsp;<TT>(alpha, 1, d, u)</TT>.
    * 
    */
   public static double inverseF (double alpha, int d, double u) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (u > 1.0 || u < 0.0)
         throw new IllegalArgumentException ("u not in [0,1]");
      if (u <= 0.0)
         return 0;
      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;
      if (d <= 0)
         throw new IllegalArgumentException ("d <= 0");
      if (d > 15)
         d = 15;
      final double EPS = Math.pow (10.0, -d);

      double sigma = GammaDist.getStandardDeviation (alpha, 1.0);
      double x = NormalDist.inverseF (alpha, sigma, u);
      double v = GammaDist.cdf (alpha, d, x);
      double xmax;
      if (alpha < 1.0)
         xmax = 100.0;
      else
         xmax = alpha + 40.0 * sigma;
      myFunc f = new myFunc (alpha, d, u);

     if (u <= 1.0e-8 || alpha <= 1.5) {
         if (v < u)
            return RootFinder.bisection (x, xmax, f, EPS);
         else
            return RootFinder.bisection (0, x, f, EPS);
      } else {
          if (v < u)
            return RootFinder.brentDekker (x, xmax, f, EPS);
         else
            return RootFinder.brentDekker (0, x, f, EPS);
      }
   }


   /**
    * Estimates the parameters 
    * <SPAN CLASS="MATH">(<I>&#945;</I>, <I>&#955;</I>)</SPAN> of the gamma distribution
    *    using the maximum likelihood method, from the <SPAN CLASS="MATH"><I>n</I></SPAN> observations
    *    <SPAN CLASS="MATH"><I>x</I>[<I>i</I>]</SPAN>, 
    * <SPAN CLASS="MATH"><I>i</I> = 0, 1,&#8230;, <I>n</I> - 1</SPAN>. The estimates are returned in a two-element
    *     array, in regular order: [<SPAN CLASS="MATH"><I>&#945;</I></SPAN>, <SPAN CLASS="MATH"><I>&#955;</I></SPAN>].
    * 
    * @param x the list of observations to use to evaluate parameters
    * 
    *    @param n the number of observations to use to evaluate parameters
    * 
    *    @return returns the parameters [
    * <SPAN CLASS="MATH">hat(&alpha;)</SPAN>, 
    * <SPAN CLASS="MATH">hat(&lambda;)</SPAN>]
    * 
    */
   public static double[] getMLE (double[] x, int n) {
      double parameters[];
      double sum = 0.0;
      double sumLn = 0.0;
      double empiricalMean;
      double alphaMME;
      double a;
      final double LN_EPS = Num.LN_DBL_MIN - Num.LN2;

      parameters = new double[2];
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      for (int i = 0; i < n; i++)
      {
         sum += x[i];
         if (x[i] <= 0.0)
            sumLn += LN_EPS;
         else
            sumLn += Math.log (x[i]);
      }
      empiricalMean = sum / (double) n;

      sum = 0.0;
      for (int i = 0; i < n; i++) {
         sum += (x[i] - empiricalMean) * (x[i] - empiricalMean);
      }

      alphaMME = (empiricalMean * empiricalMean * (double) n) / sum;
      if ((a = alphaMME - 10.0) <= 0) {
         a = 1.0e-5;
      }

      Function f = new Function (n, empiricalMean, sumLn);
      parameters[0] = RootFinder.brentDekker (a, alphaMME + 10.0, f, 1e-7);
      parameters[1] = parameters[0] / empiricalMean;

      return parameters;
   }


   /**
    * Creates a new instance of a gamma distribution with parameters <SPAN CLASS="MATH"><I>&#945;</I></SPAN> and <SPAN CLASS="MATH"><I>&#955;</I></SPAN>
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
   public static GammaDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new GammaDist (parameters[0], parameters[1]);
   }


   /**
    * Computes and returns the mean 
    * <SPAN CLASS="MATH"><I>E</I>[<I>X</I>] = <I>&#945;</I>/<I>&#955;</I></SPAN>
    *    of the gamma distribution with parameters <SPAN CLASS="MATH"><I>&#945;</I></SPAN> and <SPAN CLASS="MATH"><I>&#955;</I></SPAN>.
    * 
    * @return the mean of the gamma distribution 
    * <SPAN CLASS="MATH"><I>E</I>[<I>X</I>] = <I>&#945;</I>/<I>&#955;</I></SPAN>
    * 
    */
   public static double getMean (double alpha, double lambda) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (alpha / lambda);
   }


   /**
    * Computes and returns the variance 
    * <SPAN CLASS="MATH">Var[<I>X</I>] = <I>&#945;</I>/<I>&#955;</I><SUP>2</SUP></SPAN>
    *    of the gamma distribution with parameters <SPAN CLASS="MATH"><I>&#945;</I></SPAN> and <SPAN CLASS="MATH"><I>&#955;</I></SPAN>.
    * 
    * @return the variance of the gamma distribution 
    * <SPAN CLASS="MATH">Var[<I>X</I>] = <I>&#945;</I>/<I>&#955;</I><SUP>2</SUP></SPAN>
    * 
    */
   public static double getVariance (double alpha, double lambda) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (alpha / (lambda * lambda));
   }


   /**
    * Computes and returns the standard deviation of the gamma
    *    distribution with parameters <SPAN CLASS="MATH"><I>&#945;</I></SPAN> and <SPAN CLASS="MATH"><I>&#955;</I></SPAN>.
    * 
    * @return the standard deviation of the gamma distribution
    * 
    */
   public static double getStandardDeviation (double alpha, double lambda) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (Math.sqrt(alpha) / lambda);
   }


   /**
    * Return the parameter <SPAN CLASS="MATH"><I>&#945;</I></SPAN> for this object.
    * 
    */
   public double getAlpha() {
      return alpha;
   }


   /**
    * Return the parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN> for this object.
    * 
    */
   public double getLambda() {
      return lambda;
   }


   public void setParams (double alpha, double lambda, int d) {
      if (alpha <= 0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");

      this.alpha   = alpha;
      this.lambda  = lambda;
      this.decPrec = d;
      logFactor    = alpha * Math.log(lambda) - Num.lnGamma (alpha);
      supportA = 0.0;
    } 


   /**
    * Return a table containing the parameters of the current distribution.
    *    This table is put in regular order: [<SPAN CLASS="MATH"><I>&#945;</I></SPAN>, <SPAN CLASS="MATH"><I>&#955;</I></SPAN>].
    * 
    * 
    */
   public double[] getParams () {
      double[] retour = {alpha, lambda};
      return retour;
   }


   public String toString () {
      return getClass().getSimpleName() + " : alpha = " + alpha + ", lambda = " + lambda;
   }

}
