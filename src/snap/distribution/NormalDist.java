
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
 * Extends the class {@link ContinuousDistribution} for the <EM>normal</EM>
 * distribution (e.g.,). It has mean <SPAN CLASS="MATH"><I>&#956;</I></SPAN>
 * and variance <SPAN CLASS="MATH"><I>&#963;</I><SUP>2</SUP></SPAN>.  Its density function is
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay">
 * <I>f</I> (<I>x</I>) = <I>e</I><SUP>-(x-<I>&#956;</I>)<SUP>2</SUP>/(2<I>&#963;</I><SUP>2</SUP>)</SUP>/((2&pi;)<SUP>1/2</SUP><I>&#963;</I>)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;for  - &#8734; &lt; <I>x</I> &lt; &#8734;,
 * </DIV><P></P>
 * where 
 * <SPAN CLASS="MATH"><I>&#963;</I> &gt; 0</SPAN>.
 * When <SPAN CLASS="MATH"><I>&#956;</I> = 0</SPAN> and <SPAN CLASS="MATH"><I>&#963;</I> = 1</SPAN>, we have the <EM>standard normal</EM>
 * distribution, with corresponding distribution function
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay">
 * <I>F</I>(<I>x</I>) = <I>&#934;</I>(<I>x</I>) = &int;<SUB>-&#8734;</SUB><SUP>x</SUP><I>e</I><SUP>-t<SUP>2</SUP>/2</SUP>&nbsp;<I>dt</I>/(2&pi;)<SUP>1/2</SUP>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;for  - &#8734; &lt; <I>x</I> &lt; &#8734;.
 * </DIV><P></P>
 * The non-static methods <TT>cdf</TT>, <TT>barF</TT>, and <TT>inverseF</TT> are
 * implemented via {@link #cdf01 cdf01}, {@link #barF01 barF01},
 * and {@link #inverseF01 inverseF01}, respectively.
 * 
 */
public class NormalDist extends ContinuousDistribution {
   protected double mu;
   protected double sigma;
   protected static final double RAC2PI = 2.50662827463100050; // Sqrt(2*Pi)

   private static final double[]  AbarF = {
         6.10143081923200418E-1,
         -4.34841272712577472E-1,
         1.76351193643605501E-1,
         -6.07107956092494149E-2,
         1.77120689956941145E-2,
         -4.32111938556729382E-3,
         8.54216676887098679E-4,
         -1.27155090609162743E-4,
         1.12481672436711895E-5,
         3.13063885421820973E-7,
         -2.70988068537762022E-7,
         3.07376227014076884E-8,
         2.51562038481762294E-9,
         -1.02892992132031913E-9,
         2.99440521199499394E-11,
         2.60517896872669363E-11,
         -2.63483992417196939E-12,
         -6.43404509890636443E-13,
         1.12457401801663447E-13,
         1.7281533389986098E-14,
         -4.2641016949424E-15,
         -5.4537197788E-16,
         1.5869760776E-16,
         2.08998378E-17,
         -0.5900E-17
      };




   /**
    * Constructs a <TT>NormalDist</TT> object with default parameters <SPAN CLASS="MATH"><I>&#956;</I> = 0</SPAN>
    *  and 
    * <SPAN CLASS="MATH"><I>&#963;</I> = 1</SPAN>.
    * 
    */
   public NormalDist() {
      setParams (0.0, 1.0);
   }


   /**
    * Constructs a <TT>NormalDist</TT> object with mean <SPAN CLASS="MATH"><I>&#956;</I></SPAN> = <TT>mu</TT>
    *  and standard deviation <SPAN CLASS="MATH"><I>&#963;</I></SPAN> = <TT>sigma</TT>.
    * 
    */
   public NormalDist (double mu, double sigma) {
      setParams (mu, sigma);
   }


   public double density (double x) {
      double z = (x - mu)/sigma;
      return Math.exp (-0.5*z*z)/ (RAC2PI*sigma);
   }

   public double cdf (double x) {
      return cdf01 ((x-mu)/sigma);
   }

   public double barF (double x) {
      return barF01 ((x-mu)/sigma);
   }

   public double inverseF (double u) {
      return mu + sigma * inverseF01 (u);
   }

   public double getMean() {
      return NormalDist.getMean (mu, sigma);
   }

   public double getVariance() {
      return NormalDist.getVariance (mu, sigma);
   }

   public double getStandardDeviation() {
      return NormalDist.getStandardDeviation (mu, sigma);
   }

   /**
    * Same as {@link #density(double,double,double) density}&nbsp;<TT>(0, 1, x)</TT>.
    * 
    */
   public static double density01 (double x)  {
      return Math.exp(-0.5*x*x)/RAC2PI;
   }


   /**
    * Computes the normal density function.
    * 
    */
   public static double density (double mu, double sigma, double x)  {
      if (sigma <= 0)
         throw new IllegalArgumentException ("sigma <= 0");
      double z = (x - mu)/sigma;
      return Math.exp (-0.5*z*z)/ (RAC2PI*sigma);
   }

  /*
    * The precision of double is 16 decimals; we shall thus use coeffmax = 24
    * coefficients. But the approximation is good to 30 decimals of precision
    * with 44 coefficients.
    */
   private static final int COEFFMAX = 24;

   private static final double NORMAL2_A[] = {
    6.10143081923200417926465815756e-1,
    -4.34841272712577471828182820888e-1,
    1.76351193643605501125840298123e-1,
    -6.0710795609249414860051215825e-2,
    1.7712068995694114486147141191e-2,
    -4.321119385567293818599864968e-3,
    8.54216676887098678819832055e-4,
    -1.27155090609162742628893940e-4,
    1.1248167243671189468847072e-5,
    3.13063885421820972630152e-7,
    -2.70988068537762022009086e-7,
    3.0737622701407688440959e-8,
    2.515620384817622937314e-9,
    -1.028929921320319127590e-9,
    2.9944052119949939363e-11,
    2.6051789687266936290e-11,
    -2.634839924171969386e-12,
    -6.43404509890636443e-13,
    1.12457401801663447e-13,
    1.7281533389986098e-14,
    -4.264101694942375e-15,
    -5.45371977880191e-16,
    1.58697607761671e-16,
    2.0899837844334e-17,
    -5.900526869409e-18,
   -9.41893387554e-19
/*,     2.14977356470e-19,
    4.6660985008e-20,
    -7.243011862e-21,
    -2.387966824e-21,
    1.91177535e-22,
    1.20482568e-22,
    -6.72377e-25,
    -5.747997e-24,
    -4.28493e-25,
    2.44856e-25,
    4.3793e-26,
    -8.151e-27,
    -3.089e-27,
    9.3e-29,
    1.74e-28,
    1.6e-29,
    -8.0e-30,
    -2.0e-30
*/
   };


   /**
    * Same as {@link #cdf(double,double,double) cdf}&nbsp;<TT>(0, 1, x)</TT>.
    * 
    */
   public static double cdf01 (double x) {
   /*
    * Returns P[X < x] for the normal distribution.
    * As in J. L. Schonfelder, Math. of Computation, Vol. 32,
    * pp 1232--1240, (1978).
    */

      double t, r;

      if (x <= -XBIG)
         return 0.0;
      if (x >= XBIG)
         return 1.0;

      x = -x/Num.RAC2;
      if (x < 0) {
         x = -x;
         t = (x - 3.75) / (x + 3.75);
         r = 1.0 - 0.5 * Math.exp ( -x * x) * Num.evalCheby (NORMAL2_A, COEFFMAX, t);
      } else {
         t = (x - 3.75) / (x + 3.75);
         r = 0.5 * Math.exp ( -x * x) * Num.evalCheby (NORMAL2_A, COEFFMAX, t);
      }
      return r;
   }


   /**
    * Computes the normal distribution function with mean
    *    <SPAN CLASS="MATH"><I>&#956;</I></SPAN> and variance <SPAN CLASS="MATH"><I>&#963;</I><SUP>2</SUP></SPAN>.
    *    Uses the Chebyshev approximation ,
    * which gives 16 decimals of precision.
    * 
    */
   public static double cdf (double mu, double sigma, double x) {
      if (sigma <= 0)
         throw new IllegalArgumentException ("sigma <= 0");
      return cdf01 ((x-mu)/sigma);
   }


   /**
    * Same as
    *  {@link #barF(double,double,double) barF}&nbsp;<TT>(0, 1, x)</TT>.
    * 
    */
   public static double barF01 (double x) {
   /*
    * Returns P[X >= x] = 1 - F (x) where F is the normal distribution by
    * computing the complementary distribution directly; it is thus more
    * precise in the tail.
    */

      final double KK = 5.30330085889910643300;      // 3.75 Sqrt (2)
      double y, t;
      int neg;

      if (x >= XBIG)
         return 0.0;
      if (x <= -XBIG)
         return 1.0;

      if (x >= 0.0)
         neg = 0;
      else {
         neg = 1;
         x = -x;
      }

      t = (x - KK)/(x + KK);
      y = Num.evalCheby (AbarF, 24, t);
      y = y * Math.exp ( -x * x / 2.0) / 2.0;

      if (neg == 1)
         return 1.0 - y;
      else
         return y;
   }


   /**
    * Computes the complementary normal distribution function
    *   
    * <SPAN CLASS="MATH">bar(F)(<I>x</I>) = 1 - <I>&#934;</I>((<I>x</I> - <I>&#956;</I>)/<I>&#963;</I>)</SPAN>,
    *   with  mean <SPAN CLASS="MATH"><I>&#956;</I></SPAN> and variance <SPAN CLASS="MATH"><I>&#963;</I><SUP>2</SUP></SPAN>.
    *   Uses a Chebyshev series giving 16 decimal digits of
    *   precision.
    * 
    */
   public static double barF (double mu, double sigma, double x) {
      if (sigma <= 0)
         throw new IllegalArgumentException ("sigma <= 0");
      return barF01 ((x-mu)/sigma);
   }

    private static final double[] InvP1 = {
        0.160304955844066229311E2,
       -0.90784959262960326650E2,
        0.18644914861620987391E3,
       -0.16900142734642382420E3,
        0.6545466284794487048E2,
       -0.864213011587247794E1,
        0.1760587821390590
    };

    private static final double[] InvQ1 = {
        0.147806470715138316110E2,
       -0.91374167024260313396E2,
        0.21015790486205317714E3,
       -0.22210254121855132366E3,
        0.10760453916055123830E3,
       -0.206010730328265443E2,
        0.1E1
    };

    private static final double[] InvP2 = {
       -0.152389263440726128E-1,
        0.3444556924136125216,
       -0.29344398672542478687E1,
        0.11763505705217827302E2,
       -0.22655292823101104193E2,
        0.19121334396580330163E2,
       -0.5478927619598318769E1,
        0.237516689024448000
    };

    private static final double[] InvQ2 = {
      -0.108465169602059954E-1,
       0.2610628885843078511,
      -0.24068318104393757995E1,
       0.10695129973387014469E2,
      -0.23716715521596581025E2,
       0.24640158943917284883E2,
      -0.10014376349783070835E2,
       0.1E1
    };

    private static final double[] InvP3 = {
        0.56451977709864482298E-4,
        0.53504147487893013765E-2,
        0.12969550099727352403,
        0.10426158549298266122E1,
        0.28302677901754489974E1,
        0.26255672879448072726E1,
        0.20789742630174917228E1,
        0.72718806231556811306,
        0.66816807711804989575E-1,
       -0.17791004575111759979E-1,
        0.22419563223346345828E-2
    };

    private static final double[] InvQ3 = {
        0.56451699862760651514E-4,
        0.53505587067930653953E-2,
        0.12986615416911646934,
        0.10542932232626491195E1,
        0.30379331173522206237E1,
        0.37631168536405028901E1,
        0.38782858277042011263E1,
        0.20372431817412177929E1,
        0.1E1
    };

   /**
    * Same as
    * {@link #inverseF(double,double,double) inverseF}&nbsp;<TT>(0, 1, u)</TT>.
    * 
    */
   public static double inverseF01 (double u) {
       /*
        * Returns the inverse of the cdf of the normal distribution.
        * Rational approximations giving 16 decimals of precision.
        * J.M. Blair, C.A. Edwards, J.H. Johnson, "Rational Chebyshev
        * approximations for the Inverse of the Error Function", in
        * Mathematics of Computation, Vol. 30, 136, pp 827, (1976)
        */

      int i;
      boolean negatif;
      double y, z, v, w;
      double x = u;

      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u is not in [0, 1]");
      if (u <= 0.0)
         return Double.NEGATIVE_INFINITY;
      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;

      // Transform x as argument of InvErf
      x = 2.0 * x - 1.0;
      if (x < 0.0) {
         x = -x;
         negatif = true;
      } else {
         negatif = false;
      }

      if (x <= 0.75) {
         y = x * x - 0.5625;
         v = w = 0.0;
         for (i = 6; i >= 0; i--) {
            v = v * y + InvP1[i];
            w = w * y + InvQ1[i];
         }
         z = (v / w) * x;

      } else if (x <= 0.9375) {
         y = x * x - 0.87890625;
         v = w = 0.0;
         for (i = 7; i >= 0; i--) {
            v = v * y + InvP2[i];
            w = w * y + InvQ2[i];
         }
         z = (v / w) * x;

      } else {
         if (u > 0.5)
            y = 1.0 / Math.sqrt (-Math.log (1.0 - x));
         else
            y = 1.0 / Math.sqrt (-Math.log (2.0 * u));
         v = 0.0;
         for (i = 10; i >= 0; i--)
            v = v * y + InvP3[i];
         w = 0.0;
         for (i = 8; i >= 0; i--)
            w = w * y + InvQ3[i];
         z = (v / w) / y;
      }

      if (negatif) {
         if (u < 1.0e-105) {
            final double RACPI = 1.77245385090551602729;
            w = Math.exp (-z * z) / RACPI;  // pdf
            y = 2.0 * z * z;
            v = 1.0;
            double term = 1.0;
            // Asymptotic series for erfc(z) (apart from exp factor)
            for (i = 0; i < 6; ++i) {
               term *= -(2 * i + 1) / y;
               v += term;
            }
            // Apply 1 iteration of Newton solver to get last few decimals
            z -= u / w - 0.5 * v / z;

         }
         return -(z * Num.RAC2);

      } else
         return z * Num.RAC2;
   }


   /**
    * Computes the inverse normal distribution function
    *   with mean <SPAN CLASS="MATH"><I>&#956;</I></SPAN> and variance <SPAN CLASS="MATH"><I>&#963;</I><SUP>2</SUP></SPAN>. Uses different
    *   rational Chebyshev approximations.
    *   Returns 16 decimal digits of precision for 
    * <SPAN CLASS="MATH">2.2&#215;10<SUP>-308</SUP> &lt; <I>u</I> &lt; 1</SPAN>.
    * 
    */
   public static double inverseF (double mu, double sigma, double u)  {
      if (sigma <= 0)
         throw new IllegalArgumentException ("sigma <= 0");
      return mu + sigma*inverseF01 (u);
   }


   /**
    * Estimates the parameters 
    * <SPAN CLASS="MATH">(<I>&#956;</I>, <I>&#963;</I>)</SPAN> of the normal distribution
    *    using the maximum likelihood method, from the <SPAN CLASS="MATH"><I>n</I></SPAN> observations
    *    <SPAN CLASS="MATH"><I>x</I>[<I>i</I>]</SPAN>, 
    * <SPAN CLASS="MATH"><I>i</I> = 0, 1,&#8230;, <I>n</I> - 1</SPAN>. The estimates are returned in a two-element
    *     array, in regular order: [<SPAN CLASS="MATH">hat(&mu;)</SPAN>, 
    * <SPAN CLASS="MATH">hat(&sigma;)</SPAN>].
    * 
    * @param x the list of observations used to evaluate parameters
    * 
    *    @param n the number of observations used to evaluate parameters
    * 
    *    @return returns the parameters [<SPAN CLASS="MATH">hat(&mu;)</SPAN>, 
    * <SPAN CLASS="MATH">hat(&sigma;)</SPAN>]
    * 
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double[] parameters = new double[2];
      double sum = 0.0;
      for (int i = 0; i < n; i++)
         sum += x[i];
      parameters[0] = sum / n;

      sum = 0.0;
      for (int i = 0; i < n; i++)
         sum = sum + (x[i] - parameters[0]) * (x[i] - parameters[0]);
      parameters[1] = Math.sqrt (sum / n);

      return parameters;
   }


   /**
    * Creates a new instance of a normal distribution with parameters <SPAN CLASS="MATH"><I>&#956;</I></SPAN> and <SPAN CLASS="MATH"><I>&#963;</I></SPAN>
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
   public static NormalDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new NormalDist (parameters[0], parameters[1]);
   }


   /**
    * Computes and returns the mean 
    * <SPAN CLASS="MATH"><I>E</I>[<I>X</I>] = <I>&#956;</I></SPAN> of the normal distribution
    *    with parameters <SPAN CLASS="MATH"><I>&#956;</I></SPAN> and <SPAN CLASS="MATH"><I>&#963;</I></SPAN>.
    * 
    * @return the mean of the normal distribution 
    * <SPAN CLASS="MATH"><I>E</I>[<I>X</I>] = <I>&#956;</I></SPAN>
    * 
    */
   public static double getMean (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");

      return mu;
   }


   /**
    * Computes and returns the variance 
    * <SPAN CLASS="MATH">Var[<I>X</I>] = <I>&#963;</I><SUP>2</SUP></SPAN> of the
    *    normal distribution with parameters <SPAN CLASS="MATH"><I>&#956;</I></SPAN> and <SPAN CLASS="MATH"><I>&#963;</I></SPAN>.
    * 
    * @return the variance of the normal distribution 
    * <SPAN CLASS="MATH">Var[<I>X</I>] = <I>&#963;</I><SUP>2</SUP></SPAN>
    * 
    */
   public static double getVariance (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");

      return (sigma * sigma);
   }


   /**
    * Computes and returns the standard deviation <SPAN CLASS="MATH"><I>&#963;</I></SPAN> of the
    *    normal distribution with parameters <SPAN CLASS="MATH"><I>&#956;</I></SPAN> and <SPAN CLASS="MATH"><I>&#963;</I></SPAN>.
    * 
    * @return the standard deviation of the normal distribution
    * 
    */
   public static double getStandardDeviation (double mu, double sigma) {
      return sigma;
   }


   /**
    * Returns the parameter <SPAN CLASS="MATH"><I>&#956;</I></SPAN>.
    * 
    */
   public double getMu() {
      return mu;
   }


   /**
    * Returns the parameter <SPAN CLASS="MATH"><I>&#963;</I></SPAN>.
    * 
    */
   public double getSigma() {
      return sigma;
   }


   /**
    * Sets the parameters <SPAN CLASS="MATH"><I>&#956;</I></SPAN> and <SPAN CLASS="MATH"><I>&#963;</I></SPAN> of this object.
    * 
    */
   public void setParams (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      this.mu = mu;
      this.sigma = sigma;
   }


   /**
    * Return a table containing the parameters of the current distribution.
    *    This table is put in regular order: [<SPAN CLASS="MATH"><I>&#956;</I></SPAN>, <SPAN CLASS="MATH"><I>&#963;</I></SPAN>].
    * 
    * 
    */
   public double[] getParams () {
      double[] retour = {mu, sigma};
      return retour;
   }


   public String toString () {
      return getClass().getSimpleName() + " : mu = " + mu + ", sigma = " + sigma;
   }

}
