
package snap.distribution;



/**
 * Classes implementing continuous distributions should inherit from this base
 *  class.
 * Such distributions are characterized by a <SPAN  CLASS="textit">density</SPAN> function <SPAN CLASS="MATH"><I>f</I> (<I>x</I>)</SPAN>,
 * thus the signature of a <TT>density</TT> method is supplied here.
 * This class also provides default implementations for <SPAN CLASS="MATH">bar(F)(<I>x</I>)</SPAN>
 * and for <SPAN CLASS="MATH"><I>F</I><SUP>-1</SUP>(<I>u</I>)</SPAN>, the latter using the Brent-Dekker method to find the inverse
 * of a generic distribution function <SPAN CLASS="MATH"><I>F</I></SPAN>. 
 */
public abstract class ContinuousDistribution implements Distribution  {
   @Deprecated
   public int decPrec = 15;

   private int getDecPrec() {
      return decPrec;
   }


    // x infinity for some distributions
//    protected static final double XINF = Double.MAX_VALUE;
    protected static final double XBIG = 100.0;
    protected static final double XBIGM = 1000.0;

    // [supportA, supportB] is the support of the pdf(x)
//    protected double supportA = -XINF;
//    protected double supportB = XINF;
    protected double supportA = Double.NEGATIVE_INFINITY;
    protected double supportB = Double.POSITIVE_INFINITY;

    // EPSARRAY[j]: Epsilon required for j decimal degits of precision
    protected static final double[] EPSARRAY = {
    0.5, 0.5E-1, 0.5E-2, 0.5E-3, 0.5E-4, 0.5E-5, 0.5E-6, 0.5E-7, 0.5E-8,
    0.5E-9, 0.5E-10, 0.5E-11, 0.5E-12, 0.5E-13, 0.5E-14, 0.5E-15, 0.5E-16,
    0.5E-17, 0.5E-18, 0.5E-19, 0.5E-20, 0.5E-21, 0.5E-22, 0.5E-23, 0.5E-24,
    0.5E-25, 0.5E-26, 0.5E-27, 0.5E-28, 0.5E-29, 0.5E-30, 0.5E-31, 0.5E-32,
    0.5E-33, 0.5E-34, 0.5E-35
    };

   /**
    * Returns <SPAN CLASS="MATH"><I>f</I> (<I>x</I>)</SPAN>, the density evaluated at <SPAN CLASS="MATH"><I>x</I></SPAN>.
    * 
    * @param x value at which the density is evaluated
    * 
    *    @return density function evaluated at <TT>x</TT>
    * 
    */
   public abstract double density (double x);


   /**
    * Returns the complementary distribution function.
    *     The default implementation computes 
    * <SPAN CLASS="MATH">bar(F)(<I>x</I>) = 1 - <I>F</I>(<I>x</I>)</SPAN>.
    *     
    * @param x value at which the complementary distribution function is evaluated
    * 
    *    @return complementary distribution function evaluated at <TT>x</TT>
    * 
    */
   public double barF (double x) {
      return 1.0 - cdf (x);
   }


   private void findInterval (double u, double [] iv) {
      // Finds an interval [a, b] that certainly contains x defined as
      // u = cdf(x). The result is written in iv[0] = a and iv[1] = b.

      if (u > 1.0 || u < 0.0)
         throw new IllegalArgumentException ("u not in [0, 1]");
      final double XLIM =  Double.MAX_VALUE/2.0;
      final double B0 = 8.0;
      double b = B0;
      while (b < XLIM && u > cdf(b))
         b *= 2.0;
      if (b > B0) {
         iv[0] = b/2.0;
         iv[1] = Math.min (b, supportB);
         return;
      }

      double a = -B0;
      while (a > -XLIM && u < cdf(a))
         a *= 2.0;
      if (a < -B0) {
         iv[1] = a/2.0;
         iv[0] = Math.max (a, supportA);
         return;
      }
      iv[0] = Math.max (a, supportA);
      iv[1] = Math.min (b, supportB);
   }

   /**
    * Computes the inverse distribution function 
    * <SPAN CLASS="MATH"><I>x</I> = <I>F</I><SUP>-1</SUP>(<I>u</I>)</SPAN>,
    *     using the Brent-Dekker method. The interval <SPAN CLASS="MATH">[<I>a</I>, <I>b</I>]</SPAN> <SPAN  CLASS="textit">must</SPAN> contain
    *     the root <SPAN CLASS="MATH"><I>x</I></SPAN> such that 
    * <SPAN CLASS="MATH"><I>F</I>(<I>a</I>)&nbsp;&lt;=&nbsp;<I>u</I>&nbsp;&lt;=&nbsp;<I>F</I>(<I>b</I>)</SPAN>, where <SPAN CLASS="MATH"><I>u</I> = <I>F</I>(<I>x</I>)</SPAN>.
    *     The calculations are done with an approximate precision of <TT>tol</TT>.
    *     Returns 
    * <SPAN CLASS="MATH"><I>x</I> = <I>F</I><SUP>-1</SUP>(<I>u</I>)</SPAN>. Restrictions: 
    * <SPAN CLASS="MATH"><I>u</I>&#8712;[0, 1]</SPAN>.
    * 
    * @param a left endpoint of initial interval
    * 
    *    @param b right endpoint of initial interval
    * 
    *    @param u value at which the inverse distribution function is evaluated
    * 
    *    @param tol accuracy goal
    * 
    *    @return inverse distribution function evaluated at <TT>u</TT>
    * 
    */
   public double inverseBrent (double a, double b, double u, double tol)  {
      if (u > 1.0 || u < 0.0)
          throw new IllegalArgumentException ("u not in [0, 1]");
      if (b < a) {
         double ctemp = a;   a = b;   b = ctemp;
      }
      if (u <= 0.0) {
          System.out.println ("********** WARNING,  inverseBrent:   u = 0");
          return supportA;
      }
      if (u >= 1.0) {
          System.out.println ("********** WARNING,  inverseBrent:   u = 1");
          return supportB;
      }
      final int MAXITER = 50;      // Maximum number of iterations
      tol += EPSARRAY[decPrec] + Num.DBL_EPSILON;    // in case tol is too small
      double ua = cdf(a) - u;
      if (ua > 0.0)
          throw new IllegalArgumentException ("u < cdf(a)");
      double ub = cdf(b) - u;
      if (ub < 0.0)
          throw new IllegalArgumentException ("u > cdf(b)");

      final boolean DEBUG = false;
      if (DEBUG) {
//          String ls = System.getProperty("line.separator");
//          System.out.println (
//             "-------------------------------------------------------------"
//              + ls + "u = " + PrintfFormat.g (20, 15, u));
//          System.out.println
//           (ls + "iter           b                  c               F(x) - u" + ls);
      }
      // Initialize
      double c = a;
      double uc = ua;
      double len = b - a;
      double t = len;
      if (Math.abs(uc) < Math.abs(ub)) {
            a = b; b = c; c = a;
            ua = ub; ub = uc; uc = ua;
      }
      int i;
      for (i = 0; i < MAXITER; ++i) {
         double tol1 = tol + 4.0*Num.DBL_EPSILON*Math.abs(b);
         double xm = 0.5*(c - b);
         if (DEBUG) {
//            System.out.println (PrintfFormat.d (3, i) + "  " +
//                PrintfFormat.g (18, decPrec, b) + "  " +
//                PrintfFormat.g (18, decPrec, c) + "  " +
//                PrintfFormat.g (14, 4, ub));
         }
         if (Math.abs(ub) == 0.0 || (Math.abs(xm) <= tol1)) {
            if (b <= supportA) return supportA;
            if (b >= supportB) return supportB;
            return b;
         }

         double s, p, q, r;
         if ((Math.abs(t) >= tol1) && (Math.abs(ua) > Math.abs(ub))) {
            if (a == c) {
               // linear interpolation
               s = ub/ua;
               q = 1.0 - s;
               p = 2.0 * xm * s;
            } else {
               // quadratic interpolation
               q = ua / uc;
               r = ub / uc;
               s = ub / ua;
               p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
               q = (q - 1.0)*(r - 1.0)* (s - 1.0);
            }
            if (p > 0.0)
               q = -q;
            p = Math.abs(p);

            // Accept interpolation?
            if ((2.0*p >= (3.0*xm*q - Math.abs(q*tol1))) ||
                (p >= Math.abs(0.5*t*q))) {
               len = xm;
               t = len;
            } else {
               t = len;
               len = p/q;
            }

         } else {
            len = xm;
            t = len;
         }

         a = b;
         ua = ub;
         if (Math.abs(len) > tol1)
            b += len;
         else if (xm < 0.0)
            b -= tol1;
         else
            b += tol1;
         ub = cdf(b) - u;

         if (ub*(uc/Math.abs(uc)) > 0.0) {
            c = a;
            uc = ua;
            len = b - a;
            t = len;
         } else if (Math.abs(uc) < Math.abs(ub)) {
            a = b; b = c; c = a;
            ua = ub; ub = uc; uc = ua;
         }
      }
      if (i >= MAXITER) {
         String lineSep = System.getProperty("line.separator");
         System.out.println (lineSep +
           "*********** inverseBrent:   no convergence after " + MAXITER +
           " iterations");
      }
      return b;
   }


   /**
    * Computes and returns the inverse distribution function 
    * <SPAN CLASS="MATH"><I>x</I> = <I>F</I><SUP>-1</SUP>(<I>u</I>)</SPAN>,
    * using bisection. Restrictions: 
    * <SPAN CLASS="MATH"><I>u</I>&#8712;[0, 1]</SPAN>.
    * 
    * @param u value at which the inverse distribution function is evaluated
    * 
    *    @return the inverse distribution function evaluated at <TT>u</TT>
    *    @exception IllegalArgumentException if <SPAN CLASS="MATH"><I>u</I></SPAN> is  not in the interval <SPAN CLASS="MATH">[0, 1]</SPAN>
    * 
    * 
    */
   public double inverseBisection (double u) {
      final int MAXITER = 100;              // Maximum number of iterations
      final double EPSILON = EPSARRAY[decPrec];  // Absolute precision
      final double XLIM =  Double.MAX_VALUE/2.0;
      final boolean DEBUG = false;
      final String lineSep = System.getProperty("line.separator");

      if (u > 1.0 || u < 0.0)
          throw new IllegalArgumentException ("u not in [0, 1]");
      if (decPrec > Num.DBL_DIG)
          throw new IllegalArgumentException ("decPrec too large");
      if (decPrec <= 0)
          throw new IllegalArgumentException ("decPrec <= 0");
      if (DEBUG) {
//          System.out.println ("---------------------------" +
//              " -----------------------------" + lineSep +
//               PrintfFormat.f (10, 8, u));
      }

      double x = 0.0;
      if (u <= 0.0) {
          x = supportA;
          if (DEBUG) {
//             System.out.println (lineSep + "            x                   y" +
//                 lineSep + PrintfFormat.g (17, 2, x) + " " +
//                 PrintfFormat.f (17, decPrec, u));
          }
          return x;
      }
      if (u >= 1.0) {
          x = supportB;
          if (DEBUG) {
//             System.out.println (lineSep + "            x                   y" +
//                            lineSep + PrintfFormat.g (17, 2, x) + " " +
//                            PrintfFormat.f (17, decPrec, u));
          }
          return x;
      }

      double [] iv = new double [2];
      findInterval (u, iv);
      double xa = iv[0];
      double xb = iv[1];
      double yb = cdf(xb) - u;
      double ya = cdf(xa) - u;
      double y;

      if (DEBUG)
          System.out.println (lineSep +
             "iter              xa                   xb           F - u");

      boolean fini = false;
      int i = 0;
      while (!fini) {
          if (DEBUG)
//             System.out.println (PrintfFormat.d (3, i) + "  " +
//                                 PrintfFormat.g (18, decPrec, xa) + "  " +
//                                 PrintfFormat.g (18, decPrec, xb) + "  " +
//                                 PrintfFormat.g (14, 4, y));
          x = (xa + xb)/2.0;
          y = cdf(x) - u;
          if ((y == 0.0) ||
              (Math.abs ((xb - xa)/(x + Num.DBL_EPSILON)) <= EPSILON))  {
              fini = true;
//              if (DEBUG)
//                System.out.println (lineSep + "                x" +
//                     "                     u" + lineSep +
//                     PrintfFormat.g (20, decPrec, x) + "  " +
//                     PrintfFormat.g (18, decPrec, y+u));
          }
          else if (y*ya < 0.0)
             xb = x;
          else
             xa = x;
          ++i;

          if (i > MAXITER) {
                //System.out.println (lineSep +
                //  "** inverseF:SEARCH DOES NOT SEEM TO CONVERGE");
              fini = true;
          }
      }
      return x;
   }


   /**
    * Returns the inverse distribution function 
    * <SPAN CLASS="MATH"><I>x</I> = <I>F</I><SUP>-1</SUP>(<I>u</I>)</SPAN>.
    *  Restrictions: 
    * <SPAN CLASS="MATH"><I>u</I>&#8712;[0, 1]</SPAN>.
    * 
    * @param u value at which the inverse distribution function is evaluated
    * 
    *    @return the inverse distribution function evaluated at <TT>u</TT>
    *    @exception IllegalArgumentException if <SPAN CLASS="MATH"><I>u</I></SPAN> is  not in the interval <SPAN CLASS="MATH">[0, 1]</SPAN>
    * 
    * 
    */
   public double inverseF (double u) {
      double [] iv = new double [2];
      findInterval (u, iv);
      return inverseBrent (iv[0], iv[1], u, EPSARRAY[decPrec]);
   }


   /**
    * Returns the mean.
    * 
    * @return the mean
    * 
    */
   public double getMean() {
      throw new UnsupportedOperationException("getMean is not implemented ");
   }


   /**
    * Returns the variance.
    * 
    * @return the variance
    * 
    */
   public double getVariance() {
      throw new UnsupportedOperationException("getVariance is not implemented ");
   }


   /**
    * Returns the standard deviation.
    * 
    * @return the standard deviation
    * 
    */
   public double getStandardDeviation() {
      throw new UnsupportedOperationException (
         "getStandardDeviation is not implemented ");
   }



   /**
    * Returns <SPAN CLASS="MATH"><I>x</I><SUB>a</SUB></SPAN> such that the probability density is 0 everywhere
    *  outside the interval 
    * <SPAN CLASS="MATH">[<I>x</I><SUB>a</SUB>, <I>x</I><SUB>b</SUB>]</SPAN>.
    * 
    * @return lower limit of support
    * 
    */
   public double getXinf() {
      return supportA;
   }


   /**
    * Returns <SPAN CLASS="MATH"><I>x</I><SUB>b</SUB></SPAN> such that the probability density is 0 everywhere
    *  outside the interval 
    * <SPAN CLASS="MATH">[<I>x</I><SUB>a</SUB>, <I>x</I><SUB>b</SUB>]</SPAN>.
    * 
    * @return upper limit of support
    * 
    */
   public double getXsup() {
      return supportB;
   }



   /**
    * Sets the value <SPAN CLASS="MATH"><I>x</I><SUB>a</SUB> =</SPAN> <TT>xa</TT>, such that the probability
    *  density is 0 everywhere outside the interval 
    * <SPAN CLASS="MATH">[<I>x</I><SUB>a</SUB>, <I>x</I><SUB>b</SUB>]</SPAN>.
    * 
    * @param xa lower limit of support
    * 
    * 
    */
   public void setXinf (double xa) {
      supportA = xa;
   }


   /**
    * Sets the value <SPAN CLASS="MATH"><I>x</I><SUB>b</SUB> =</SPAN> <TT>xb</TT>, such that the probability
    *  density is 0 everywhere outside the interval 
    * <SPAN CLASS="MATH">[<I>x</I><SUB>a</SUB>, <I>x</I><SUB>b</SUB>]</SPAN>.
    * 
    * @param xb upper limit of support
    * 
    */
   public void setXsup (double xb) {
      supportB = xb;
   }
}
