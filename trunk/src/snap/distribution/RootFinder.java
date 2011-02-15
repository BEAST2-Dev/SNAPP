
package snap.distribution;


/**
 * This class provides methods to solve non-linear equations.
 * 
 */
public class RootFinder {
   private RootFinder() {}


   /**
    * Computes a root <SPAN CLASS="MATH"><I>x</I></SPAN> of the function in <TT>f</TT> using the
    *     Brent-Dekker method. The interval <SPAN CLASS="MATH">[<I>a</I>, <I>b</I>]</SPAN> must contain the root <SPAN CLASS="MATH"><I>x</I></SPAN>.
    *     The calculations are done with an approximate relative precision
    *     <TT>tol</TT>.  Returns <SPAN CLASS="MATH"><I>x</I></SPAN> such that <SPAN CLASS="MATH"><I>f</I> (<I>x</I>) = 0</SPAN>.
    *  
    * @param a left endpoint of initial interval
    * 
    *    @param b right endpoint of initial interval
    * 
    *    @param f the function which is evaluated
    * 
    *    @param tol accuracy goal
    * 
    *    @return the root <SPAN CLASS="MATH"><I>x</I></SPAN>
    * 
    */
   public static double brentDekker (double a, double b,
                                     MathFunction f, double tol) {
      final double EPS = 0.5E-15;
      final int MAXITER = 120;    // Maximum number of iterations
      double c, d, e;
      double fa, fb, fc;

      // Special case I = [b, a]
      if (b < a) {
         double ctemp = a;
         a = b;
         b = ctemp;
      }
      // Initialization
      fa = f.evaluate (a);
      fb = f.evaluate (b);
      c = a;
      fc = fa;
      d = e = b - a;
      tol += EPS + Num.DBL_EPSILON; // in case tol is too small

      if (Math.abs (fc) < Math.abs (fb)) {
         a = b;
         b = c;
         c = a;
         fa = fb;
         fb = fc;
         fc = fa;
      }

      for (int i = 0; i < MAXITER; i++) {
         double s, p, q, r;
         double tol1 = tol + 4.0 * Num.DBL_EPSILON * Math.abs (b);
         double xm = 0.5 * (c - b);

         if ((Math.abs (fb) == 0.0) || (Math.abs (xm) <= tol1))
            return b;

         if ((Math.abs (e) >= tol1) && (Math.abs (fa) > Math.abs (fb))) {
            if (a != c) {
               // Inverse quadratic interpolation
               q = fa / fc;
               r = fb / fc;
               s = fb / fa;
               p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
               q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            } else {
               // Linear interpolation
               s = fb / fa;
               p = 2.0 * xm * s;
               q = 1.0 - s;
            }

            // Adjust signs
            if (p > 0.0)
               q = -q;
            p = Math.abs (p);

            // Is interpolation acceptable ?
            if (((2.0 * p) >= (3.0 * xm * q - Math.abs (tol1 * q)))
                  || (p >= Math.abs (0.5 * e * q))) {
               d = xm;
               e = d;
            } else {
               e = d;
               d = p / q;
            }
         } else {
            // Bisection necessary
            d = xm;
            e = d;
         }

         a = b;
         fa = fb;
         if (Math.abs (d) > tol1)
            b += d;
         else if (xm < 0.0)
            b -= tol1;
         else
            b += tol1;
         fb = f.evaluate (b);
         if ((fb * (fc / Math.abs (fc))) > 0.0) {
            c = a;
            fc = fa;
            d = e = b - a;
         } else {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
         }
      }

      return b;
   }


   /**
    * Computes a root <SPAN CLASS="MATH"><I>x</I></SPAN> of the function in <TT>f</TT> using the
    *     <SPAN  CLASS="textit">bisection</SPAN> method. The interval <SPAN CLASS="MATH">[<I>a</I>, <I>b</I>]</SPAN> must contain the root <SPAN CLASS="MATH"><I>x</I></SPAN>.
    *     The calculations are done with an approximate relative precision
    *     <TT>tol</TT>.  Returns <SPAN CLASS="MATH"><I>x</I></SPAN> such that <SPAN CLASS="MATH"><I>f</I> (<I>x</I>) = 0</SPAN>.
    *  
    * @param a left endpoint of initial interval
    * 
    *    @param b right endpoint of initial interval
    * 
    *    @param f the function which is evaluated
    * 
    *    @param tol accuracy goal
    * 
    *    @return the root <SPAN CLASS="MATH"><I>x</I></SPAN>
    * 
    */
   public static double bisection (double a, double b,
                                   MathFunction f, double tol) {
      // Case I = [b, a]
      if (b < a) {
         double ctemp = a;
         a = b;
         b = ctemp;
      }
      double xa = a;
      double xb = b;
      double yb = f.evaluate (b);
      double ya = f.evaluate (a);
      double x = 0, y = 0;
      final int MAXITER = 1200;   // Maximum number of iterations
      final boolean DEBUG = false;
      final double myMIN_NORMAL = 2.25e-308;

      if (DEBUG)
         System.out.println
         ("\niter              xa                   xb              f(x)");

      boolean fini = false;
      int i = 0;
      while (!fini) {
         x = (xa + xb) / 2.0;
         y = f.evaluate (x);
         if ((Math.abs (y) <= myMIN_NORMAL) ||
             (Math.abs (xb - xa) <= tol * Math.abs (x)) ||
             (Math.abs (xb - xa) <= 0)) {
            return x;
         }
         if (y * ya < 0.0)
            xb = x;
         else
            xa = x;
         ++i;
         if (DEBUG)
            System.out.printf("%3d    %18.12g     %18.12g    %14.4g%n",
                              i, xa, xb, y); 
         if (i > MAXITER) {
            System.out.println ("***** bisection:  SEARCH DOES NOT CONVERGE");
            fini = true;
         }
      }
      return x;
   }

}
