
package snap.distribution;
//import umontreal.iro.lecuyer.util.Misc;
//import cern.jet.math.Bessel;

/**
 * This class provides a few constants and some methods to compute numerical
 * quantities such as factorials, combinations, gamma functions, and so on.
 * 
 */
public class Num {
   protected static final double XBIG = 40.0;
   private static final double SQPI_2 = 0.88622692545275801; // Sqrt(Pi)/2
   private static final double LOG_SQPI_2 = -0.1207822376352453; // Ln(Sqrt(Pi)/2)
   private static final double LOG4 = 1.3862943611198906;   // Ln(4)
   private static final double LOG_PI = 1.14472988584940017413; // Ln(Pi)
   private static final double PIsur2 = Math.PI/2.0;

   private static final double UNSIX = 1.0/6.0;
   private static final double QUARAN = 1.0/42.0;
   private static final double UNTRENTE = 1.0 / 30.0;
   private static final double DTIERS = 2.0 / 3.0;
   private static final double CTIERS = 5.0 / 3.0;
   private static final double STIERS = 7.0 / 3.0;
   private static final double QTIERS = 14.0 / 3.0;

   private static final double[] AERF = {
      // used in erf(x)
      1.4831105640848036E0,
      -3.0107107338659494E-1,
      6.8994830689831566E-2,
      -1.3916271264722188E-2,
      2.4207995224334637E-3,
      -3.6586396858480864E-4,
      4.8620984432319048E-5,
      -5.7492565580356848E-6,
      6.1132435784347647E-7,
      -5.8991015312958434E-8,
      5.2070090920686482E-9,
      -4.2329758799655433E-10,
      3.188113506649174974E-11,
      -2.2361550188326843E-12,
      1.46732984799108492E-13,
      -9.044001985381747E-15,
      5.25481371547092E-16
   };

   private static final double[] AERFC = {
      // used in erfc(x)
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


   private static final double[] AlnGamma = {
      /* Chebyshev coefficients for lnGamma (x + 3), 0 <= x <= 1 In Yudell
         Luke: The special functions and their approximations, Vol. II,
         Academic Press, p. 301, 1969. There is an error in the additive
         constant in the formula: (Ln (2)). */
      0.52854303698223459887,
      0.54987644612141411418,
      0.02073980061613665136,
     -0.00056916770421543842,
      0.00002324587210400169,
     -0.00000113060758570393,
      0.00000006065653098948,
     -0.00000000346284357770,
      0.00000000020624998806,
     -0.00000000001266351116,
      0.00000000000079531007,
     -0.00000000000005082077,
      0.00000000000000329187,
     -0.00000000000000021556,
      0.00000000000000001424,
     -0.00000000000000000095
   };

   private Num() {}


   /**
    * Difference between 1.0 and the smallest <TT>double</TT> greater than 1.0.
    * 
    */
   public static final double DBL_EPSILON = 2.2204460492503131e-16;


   /**
    * Largest <TT>int</TT> <SPAN CLASS="MATH"><I>x</I></SPAN> such that <SPAN CLASS="MATH">2<SUP>x-1</SUP></SPAN> is representable
    *   (approximately) as a <TT>double</TT>.
    * 
    */
   public static final int DBL_MAX_EXP = 1024;


   /**
    * Smallest <TT>int</TT> <SPAN CLASS="MATH"><I>x</I></SPAN> such that <SPAN CLASS="MATH">2<SUP>x-1</SUP></SPAN> is representable
    *   (approximately) as a normalised <TT>double</TT>.
    * 
    */
   public static final int DBL_MIN_EXP = -1021;


   /**
    * Largest <TT>int</TT> <SPAN CLASS="MATH"><I>x</I></SPAN> such that <SPAN CLASS="MATH">10<SUP>x</SUP></SPAN> is representable
    *    (approximately) as a <TT>double</TT>.
    * 
    */
   public static final int DBL_MAX_10_EXP = 308;


   /**
    * Smallest normalized positive floating-point <TT>double</TT>.
    * 
    */
   public static final double DBL_MIN = 2.2250738585072014e-308;


   /**
    * Natural logarithm of <TT>DBL_MIN</TT>.
    * 
    */
   public static final double LN_DBL_MIN = -708.3964185322641;


   /**
    * Number of decimal digits of precision in a <TT>double</TT>.
    * 
    */
   public static final int DBL_DIG = 15;


   /**
    * The constant <SPAN CLASS="MATH"><I>e</I></SPAN>.
    * 
    */
   public static final double EBASE = 2.7182818284590452354;


   /**
    * The Euler-Mascheroni constant.
    * 
    */
   public static final double EULER = 0.57721566490153286;


   /**
    * The value of <SPAN CLASS="MATH">(2)<SUP>1/2</SUP></SPAN>.
    * 
    */
   public static final double RAC2 = 1.41421356237309504880;


   /**
    * The value of 
    * <SPAN CLASS="MATH">1/(2)<SUP>1/2</SUP></SPAN>.
    * 
    */
   public static final double IRAC2 = 0.70710678118654752440;


   /**
    * The values of <SPAN CLASS="MATH">ln&nbsp;2</SPAN>.
    * 
    */
   public static final double LN2 = 0.69314718055994530941;


   /**
    * The values of <SPAN CLASS="MATH">1/ln&nbsp;2</SPAN>.
    * 
    */
   public static final double ILN2 = 1.44269504088896340737;


   /**
    * Largest integer 
    * <SPAN CLASS="MATH"><I>n</I><SUB>0</SUB> = 2<SUP>53</SUP></SPAN> such that any integer
    *   <SPAN CLASS="MATH"><I>n</I>&nbsp;&lt;=&nbsp;<I>n</I><SUB>0</SUB></SPAN> is represented  exactly as a <TT>double</TT>.
    * 
    */
   public static final double MAXINTDOUBLE = 9007199254740992.0;


   /**
    * Powers of 2 up to <TT>MAXTWOEXP</TT> are stored exactly
    *     in the array <TT>TWOEXP</TT>.
    * 
    */
   public static final double MAXTWOEXP = 64;


   /**
    * Contains the precomputed positive powers of 2.
    *    One has <TT>TWOEXP[j]</TT><SPAN CLASS="MATH">= 2<SUP>j</SUP></SPAN>, for 
    * <SPAN CLASS="MATH"><I>j</I> = 0,..., 64</SPAN>.
    * 
    */
   public static final double TWOEXP[] = {
      1.0, 2.0, 4.0, 8.0, 1.6e1, 3.2e1,
      6.4e1, 1.28e2, 2.56e2, 5.12e2, 1.024e3,
      2.048e3, 4.096e3, 8.192e3, 1.6384e4, 3.2768e4,
      6.5536e4, 1.31072e5, 2.62144e5, 5.24288e5,
      1.048576e6, 2.097152e6, 4.194304e6, 8.388608e6,
      1.6777216e7, 3.3554432e7, 6.7108864e7,
      1.34217728e8, 2.68435456e8, 5.36870912e8,
      1.073741824e9, 2.147483648e9, 4.294967296e9,
      8.589934592e9, 1.7179869184e10, 3.4359738368e10,
      6.8719476736e10, 1.37438953472e11, 2.74877906944e11,
      5.49755813888e11, 1.099511627776e12, 2.199023255552e12,
      4.398046511104e12, 8.796093022208e12,
      1.7592186044416e13, 3.5184372088832e13,
      7.0368744177664e13, 1.40737488355328e14,
      2.81474976710656e14, 5.62949953421312e14,
      1.125899906842624e15, 2.251799813685248e15,
      4.503599627370496e15, 9.007199254740992e15,
      1.8014398509481984e16, 3.6028797018963968e16,
      7.2057594037927936e16, 1.44115188075855872e17,
      2.88230376151711744e17, 5.76460752303423488e17,
      1.152921504606846976e18, 2.305843009213693952e18,
      4.611686018427387904e18, 9.223372036854775808e18,
      1.8446744073709551616e19
     };



   /**
    * Contains the precomputed negative powers of 10.
    *    One has <TT>TEN_NEG_POW[j]</TT><SPAN CLASS="MATH">= 10<SUP>-j</SUP></SPAN>, for 
    * <SPAN CLASS="MATH"><I>j</I> = 0,&#8230;, 16</SPAN>.
    * 
    */
   public static final double TEN_NEG_POW[] = {
      1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8,
      1.0e-9, 1.0e-10, 1.0e-11, 1.0e-12, 1.0e-13, 1.0e-14, 1.0e-15, 1.0e-16
     };



   /**
    * Returns the greatest common divisor (gcd) of <SPAN CLASS="MATH"><I>x</I></SPAN> and <SPAN CLASS="MATH"><I>y</I></SPAN>.
    *  
    *  @param x integer
    * 
    *        @param y integer
    * 
    *        @return the GCD of <SPAN CLASS="MATH"><I>x</I></SPAN> and <SPAN CLASS="MATH"><I>y</I></SPAN>
    *  
    */
   public static int gcd (int x, int y) {
      if (x < 0) x = -x;
      if (y < 0) y = -y;
      int r;
      while (y != 0) {
         r = x % y;
         x = y;
         y = r;
      }
      return x;
   }


   /**
    * Returns the greatest common divisor (gcd) of <SPAN CLASS="MATH"><I>x</I></SPAN> and <SPAN CLASS="MATH"><I>y</I></SPAN>.
    *  
    *  @param x integer
    * 
    *        @param y integer
    * 
    *        @return the GCD of <SPAN CLASS="MATH"><I>x</I></SPAN> and <SPAN CLASS="MATH"><I>y</I></SPAN>
    *  
    */
   public static long gcd (long x, long y) {
      if (x < 0) x = -x;
      if (y < 0) y = -y;
      long r;
      while (y != 0) {
         r = x % y;
         x = y;
         y = r;
      }
      return x;
   }


   /**
    * Returns the number of different combinations
    *    of <SPAN CLASS="MATH"><I>s</I></SPAN> objects amongst <SPAN CLASS="MATH"><I>n</I></SPAN>.  Uses an algorithm that prevents overflows
    *    (when computing factorials), if possible.
    *  
    *  @param n total number of objects
    * 
    *       @param s number of chosen objects on a combination
    * 
    *       @return the number of different combinations of <SPAN CLASS="MATH"><I>s</I></SPAN> objects amongst <SPAN CLASS="MATH"><I>n</I></SPAN>
    *  
    */
   public static double combination (int n, int s) {
      final int NLIM = 100;      // pour eviter les debordements
      int i;
      if (s == 0 || s == n)
         return 1;
      if (s < 0) {
         System.err.println ("combination:   s < 0");
         return 0;
      }
      if (s > n) {
         System.err.println ("combination:   s > n");
         return 0;
      }
      if (s > (n/2))
         s = n - s;
      if (n <= NLIM) {
         double Res = 1.0;
         int Diff = n - s;
         for (i = 1; i <= s; i++) {
            Res *= (double)(Diff + i)/(double)i;
         }
         return Res;
      }
      else {
         double Res = (lnFactorial (n) - lnFactorial (s))
            - lnFactorial (n - s);
         return Math.exp (Res);
      }
   }



   /**
    * Returns the value of factorial <SPAN CLASS="MATH"><I>n</I></SPAN>.
    *  
    *  @param n the integer for which the factorial must be computed
    * 
    *        @return the value of <SPAN CLASS="MATH"><I>n</I>!</SPAN>
    *  
    */
   public static double factorial (int n) {
      if (n < 0)
        throw new IllegalArgumentException ("factorial:   n < 0");
      double T = 1;
      for (int j = 2; j <= n; j++)
         T *= j;
      return T;
   }


   /**
    * Returns the value of the natural logarithm of
    *   factorial <SPAN CLASS="MATH"><I>n</I></SPAN>. Gives 16 decimals of precision
    *   (relative error 
    * <SPAN CLASS="MATH">&lt; 0.5&#215;10<SUP>-15</SUP></SPAN>).
    *  
    *  @param n the integer for which the log-factorial has to be computed
    * 
    *        @return the natural logarithm of <SPAN CLASS="MATH"><I>n</I></SPAN> factorial
    *  
    */
   public static double lnFactorial (int n) {
      final int NLIM = 14;

      if (n < 0)
        throw new IllegalArgumentException ("lnFactorial:   n < 0");

      if (n == 0 || n == 1)
         return 0.0;
      if (n <= NLIM) {
         long z = 1;
         long x = 1;
         for (int i = 2; i <= n; i++) {
            ++x;
            z *= x;
         }
         return Math.log (z);
      }
      else {
         double x = (double)(n + 1);
         double y = 1.0/(x*x);
         double z = ((-(5.95238095238E-4*y) + 7.936500793651E-4)*y -
            2.7777777777778E-3)*y + 8.3333333333333E-2;
         z = ((x - 0.5)*Math.log (x) - x) + 9.1893853320467E-1 + z/x;
         return z;
      }
   }


   /**
    * Computes and returns the Stirling numbers of the second kind
    * 
    * @param m number of rows of the allocated matrix
    * 
    *        @param n number of columns of the allocated matrix
    * 
    *        @return the matrix of Stirling numbers
    *  
    */
   public static double[][] calcMatStirling (int m, int n) {
      int i, j, k;
      double[][] M = new double[m+1][n+1];

      for (i = 0; i <= m; i++)
         for (j = 0; j <= n; j++)
            M[i][j] = 0.0;

      M[0][0] = 1.0;
      for (j = 1; j <= n; j++) {
         M[0][j] = 0.0;
         if (j <= m) {
            k = j - 1;
            M[j][j] = 1.0;
         }
         else
            k = m;
         for (i = 1; i <= k; i++)
            M[i][j] = (double)i*M[i][j - 1] + M[i - 1][j - 1];
      }
      return M;
   }


   /**
    * Returns <SPAN CLASS="MATH">log<SUB>2</SUB>(</SPAN><TT>x</TT><SPAN CLASS="MATH">)</SPAN>.
    *  
    *  @param x the value for which the logarithm must be computed
    * 
    *        @return the value of <SPAN CLASS="MATH">log<SUB>2</SUB>(</SPAN><TT>x</TT><SPAN CLASS="MATH">)</SPAN>
    *  
    */
   public static double log2 (double x) {
     return ILN2*Math.log (x);
   }


   /**
    * Returns the natural logarithm of the gamma function <SPAN CLASS="MATH"><I>&#915;</I>(<I>x</I>)</SPAN>
    *    evaluated at <TT>x</TT>.
    *    Gives 16 decimals of precision, but is implemented only for <SPAN CLASS="MATH"><I>x</I> &gt; 0</SPAN>.
    *   
    *   @param x the value for which the lnGamma function must be computed
    * 
    *        @return the natural logarithm of the gamma function
    *  
    */
   public static double lnGamma (double x) {
      if (x <= 0.0)
         throw new IllegalArgumentException ("lnGamma:   x <= 0");

      final double XLIMBIG = 1.0/DBL_EPSILON;
      final double XLIM1 = 18.0;
      final double DK2 = 0.91893853320467274177;     // Ln (sqrt (2 Pi))
      final double DK1 = 0.9574186990510627;
      final int N = 15;              // Degree of Chebyshev polynomial
      double y, z;
      int i, k;

/*
      if (x <= 0.0) {
         double f = (1.0 - x) - Math.floor (1.0 - x);
         return LOG_PI - lnGamma (1.0 - x) - Math.log (Math.sin (Math.PI * f));
      }
*/
      if (x > XLIM1) {
         if (x > XLIMBIG)
            y = 0.0;
         else
            y = 1.0/(x*x);
         z = ((-(5.95238095238E-4*y) + 7.936500793651E-4)*y -
            2.7777777777778E-3)*y + 8.3333333333333E-2;
         z = ((x - 0.5)*Math.log (x) - x) + DK2 + z/x;
         return z;

      } else if (x > 4.0) {
         k = (int)x;
         z = x - k;
         y = 1.0;
         for (i = 3; i < k; i++)
            y *= z + i;
         y = Math.log (y);

      } else if (x <= 0.0)
         return Double.MAX_VALUE;

      else if (x < 3.0) {
         k = (int)x;
         z = x - k;
         y = 1.0;
         for (i = 2; i >= k; i--)
            y *= z + i;
         y = -Math.log (y);
      }
      else {           // 3 <= x <= 4
         z = x - 3.0;
         y = 0.0;
      }
      z = evalCheby (AlnGamma, N, 2.0*z - 1.0);
      return z + DK1 + y;
   }


   /**
    * Computes the natural logarithm of the Beta function
    *  
    * <SPAN CLASS="MATH"><I>B</I>(<I>&#955;</I>, <I>&#957;</I>)</SPAN>.  It is defined in terms of the Gamma function as
    *  <P>
    * </P>
    * <DIV ALIGN="CENTER" CLASS="mathdisplay">
    * <I>B</I>(<I>&#955;</I>, <I>&#957;</I>) = <IMG
    *  ALIGN="MIDDLE" BORDER="0" SRC="Numimg1.png"
    *  ALT="$\displaystyle {\frac{{\Gamma(\lambda)\Gamma(\nu)}}{{\Gamma(\lambda + \nu)}}}$">
    * </DIV><P></P>
    * with <TT>lam</TT> <SPAN CLASS="MATH">= <I>&#955;</I></SPAN> and  <TT>nu</TT> <SPAN CLASS="MATH">= <I>&#957;</I></SPAN>.
    * 
    */
   public static double lnBeta (double lam, double nu) {
      return lnGamma (lam) + lnGamma (nu) - lnGamma (lam + nu);
   }


   /**
    * Returns the value of the logarithmic derivative of the Gamma function
    *    
    * <SPAN CLASS="MATH"><I>&#968;</I>(<I>x</I>) = <I>&#915;</I>'(<I>x</I>)/<I>&#915;</I>(<I>x</I>)</SPAN>.
    * 
    */
   public static double digamma (double x) {
      final double C7[][] = {
       {1.3524999667726346383e4, 4.5285601699547289655e4, 4.5135168469736662555e4,
        1.8529011818582610168e4, 3.3291525149406935532e3, 2.4068032474357201831e2,
        5.1577892000139084710, 6.2283506918984745826e-3},
       {6.9389111753763444376e-7, 1.9768574263046736421e4, 4.1255160835353832333e4,
          2.9390287119932681918e4, 9.0819666074855170271e3,
          1.2447477785670856039e3, 6.7429129516378593773e1, 1.0}
      };
      final double C4[][] = {
       {-2.728175751315296783e-15, -6.481571237661965099e-1, -4.486165439180193579,
        -7.016772277667586642, -2.129404451310105168},
       {7.777885485229616042, 5.461177381032150702e1,
        8.929207004818613702e1, 3.227034937911433614e1, 1.0}
      };

      double prodPj = 0.0;
      double prodQj = 0.0;
      double digX = 0.0;

      if (x >= 3.0) {
         double x2 = 1.0 / (x * x);
         for (int j = 4; j >= 0; j--) {
            prodPj = prodPj * x2 + C4[0][j];
            prodQj = prodQj * x2 + C4[1][j];
         }
         digX = Math.log (x) - (0.5 / x) + (prodPj / prodQj);

      } else if (x >= 0.5) {
         final double X0 = 1.46163214496836234126;
         for (int j = 7; j >= 0; j--) {
            prodPj = x * prodPj + C7[0][j];
            prodQj = x * prodQj + C7[1][j];
         }
         digX = (x - X0) * (prodPj / prodQj);

      } else {
         double f = (1.0 - x) - Math.floor (1.0 - x);
         digX = digamma (1.0 - x) + Math.PI / Math.tan (Math.PI * f);
      }

      return digX;
   }


   /**
    * Returns the value of the trigamma function 
    * <SPAN CLASS="MATH"><I>d&#968;</I>(<I>x</I>)/<I>dx</I></SPAN>, the derivative of
    *    the digamma function, evaluated at <SPAN CLASS="MATH"><I>x</I></SPAN>.
    * 
    */
   public static double trigamma (double x) {
      double y, sum;

      if (x < 0.5) {
         y = (1.0 - x) - Math.floor (1.0 - x);
         sum = Math.PI / Math.sin (Math.PI * y);
         return  sum * sum - trigamma (1.0 - x);
      }

      if (x >= 40.0) {
         // Asymptotic series
         y = 1.0/(x*x);
         sum = 1.0 + y*(1.0/6.0 - y*(1.0/30.0 - y*(1.0/42.0 - 1.0/30.0*y)));
         sum += 0.5/x;
         return sum/x;
      }

      int i;
      int p = (int) x;
      y = x - p;
      sum = 0.0;

      if (p > 3) {
         for (i = 3; i < p; i++)
            sum -= 1.0/((y + i)*(y + i));

      } else if (p < 3) {
         for (i = 2; i >= p; i--)
            sum += 1.0/((y + i)*(y + i));
      }

      /* Chebyshev coefficients for trigamma (x + 3), 0 <= x <= 1. In Yudell
         Luke: The special functions and their approximations, Vol. II,
         Academic Press, p. 301, 1969. */
      final int N = 15;
      final double A[] = { 2.0*0.33483869791094938576, -0.05518748204873009463,
         0.00451019073601150186, -0.00036570588830372083,
         2.943462746822336e-5, -2.35277681515061e-6, 1.8685317663281e-7,
         -1.475072018379e-8, 1.15799333714e-9, -9.043917904e-11,
         7.029627e-12, -5.4398873e-13, 0.4192525e-13, -3.21903e-15, 0.2463e-15,
        -1.878e-17, 0., 0. };

      return sum + evalChebyStar (A, N, y);
   }


   /**
    * Returns the value of the tetragamma function 
    * <SPAN CLASS="MATH"><I>d</I><SUP>2</SUP><I>&#968;</I>(<I>x</I>)/<I>d</I><SUP>2</SUP><I>x</I></SPAN>, the second
    *    derivative of the digamma function, evaluated at <SPAN CLASS="MATH"><I>x</I></SPAN>.
    * 
    */
   public static double tetragamma (double x) {
      double y, sum;

      if (x < 0.5) {
         y = (1.0 - x) - Math.floor (1.0 - x);
         sum = Math.PI / Math.sin (Math.PI * y);
         return 2.0 * Math.cos (Math.PI * y) * sum * sum * sum +
               tetragamma (1.0 - x);
      }

      if (x >= 20.0) {
         // Asymptotic series
         y = 1.0/(x*x);
         sum = y*(0.5 - y*(1.0/6.0 - y*(1.0/6.0 - y*(0.3 - 5.0/6.0*y))));
         sum += 1.0 + 1.0/x;
         return -sum*y;
      }

      int i;
      int p = (int) x;
      y = x - p;
      sum = 0.0;

      if (p > 3) {
         for (i = 3; i < p; i++)
            sum += 1.0 / ((y + i) * (y + i) * (y + i));

      } else if (p < 3) {
         for (i = 2; i >= p; i--)
            sum -= 1.0 / ((y + i) * (y + i) * (y + i));
      }

      /* Chebyshev coefficients for tetragamma (x + 3), 0 <= x <= 1. In Yudell
         Luke: The special functions and their approximations, Vol. II,
         Academic Press, p. 301, 1969. */
      final int N = 16;
      final double A[] = { -0.11259293534547383037*2.0, 0.03655700174282094137,
         -0.00443594249602728223, 0.00047547585472892648,
         -4.747183638263232e-5, 4.52181523735268e-6, -4.1630007962011e-7,
         3.733899816535e-8, -3.27991447410e-9, 2.8321137682e-10,
         -2.410402848e-11, 2.02629690e-12, -1.6852418e-13, 1.388481e-14,
         -1.13451e-15, 9.201e-17, -7.41e-18, 5.9e-19, -5.0e-20 };

      return 2.0 * sum + evalChebyStar (A, N, y);
   }


   /**
    * Returns the value of the ratio 
    * <SPAN CLASS="MATH"><I>&#915;</I>(<I>x</I> + 1/2)/<I>&#915;</I>(<I>x</I>)</SPAN> of two gamma
    * functions, evaluated in a numerically stable way. Restriction: <SPAN CLASS="MATH"><I>x</I> &gt; 0</SPAN>.
    * 
    */
   public static double gammaRatioHalf (double x) {
      if (x <= 0.0)
         throw new IllegalArgumentException ("gammaRatioHalf:   x <= 0");

      if (x <= 10.0) {
         double y = lnGamma (x + 0.5) - lnGamma (x);
         return Math.exp(y);
      }

      double sum;
      if (x <= 300.0) {
         // The sum converges very slowly for small x, but faster as x increases
         final double EPSILON = 1.0e-15;
         double term = 1.0;
         sum = 1.0;
         int i = 1;
         while (term > EPSILON*sum) {
            term *= (i - 1.5)*(i - 1.5) /(i*(x + i - 1.5));
            sum += term;
            i++;
         }
         return Math.sqrt ((x - 0.5)*sum);
      }

      // Asymptotic series for Gamma(x + 0.5) / (Gamma(x) * Sqrt(x))
      // Comparer la vitesse de l'asymptotique avec la somme ci-dessus !!!
      double y = 1.0 / (8.0*x);
      sum = 1.0 + y*(-1.0 + y*(0.5 + y*(2.5 - y*(2.625 + 49.875*y))));
      return sum*Math.sqrt(x);
   }


   /**
    * Computes the <SPAN CLASS="MATH"><I>n</I></SPAN>-th harmonic number 
    * <SPAN CLASS="MATH"><I>H</I><SUB>n</SUB> = &sum;<SUB>j=1</SUB><SUP>n</SUP>1/<I>j</I></SPAN>.
    * 
    */
   public static double harmonic (long n) {
      if (n < 1)
         throw new IllegalArgumentException ("n < 1");
      return digamma(n + 1) + EULER;
   }


   /**
    * .
    * 
    * Computes the sum
    * <P>
    * </P>
    * <DIV ALIGN="CENTER" CLASS="mathdisplay">
    * <IMG
    *  ALIGN="MIDDLE" BORDER="0" SRC="Numimg2.png"
    *  ALT="$\displaystyle \sideset{}{'}\htsum_{{-n/2&lt;j\le n/2}}^{}$"> &nbsp;<IMG
    *  ALIGN="MIDDLE" BORDER="0" SRC="Numimg3.png"
    *  ALT="$\displaystyle {\frac{1}{{\vert j\vert}}}$">,
    * </DIV><P></P>
    * where the symbol 
    * <SPAN CLASS="MATH">&sum;<SUP>&#8242;</SUP></SPAN> means that the term with <SPAN CLASS="MATH"><I>j</I> = 0</SPAN> is excluded
    *  from the sum.
    * 
    */
   public static double harmonic2 (long n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (1 == n)
         return 0.0;
      if (2 == n)
         return 1.0;

      long k = n / 2;
      if ((n & 1) == 1)
         return  2.0*harmonic(k);         // n odd
      return  1.0/k + 2.0*harmonic(k-1);  // n even
   }


   /**
    * Returns the volume <SPAN CLASS="MATH"><I>V</I></SPAN> of a sphere of radius 1 in <SPAN CLASS="MATH"><I>t</I></SPAN> dimensions
    *   using the norm <SPAN CLASS="MATH"><I>L</I><SUB>p</SUB></SPAN>. It is given by the formula
    * <P>
    * </P>
    * <DIV ALIGN="CENTER" CLASS="mathdisplay">
    * <I>V</I> = ([2<I>&#915;</I>(1 + 1/<I>p</I>)]<SUP>t</SUP>)/<I>&#915;</I>(1 + <I>t</I>/<I>p</I>),&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<I>p</I> &gt; 0,
    * </DIV><P></P>
    * where <SPAN CLASS="MATH"><I>&#915;</I></SPAN> is the gamma function.
    *   The case of the sup norm <SPAN CLASS="MATH"><I>L</I><SUB>&#8734;</SUB></SPAN> is  obtained by choosing <SPAN CLASS="MATH"><I>p</I> = 0</SPAN>.
    *   Restrictions: <SPAN CLASS="MATH"><I>p</I>&nbsp;&gt;=&nbsp; 0</SPAN> and <SPAN CLASS="MATH"><I>t</I>&nbsp;&gt;=&nbsp;1</SPAN>.
    *   
    *   @param p index of the used norm
    * 
    *        @param t number of dimensions
    * 
    *        @return the volume of a sphere
    *  
    */
   public static double volumeSphere (double p, int t) {
      final double EPS = 2.0*DBL_EPSILON;
      int pLR = (int)p;
      double kLR = (double)t;
      double Vol;
      int s;

      if (p < 0)
         throw new IllegalArgumentException ("volumeSphere:   p < 0");

      if (Math.abs (p - pLR) <= EPS) {
         switch (pLR) {
         case 0:
            return TWOEXP[t];
         case 1:
            return TWOEXP[t]/(double)factorial (t);
         case 2:
            if ((t % 2) == 0)
               return Math.pow (Math.PI, kLR/2.0)/(double)factorial (t/2);
            else {
               s = (t + 1)/2;
               return Math.pow (Math.PI, (double)s - 1.0)*factorial (s)*
                  TWOEXP[2*s]/(double)factorial (2*s);
            }
          default:
         }
      }
      Vol = kLR*(LN2 + lnGamma (1.0 + 1.0/p)) -
      lnGamma (1.0 + kLR/p);
      return Math.exp (Vol);
   }


   /**
    * Evaluates the Bernoulli polynomial <SPAN CLASS="MATH"><I>B</I><SUB>n</SUB>(<I>x</I>)</SPAN> of degree <SPAN CLASS="MATH"><I>n</I></SPAN>
    *   at <SPAN CLASS="MATH"><I>x</I></SPAN>. Only degrees <SPAN CLASS="MATH"><I>n</I>&nbsp;&lt;=&nbsp;8</SPAN> are programmed for now.
    *  The first Bernoulli polynomials of even degree are:
    * <BR>
    * <DIV ALIGN="CENTER" CLASS="mathdisplay"><A NAME="bernoulli"></A>
    * <TABLE CELLPADDING="0" ALIGN="CENTER" WIDTH="100%">
    * <TR VALIGN="MIDDLE"><TD NOWRAP WIDTH="50%" ALIGN="RIGHT"><I>B</I><SUB>0</SUB>(<I>x</I>)</TD>
    * <TD WIDTH="10" ALIGN="CENTER" NOWRAP>=</TD>
    * <TD ALIGN="LEFT" NOWRAP WIDTH="50%">1</TD>
    * <TD CLASS="eqno" WIDTH=10 ALIGN="RIGHT">
    * &nbsp;</TD></TR>
    * <TR VALIGN="MIDDLE"><TD NOWRAP WIDTH="50%" ALIGN="RIGHT"><I>B</I><SUB>2</SUB>(<I>x</I>)</TD>
    * <TD WIDTH="10" ALIGN="CENTER" NOWRAP>=</TD>
    * <TD ALIGN="LEFT" NOWRAP WIDTH="50%"><I>x</I><SUP>2</SUP> - <I>x</I> + 1/6</TD>
    * <TD CLASS="eqno" WIDTH=10 ALIGN="RIGHT">
    * &nbsp;</TD></TR>
    * <TR VALIGN="MIDDLE"><TD NOWRAP WIDTH="50%" ALIGN="RIGHT"><I>B</I><SUB>4</SUB>(<I>x</I>)</TD>
    * <TD WIDTH="10" ALIGN="CENTER" NOWRAP>=</TD>
    * <TD ALIGN="LEFT" NOWRAP WIDTH="50%"><I>x</I><SUP>4</SUP> -2<I>x</I><SUP>3</SUP> + <I>x</I><SUP>2</SUP> - 1/30</TD>
    * <TD CLASS="eqno" WIDTH=10 ALIGN="RIGHT">
    * &nbsp;</TD></TR>
    * <TR VALIGN="MIDDLE"><TD NOWRAP WIDTH="50%" ALIGN="RIGHT"><I>B</I><SUB>6</SUB>(<I>x</I>)</TD>
    * <TD WIDTH="10" ALIGN="CENTER" NOWRAP>=</TD>
    * <TD ALIGN="LEFT" NOWRAP WIDTH="50%"><I>x</I><SUP>6</SUP> -3<I>x</I><SUP>5</SUP> +5<I>x</I><SUP>4</SUP>/2 - <I>x</I><SUP>2</SUP>/2 + 1/42</TD>
    * <TD CLASS="eqno" WIDTH=10 ALIGN="RIGHT">
    * &nbsp;</TD></TR>
    * <TR VALIGN="MIDDLE"><TD NOWRAP WIDTH="50%" ALIGN="RIGHT"><I>B</I><SUB>8</SUB>(<I>x</I>)</TD>
    * <TD WIDTH="10" ALIGN="CENTER" NOWRAP>=</TD>
    * <TD ALIGN="LEFT" NOWRAP WIDTH="50%"><I>x</I><SUP>8</SUP> -4<I>x</I><SUP>7</SUP> +14<I>x</I><SUP>6</SUP>/3 - 7<I>x</I><SUP>4</SUP>/3 + 2<I>x</I><SUP>2</SUP>/3 - 1/30.</TD>
    * <TD CLASS="eqno" WIDTH=10 ALIGN="RIGHT">
    * &nbsp;</TD></TR>
    * </TABLE></DIV>
    * <BR CLEAR="ALL">
    * 
    */
   public static double bernoulliPoly (int n, double x)  {
      switch (n) {
      case 0:
         return 1.0;
      case 1:
         return x - 0.5;
      case 2:
         return x*(x - 1.0) + UNSIX;
      case 3:
         return ((2.0*x - 3.0) * x + 1.0)*x*0.5;
      case 4:
         return ((x - 2.0) * x + 1.0)*x*x - UNTRENTE;
      case 5:
         return (((x - 2.5) * x + CTIERS) *x*x - UNSIX) * x;
      case 6:
         return (((x - 3.0) * x + 2.5) * x*x - 0.5) * x*x + QUARAN;
      case 7:
         return ((((x - 3.5) * x + 3.5) *x*x - 7.0/6.0) * x*x + UNSIX) * x;
      case 8:
         return ((((x - 4.0) * x +
                  QTIERS) * x*x - STIERS) * x*x + DTIERS) * x*x - UNTRENTE;
      default:
         throw new IllegalArgumentException("n must be <= 8");
      }
    //  return 0;
    }


   /**
    * Evaluates a series of Chebyshev polynomials <SPAN CLASS="MATH"><I>T</I><SUB>j</SUB></SPAN> at
    *   <SPAN CLASS="MATH"><I>x</I></SPAN> over the basic interval <SPAN CLASS="MATH">[- 1, &nbsp;1]</SPAN>. It uses
    *    the method of Clenshaw, i.e., computes and  returns
    *   <P>
    * </P>
    * <DIV ALIGN="CENTER" CLASS="mathdisplay">
    * <I>y</I> = <IMG
    *  ALIGN="MIDDLE" BORDER="0" SRC="Numimg4.png"
    *  ALT="$\displaystyle {\frac{{a_0}}{2}}$"> + &sum;<SUB>j=1</SUB><SUP>n</SUP><I>a</I><SUB>j</SUB><I>T</I><SUB>j</SUB>(<I>x</I>).
    * </DIV><P></P>
    * @param a coefficients of the polynomials
    * 
    *        @param n largest degree of polynomials
    * 
    *        @param x the parameter of the <SPAN CLASS="MATH"><I>T</I><SUB>j</SUB></SPAN> functions
    * 
    *        @return  the value of a series of Chebyshev polynomials <SPAN CLASS="MATH"><I>T</I><SUB>j</SUB></SPAN>.
    * 
    */
   public static double evalCheby (double a[], int n, double x)  {
      if (Math.abs (x) > 1.0)
         System.err.println ("Chebychev polynomial evaluated "+
                               "at x outside [-1, 1]");
      final double xx = 2.0*x;
      double b0 = 0.0, b1 = 0.0, b2 = 0.0;
      for (int j = n; j >= 0; j--) {
         b2 = b1;
         b1 = b0;
         b0 = (xx*b1 - b2) + a[j];
      }
      return (b0 - b2)/2.0;
   }


   /**
    * Evaluates a series of shifted Chebyshev polynomials <SPAN CLASS="MATH"><I>T</I><SUB>j</SUB><SUP>*</SUP></SPAN>
    *    at <SPAN CLASS="MATH"><I>x</I></SPAN> over the basic interval <SPAN CLASS="MATH">[0, &nbsp;1]</SPAN>. It uses
    *    the method of Clenshaw, i.e., computes and  returns
    *   <P>
    * </P>
    * <DIV ALIGN="CENTER" CLASS="mathdisplay">
    * <I>y</I> = [tex2html_wrap_indisplay679] + &sum;<SUB>j=1</SUB><SUP>n</SUP><I>a</I><SUB>j</SUB><I>T</I><SUB>j</SUB><SUP>*</SUP>(<I>x</I>).
    * </DIV><P></P>
    * @param a coefficients of the polynomials
    * 
    *        @param n largest degree of polynomials
    * 
    *        @param x the parameter of the <SPAN CLASS="MATH"><I>T</I><SUB>j</SUB><SUP>*</SUP></SPAN> functions
    * 
    *        @return  the value of a series of Chebyshev polynomials <SPAN CLASS="MATH"><I>T</I><SUB>j</SUB><SUP>*</SUP></SPAN>.
    * 
    */
   public static double evalChebyStar (double a[], int n, double x)  {
      if (x > 1.0 || x < 0.0)
         System.err.println ("Shifted Chebychev polynomial evaluated " +
                             "at x outside [0, 1]");
      final double xx = 2.0*(2.0*x - 1.0);
      double b0 = 0.0, b1 = 0.0, b2 = 0.0;
      for (int j = n; j >= 0; j--) {
         b2 = b1;
         b1 = b0;
         b0 = xx*b1 - b2 + a[j];
      }
      return (b0 - b2)/2.0;
   }


   /**
    * Returns the value of 
    * <SPAN CLASS="MATH"><I>K</I><SUB>1/4</SUB>(<I>x</I>)</SPAN>, where <SPAN CLASS="MATH"><I>K</I><SUB>a</SUB></SPAN> is the modified
    *   Bessel's function of the second kind.
    *   The relative error on the returned value is less than
    *   
    * <SPAN CLASS="MATH">0.5&#215;10<SUP>-6</SUP></SPAN> for 
    * <SPAN CLASS="MATH"><I>x</I> &gt; 10<SUP>-300</SUP></SPAN>.
    *  
    * @param x value at which the function is calculated
    * 
    *        @return the value of 
    * <SPAN CLASS="MATH"><I>K</I><SUB>1/4</SUB>(<I>x</I>)</SPAN>
    * 
    */
   public static double besselK025 (double x) {
      final int DEG = 6;
      double rac;
      double xx;
      double temp;
      double Res;
      double C;
      double B;
      int j;
      final double c[] = {
         32177591145.0,
         2099336339520.0,
         16281990144000.0,
         34611957596160.0,
         26640289628160.0,
         7901666082816.0,
         755914244096.0
      };

      final double b[] = {
         75293843625.0,
         2891283595200.0,
         18691126272000.0,
         36807140966400.0,
         27348959232000.0,
         7972533043200.0,
         755914244096.0
      };

      if (x < 1.E-300)
         return Double.MIN_VALUE;

      /*------------------------------------------------------------------
       * x > 0.6 => approximation asymptotique rationnelle dans Luke:
       * Yudell L.Luke "Mathematical functions and their approximations",
       * Academic Press Inc. New York, 1975, p.371
       *------------------------------------------------------------------*/
      if (x >= 0.6) {
         B = b[DEG];
         C = c[DEG];
         for (j = DEG; j >= 1; j--) {
            B = B * x + b[j - 1];
            C = C * x + c[j - 1];
         }
         Res = Math.sqrt (Math.PI / (2.0 * x)) * Math.exp (-x) * (C / B);
         return Res;
      }

      /*------------------------------------------------------------------
       * x < 0.6 => la serie de K_{1/4} = Pi/Sqrt (2) [I_{-1/4} - I_{1/4}]
       *------------------------------------------------------------------*/
      xx = x * x;
      rac = Math.pow (x/2.0, 0.25);
      Res = (((xx/1386.0 + 1.0/42.0)*xx + 1.0/3.0)*xx + 1.0)/
              (1.225416702465177*rac);
      temp = (((xx/3510.0 + 1.0/90.0)*xx + 0.2)*xx + 1.0)*rac/
                     0.906402477055477;
      Res = Math.PI*(Res - temp)/RAC2;
      return Res;
   }


   /**
    * Returns the value of 
    * <SPAN CLASS="MATH"><I>e</I><SUP>x</SUP><I>K</I><SUB>1</SUB>(<I>y</I>)</SPAN>, where <SPAN CLASS="MATH"><I>K</I><SUB>1</SUB></SPAN> is the modified Bessel
    * function of the second kind of order 1. Restriction: <SPAN CLASS="MATH"><I>y</I> &gt; 0</SPAN>.
    * 
    */
//   public static double expBesselK1 (double x, double y) {
//      if (y > 500.0) {
//           double sum = 1 + 3.0/(8.0*y) - 15.0 / (128.0*y*y) + 105.0 /(1024.0*y*y*y);
//           return Math.sqrt(PIsur2/ y) * sum * Math.exp(x - y);
//      } else if (Math.abs(x) > 500.0) {
//         double b = Bessel.k1(y);
//         return Math.exp(x + Math.log(b));
//      } else {
//         return Math.exp(x) * Bessel.k1(y);
//      }
//   }


   /**
    * Returns the value of erf(<SPAN CLASS="MATH"><I>x</I></SPAN>), the error function. It is defined as
    * <P>
    * </P>
    * <DIV ALIGN="CENTER" CLASS="mathdisplay">
    * erf(<I>x</I>) = 2/[(&pi;)<SUP>1/2</SUP>]&int;<SUB>0</SUB><SUP>x</SUP><I>dt</I>&nbsp;<I>e</I><SUP>-t<SUP>2</SUP></SUP>.
    * </DIV><P></P>
    * @param x value at which the function is calculated
    * 
    *        @return the value of erf<SPAN CLASS="MATH">(<I>x</I>)</SPAN>
    * 
    */
   public static double erf (double x) {
      if (x < 0.0)
         return -erf(-x);
      if (x >= 6.0)
         return 1.0;
      if (x >= 2.0)
         return 1.0 - erfc(x);

      double t = 0.5*x*x - 1.0;
      double y = Num.evalCheby (AERF, 16, t);
      return x*y;
   }


   /**
    * Returns the value of erfc(<SPAN CLASS="MATH"><I>x</I></SPAN>), the complementary error function.
    * It is defined as
    * <P>
    * </P>
    * <DIV ALIGN="CENTER" CLASS="mathdisplay">
    * erfc(<I>x</I>) = 2/[(&pi;)<SUP>1/2</SUP>]&int;<SUB>x</SUB><SUP>&#8734;</SUP><I>dt</I>&nbsp;<I>e</I><SUP>-t<SUP>2</SUP></SUP>.
    * </DIV><P></P>
    * 
    * @param x value at which the function is calculated
    * 
    *        @return the value of erfc<SPAN CLASS="MATH">(<I>x</I>)</SPAN>
    * 
    */
   public static double erfc (double x) {
      if (x < 0.0)
         return 2.0 - erfc(-x);
      if (x >= XBIG)
         return 0.0;
      double t = (x - 3.75)/(x + 3.75);
      double y = Num.evalCheby (AERFC, 24, t);
      y *= Math.exp(-x*x);
      return y;
   }


    private static final double[] InvP1 = {
        0.160304955844066229311e2,
       -0.90784959262960326650e2,
        0.18644914861620987391e3,
       -0.16900142734642382420e3,
        0.6545466284794487048e2,
       -0.864213011587247794e1,
        0.1760587821390590
    };

    private static final double[] InvQ1 = {
        0.147806470715138316110e2,
       -0.91374167024260313396e2,
        0.21015790486205317714e3,
       -0.22210254121855132366e3,
        0.10760453916055123830e3,
       -0.206010730328265443e2,
        0.1e1
    };

    private static final double[] InvP2 = {
       -0.152389263440726128e-1,
        0.3444556924136125216,
       -0.29344398672542478687e1,
        0.11763505705217827302e2,
       -0.22655292823101104193e2,
        0.19121334396580330163e2,
       -0.5478927619598318769e1,
        0.237516689024448000
    };

    private static final double[] InvQ2 = {
      -0.108465169602059954e-1,
       0.2610628885843078511,
      -0.24068318104393757995e1,
       0.10695129973387014469e2,
      -0.23716715521596581025e2,
       0.24640158943917284883e2,
      -0.10014376349783070835e2,
       0.1e1
    };

    private static final double[] InvP3 = {
        0.56451977709864482298e-4,
        0.53504147487893013765e-2,
        0.12969550099727352403,
        0.10426158549298266122e1,
        0.28302677901754489974e1,
        0.26255672879448072726e1,
        0.20789742630174917228e1,
        0.72718806231556811306,
        0.66816807711804989575e-1,
       -0.17791004575111759979e-1,
        0.22419563223346345828e-2
    };

    private static final double[] InvQ3 = {
        0.56451699862760651514e-4,
        0.53505587067930653953e-2,
        0.12986615416911646934,
        0.10542932232626491195e1,
        0.30379331173522206237e1,
        0.37631168536405028901e1,
        0.38782858277042011263e1,
        0.20372431817412177929e1,
        0.1e1
    };

   /**
    * Returns the value of 
    * <SPAN CLASS="MATH"><I>y</I> = erf<SUP>-1</SUP>(<I>x</I>)</SPAN>, the inverse  of  the error
    * function such that  
    * <SPAN CLASS="MATH"><I>x</I> = erf(<I>y</I>)</SPAN>.
    * 
    * @param x value at which the function is calculated
    * 
    *        @return the value of erfInv<SPAN CLASS="MATH">(<I>x</I>)</SPAN>
    * 
    */
   public static double erfInv (double x)  {
      if (x < 0.0)
         return -erfInv (-x);

      if (x > 1.0)
         throw new IllegalArgumentException ("x is not in [-1, 1]");
      if (x >= 1.0)
         return Double.POSITIVE_INFINITY;

      int i;
      double y, z, v, w;

      if (x <= 0.75) {
         y = x * x - 0.5625;
         v = Misc.evalPoly (InvP1, 6, y);
         w = Misc.evalPoly (InvQ1, 6, y);
         z = (v / w) * x;

      } else if (x <= 0.9375) {
         y = x * x - 0.87890625;
         v = Misc.evalPoly (InvP2, 7, y);
         w = Misc.evalPoly (InvQ2, 7, y);
         z = (v / w) * x;

      } else {
         y = 1.0 / Math.sqrt (-Math.log (1.0 - x));
         v = Misc.evalPoly (InvP3, 10, y);
         w = Misc.evalPoly (InvQ3, 8, y);
         z = (v / w) / y;
      }

      return z;
   }

}
