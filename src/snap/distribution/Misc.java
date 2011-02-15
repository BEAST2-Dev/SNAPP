
package snap.distribution;


/**
 * This class provides miscellaneous functions that are hard to classify.
 * Some may be moved to another class in the future.
 * 
 */
public class Misc {
   private Misc() {}


   /**
    * Returns the 
    * <SPAN CLASS="MATH"><I>k</I><SUP>th</SUP></SPAN> smallest item of the array <SPAN CLASS="MATH"><I>t</I></SPAN> of size <SPAN CLASS="MATH"><I>n</I></SPAN>.
    *    Array <SPAN CLASS="MATH"><I>t</I></SPAN> is unchanged by the method.
    * 
    * @param t the array which contain the items
    * 
    *    @param n the number of items in the array
    * 
    *    @param k the index of the smallest item
    * 
    *    @return the kth smallest item of the array <SPAN CLASS="MATH"><I>t</I></SPAN>
    * 
    */
   public static double quickSelect (double[] t, int n, int k) {
      double[] U = new double[n];
      double[] V = new double[n];
      double p = t[k - 1];
      int u = 0;
      int v = 0;
      int indV = 0;

      for (int i = 0; i < n; i++) {
         if (t[i] <= p) {
            v++;
            if (t[i] != p) {
               U[u++] = t[i];
            }
         } else
            V[indV++] = t[i];
      }

      if (k <= u)
         return quickSelect (U, u, k);
      else if (k > v)
         return quickSelect (V, indV, k - v);
      else return p;
   }


   /**
    * Returns the 
    * <SPAN CLASS="MATH"><I>k</I><SUP>th</SUP></SPAN> smallest item of the array <SPAN CLASS="MATH"><I>t</I></SPAN> of size <SPAN CLASS="MATH"><I>n</I></SPAN>.
    *    Array <SPAN CLASS="MATH"><I>t</I></SPAN> is unchanged by the method.
    * 
    * @param t the array which contain the items
    * 
    *    @param n the number of items in the array
    * 
    *    @param k the index of the smallest item
    * 
    *    @return the kth smallest item of the array <SPAN CLASS="MATH"><I>t</I></SPAN>
    * 
    */
   public static int quickSelect (int[] t, int n, int k) {
      int[] U = new int[n];
      int[] V = new int[n];
      int p = t[k - 1];
      int u = 0;
      int v = 0;
      int indV = 0;

      for (int i = 0; i < n; i++) {
         if (t[i] <= p) {
            v++;
            if (t[i] != p) {
               U[u++] = t[i];
            }
         } else
            V[indV++] = t[i];
      }

      if (k <= u)
         return quickSelect (U, u, k);
      else if (k > v)
         return quickSelect (V, indV, k - v);
      else return p;
   }


   /**
    * Returns the index of the time interval corresponding to time <TT>t</TT>.
    *  Let 
    * <SPAN CLASS="MATH"><I>t</I><SUB>0</SUB>&nbsp;&lt;=&nbsp;<SUP> ... </SUP>&nbsp;&lt;=&nbsp;<I>t</I><SUB>n</SUB></SPAN> be simulation times stored in a subset of
    *  <TT>times</TT>.  This method uses binary search to determine the
    *  smallest value <SPAN CLASS="MATH"><I>i</I></SPAN> for which 
    * <SPAN CLASS="MATH"><I>t</I><SUB>i</SUB>&nbsp;&lt;=&nbsp;<I>t</I> &lt; <I>t</I><SUB>i+1</SUB></SPAN>, and returns <SPAN CLASS="MATH"><I>i</I></SPAN>.
    *  The value of <SPAN CLASS="MATH"><I>t</I><SUB>i</SUB></SPAN> is stored in <TT>times[start+i]</TT> whereas
    *  <SPAN CLASS="MATH"><I>n</I></SPAN> is defined as <TT>end - start</TT>.
    *  If <SPAN CLASS="MATH"><I>t</I> &lt; <I>t</I><SUB>0</SUB></SPAN>, this returns <SPAN CLASS="MATH">-1</SPAN>.  If <SPAN CLASS="MATH"><I>t</I>&nbsp;&gt;=&nbsp;<I>t</I><SUB>n</SUB></SPAN>, this returns <SPAN CLASS="MATH"><I>n</I></SPAN>.
    *  Otherwise, the returned value is greater than or equal to 0, and
    *  smaller than or equal to <SPAN CLASS="MATH"><I>n</I> - 1</SPAN>. <TT>start</TT> and <TT>end</TT> are only used
    *  to set lower and upper limits of the search in the <TT>times</TT>
    *  array; the index space of the returned value always starts at 0.
    *  Note that if the elements of <TT>times</TT> with indices <TT>start</TT>,
    *  ..., <TT>end</TT> are not sorted in non-decreasing order,
    *  the behavior of this method is undefined.
    * 
    * @param times an array of simulation times.
    * 
    *    @param start the first index in the array to consider.
    * 
    *    @param end the last index (inclusive) in the array to consider.
    * 
    *    @param t the queried simulation time.
    * 
    *    @return the index of the interval.
    *    @exception NullPointerException if <TT>times</TT> is <TT>null</TT>.
    * 
    *    @exception IllegalArgumentException if <TT>start</TT> is negative,
    *     or if <TT>end</TT> is smaller than <TT>start</TT>.
    * 
    *    @exception ArrayIndexOutOfBoundsException if <TT>start + end</TT>
    *     is greater than or equal to the length of <TT>times</TT>.
    * 
    * 
    */
   public static int getTimeInterval (double[] times, int start, int end,
                                      double t) {
      if (start < 0)
         throw new IllegalArgumentException
            ("The starting index must not be negative");
      int n = end - start;
      if (n < 0)
         throw new IllegalArgumentException
            ("The ending index must be greater than or equal to the starting index");
      if (t < times[start])
         return -1;
      if (t >= times[end])
         return n;

      int start0 = start;
      // Perform binary search to find the interval index
      int mid = (start + end)/2;
      // Test if t is inside the interval mid.
      // The interval mid starts at times[mid],
      // and the interval mid+1 starts at times[mid + 1].
      while (t < times[mid] || t >= times[mid + 1]) {
         if (start == end)
            // Should not happen, safety check to avoid infinite loops.
            throw new IllegalStateException();
         if (t < times[mid])
            // time corresponds to an interval before mid.
            end = mid - 1;
         else
            // time corresponds to an interval after mid.
            start = mid + 1;
         mid = (start + end)/2;
      }
      return mid - start0;
   }


   /**
    * Computes the Newton interpolating polynomial.  Given the <SPAN CLASS="MATH"><I>n</I> + 1</SPAN>
    *  real distinct points 
    * <SPAN CLASS="MATH">(<I>x</I><SUB>0</SUB>, <I>y</I><SUB>0</SUB>),</SPAN> 
    * <SPAN CLASS="MATH">(<I>x</I><SUB>1</SUB>, <I>y</I><SUB>1</SUB>),&#8230;,(<I>x</I><SUB>n</SUB>, <I>y</I><SUB>n</SUB>)</SPAN>,
    *   with <TT>X[i]</TT> <SPAN CLASS="MATH">= <I>x</I><SUB>i</SUB></SPAN>, <TT>Y[i]</TT> <SPAN CLASS="MATH">= <I>y</I><SUB>i</SUB></SPAN>, this function computes
    *   the <SPAN CLASS="MATH"><I>n</I> + 1</SPAN> coefficients <TT>C[i]</TT> <SPAN CLASS="MATH">= <I>c</I><SUB>i</SUB></SPAN> of the Newton
    *   interpolating polynomial <SPAN CLASS="MATH"><I>P</I>(<I>x</I>)</SPAN> of degree <SPAN CLASS="MATH"><I>n</I></SPAN> passing through these points,
    *   i.e. such that 
    * <SPAN CLASS="MATH"><I>y</I><SUB>i</SUB> = <I>P</I>(<I>x</I><SUB>i</SUB>)</SPAN>, given by
    * <P><A NAME="eq.newton.interpol"></A>
    * </P>
    * <DIV ALIGN="CENTER" CLASS="mathdisplay"><A NAME="eq.newton.interpol"></A>
    * &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<I>P</I>(<I>x</I>) = <I>c</I><SUB>0</SUB> + <I>c</I><SUB>1</SUB>(<I>x</I> - <I>x</I><SUB>0</SUB>) + <I>c</I><SUB>2</SUB>(<I>x</I> - <I>x</I><SUB>0</SUB>)(<I>x</I> - <I>x</I><SUB>1</SUB>) + <SUP> ... </SUP> + <I>c</I><SUB>n</SUB>(<I>x</I> - <I>x</I><SUB>0</SUB>)(<I>x</I> - <I>x</I><SUB>1</SUB>)<SUP> ... </SUP>(<I>x</I> - <I>x</I><SUB>n-1</SUB>).
    * </DIV><P></P>
    * 
    * @param n degree of the interpolating polynomial
    * 
    *    @param X <SPAN CLASS="MATH"><I>x</I></SPAN>-coordinates of points
    * 
    *    @param Y <SPAN CLASS="MATH"><I>y</I></SPAN>-coordinates of points
    * 
    *    @param C Coefficients of the interpolating polynomial
    * 
    * 
    */
   public static void interpol (int n, double[] X, double[] Y, double[] C) {
      int j;
      // Compute divided differences for the Newton interpolating polynomial
      for (j = 0; j <= n; ++j)
         C[j] = Y[j];
      for (int i = 1; i <= n; ++i)
         for (j = n; j >= i; --j) {
            if (X[j] == X[j-i])
               C[j] = 0;
            else
               C[j] = (C[j] - C[j-1]) / (X[j] - X[j-i]);
         }
   }


   /**
    * Given <SPAN CLASS="MATH"><I>n</I></SPAN>, <SPAN CLASS="MATH"><I>X</I></SPAN> and <SPAN CLASS="MATH"><I>C</I></SPAN> as described in
    *  {@link #interpol(int,double[],double[],double[]) interpol}<TT>(n, X, Y, C)</TT>, this
    * function returns the value of the interpolating polynomial <SPAN CLASS="MATH"><I>P</I>(<I>z</I>)</SPAN> evaluated
    *  at <SPAN CLASS="MATH"><I>z</I></SPAN> (see eq. ).
    * 
    * @param n degree of the interpolating polynomial
    * 
    *    @param X <SPAN CLASS="MATH"><I>x</I></SPAN>-coordinates of points
    * 
    *    @param C Coefficients of the interpolating polynomial
    * 
    *    @param z argument where polynomial is evaluated
    * 
    *    @return Value of the interpolating polynomial <SPAN CLASS="MATH"><I>P</I>(<I>z</I>)</SPAN>
    * 
    */
   public static double evalPoly (int n, double[] X, double[] C, double z)  {
      double v = C[n];
      for (int j = n-1; j >= 0; --j)
         v = v*(z - X[j]) + C[j];
      return v;
   }


   /**
    * Evaluates the polynomial 
    * <P><A NAME="eq.horner"></A>
    * </P>
    * <DIV ALIGN="CENTER" CLASS="mathdisplay"><A NAME="eq.horner"></A>
    * &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<I>P</I>(<I>x</I>) = <I>c</I><SUB>0</SUB> + <I>c</I><SUB>1</SUB><I>x</I> + <I>c</I><SUB>2</SUB><I>x</I><SUP>2</SUP> + <SUP> ... </SUP> + <I>c</I><SUB>n</SUB><I>x</I><SUP>n</SUP>
    * </DIV><P></P>
    * of degree <SPAN CLASS="MATH"><I>n</I></SPAN> with coefficients <SPAN CLASS="MATH"><I>c</I><SUB>j</SUB> =</SPAN> <TT>C[j]</TT> at <SPAN CLASS="MATH"><I>x</I></SPAN>. 
    * 
    * @param C Coefficients of the polynomial
    * 
    *    @param n degree of the polynomial
    * 
    *    @param x argument where polynomial is evaluated
    * 
    *    @return Value of the polynomial <SPAN CLASS="MATH"><I>P</I>(<I>x</I>)</SPAN>
    * 
    */
   public static double evalPoly (double[] C, int n, double x)  {
      double v = C[n];
      for (int j = n-1; j >= 0; --j)
         v = v*x + C[j];
      return v;
   }

}
