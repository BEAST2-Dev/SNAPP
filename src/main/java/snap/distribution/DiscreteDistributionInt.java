
package snap.distribution;


/**
 * Classes implementing discrete distributions over the integers should
 * inherit from this class.
 * It specifies the signatures of methods for computing the mass function
 * (or probability) 
 * <SPAN CLASS="MATH"><I>p</I>(<I>x</I>) = <I>P</I>[<I>X</I> = <I>x</I>]</SPAN>, distribution function <SPAN CLASS="MATH"><I>F</I>(<I>x</I>)</SPAN>, 
 * complementary distribution function <SPAN CLASS="MATH">bar(F)(<I>x</I>)</SPAN>,
 * and inverse distribution function <SPAN CLASS="MATH"><I>F</I><SUP>-1</SUP>(<I>u</I>)</SPAN>, for a random variable <SPAN CLASS="MATH"><I>X</I></SPAN> 
 * with a discrete distribution over the integers.
 * 
 * <P>
 * <SPAN  CLASS="textit">WARNING:</SPAN> the complementary distribution function is defined as 
 * 
 * <SPAN CLASS="MATH">bar(F)(<I>j</I>) = <I>P</I>[<I>X</I>&nbsp;&gt;=&nbsp;<I>j</I>]</SPAN> (for integers <SPAN CLASS="MATH"><I>j</I></SPAN>, so that for discrete distributions
 * in SSJ, 
 * <SPAN CLASS="MATH"><I>F</I>(<I>j</I>) + bar(F)(<I>j</I>)&#8800;1</SPAN> since both include the term <SPAN CLASS="MATH"><I>P</I>[<I>X</I> = <I>j</I>]</SPAN>.
 * 
 * <P>
 * The implementing classes provide both static and non-static methods
 * to compute the above functions. 
 * The non-static methods require the creation of an object of
 * class {@link snap.distribution.DiscreteDistributionInt DiscreteDistributionInt}; 
 * all the non-negligible terms of the mass and distribution functions will be
 * precomputed by the constructor and kept in arrays. Subsequent accesses
 * will be very fast.
 * The static methods do not require the construction of an object.
 * These static methods are not specified in this abstract class because
 * the number and types of their parameters depend on the distribution.
 * When methods have to be called several times
 * with the same parameters  for the distributions, 
 * it is usually more efficient to create an object and use its non-static
 * methods instead of the static ones.
 * This trades memory for speed.
 * 
 */
public abstract class DiscreteDistributionInt implements Distribution {

   /**
    * Environment variable that determines what probability terms can
    *   be considered as negligible when building precomputed tables for
    *   distribution and mass functions.  Probabilities smaller than <TT>EPSILON</TT>
    *   are not stored in the 
    *   {@link snap.distribution.DiscreteDistribution DiscreteDistribution} objects 
    *   (such as those of class {@link PoissonDist}, etc.), but are computed 
    *   directly each time they are needed (which should be very seldom).
    *   The default value is set to <SPAN CLASS="MATH">10<SUP>-16</SUP></SPAN>.
    * 
    */
   public static double EPSILON = 1.0e-16;
  /*
     For better precision in the tails, we keep the cumulative probabilities
     (F) in cdf[x] for x <= xmed (i.e. cdf[x] is the sum off all the probabi-
     lities pdf[i] for i <= x),
     and the complementary cumulative probabilities (1 - F) in cdf[x] for
     x > xmed (i.e. cdf[x] is the sum off all the probabilities pdf[i]
     for i >= x).
  */ 
   protected final static double EPS_EXTRA = 1 / 100.0;
   protected double cdf[] = null;    // cumulative probabilities
   protected double pdf[] = null;    // probability terms or mass distribution
   protected int xmin = 0;           // pdf[x] = 0 for x < xmin
   protected int xmax = 0;           // pdf[x] = 0 for x > xmax
   protected int xmed = 0;           // cdf[x] = F (x) for x <= xmed, and 
                                        // cdf[x] = bar_F (x) for x > xmed
   protected int supportA = Integer.MIN_VALUE;
   protected int supportB = Integer.MAX_VALUE;


   /**
    * Returns <SPAN CLASS="MATH"><I>p</I>(<I>x</I>)</SPAN>, the probability of <SPAN CLASS="MATH"><I>x</I></SPAN>,
    *    which should be a real number in the interval <SPAN CLASS="MATH">[0, 1]</SPAN>.
    * 
    * @param x value at which the mass function must be evaluated
    * 
    *    @return the mass function evaluated at <TT>x</TT>
    * 
    */
   public abstract double prob (int x);


   /**
    * Returns the distribution function <SPAN CLASS="MATH"><I>F</I></SPAN> evaluated at <TT>x</TT>
    * (see).
    *   Calls the {@link #cdf(int) cdf}<TT>(int)</TT> method.
    * 
    * @param x value at which the distribution function must be evaluated
    * 
    *    @return the distribution function evaluated at <TT>x</TT>
    * 
    */
   public double cdf (double x) {
     return cdf ((int) x);
   }


   /**
    * Returns the distribution function <SPAN CLASS="MATH"><I>F</I></SPAN> evaluated at <TT>x</TT>
    * (see).
    * 
    * @param x value at which the distribution function must be evaluated
    * 
    *    @return the distribution function evaluated at <TT>x</TT>
    * 
    */
   public abstract double cdf (int x);


   /**
    * Returns <SPAN CLASS="MATH">bar(F)(<I>x</I>)</SPAN>, the complementary distribution function.
    *   Calls the {@link #barF(int) barF}<TT>(int)</TT> method.
    * 
    * @param x value at which the complementary distribution function
    *     must be evaluated
    * 
    *    @return the complementary distribution function evaluated at <TT>x</TT>
    * 
    */
   public double barF (double x) {
      return barF ((int) x);
   }

   
   /**
    * Returns <SPAN CLASS="MATH">bar(F)(<I>x</I>)</SPAN>, the complementary 
    *    distribution function. <SPAN  CLASS="textit">See the WARNING above</SPAN>.
    * 
    * @param x value at which the complementary distribution function
    *     must be evaluated
    * 
    *    @return the complementary distribution function evaluated at <TT>x</TT>
    * 
    */
   public double barF (int x) {
      return 1.0 - cdf (x - 1);
   }


   /**
    * Returns the lower limit <SPAN CLASS="MATH"><I>x</I><SUB>a</SUB></SPAN> of the support of the probability
    *  mass function. The probability is 0 for all <SPAN CLASS="MATH"><I>x</I> &lt; <I>x</I><SUB>a</SUB></SPAN>.
    * 
    * @return <SPAN CLASS="MATH"><I>x</I></SPAN> lower limit of support
    * 
    */
   public int getXinf() {
      return supportA;
   }


   /**
    * Returns the upper limit <SPAN CLASS="MATH"><I>x</I><SUB>b</SUB></SPAN> of the support of the  probability
    *  mass function. The probability is 0 for all <SPAN CLASS="MATH"><I>x</I> &gt; <I>x</I><SUB>b</SUB></SPAN>.
    * 
    * @return <SPAN CLASS="MATH"><I>x</I></SPAN> upper limit of support
    * 
    */
   public int getXsup() {
      return supportB;
   }


   /**
    * Returns the inverse distribution function
    *   <SPAN CLASS="MATH"><I>F</I><SUP>-1</SUP>(<I>u</I>)</SPAN>, where 
    * <SPAN CLASS="MATH">0&nbsp;&lt;=&nbsp;<I>u</I>&nbsp;&lt;=&nbsp;1</SPAN>. Calls the <TT>inverseFInt</TT> method.
    * 
    * @param u value in the interval <SPAN CLASS="MATH">(0, 1)</SPAN> for which
    *              the inverse distribution function is evaluated
    * 
    *    @return the inverse distribution function evaluated at <TT>u</TT>
    *    @exception IllegalArgumentException if <SPAN CLASS="MATH"><I>u</I></SPAN> is  not in the interval <SPAN CLASS="MATH">(0, 1)</SPAN>
    * 
    *    @exception ArithmeticException if the inverse cannot be computed,
    *      for example if it would give infinity in a theoritical context
    * 
    * 
    */
   public double inverseF (double u) {
      return inverseFInt (u);
   }


   /**
    * Returns the inverse distribution function
    *   <SPAN CLASS="MATH"><I>F</I><SUP>-1</SUP>(<I>u</I>)</SPAN>, where 
    * <SPAN CLASS="MATH">0&nbsp;&lt;=&nbsp;<I>u</I>&nbsp;&lt;=&nbsp;1</SPAN>.
    *   The default implementation uses binary search.
    * 
    * @param u value in the interval <SPAN CLASS="MATH">(0, 1)</SPAN> for which
    *              the inverse distribution function is evaluated
    * 
    *    @return the inverse distribution function evaluated at <TT>u</TT>
    *    @exception IllegalArgumentException if <SPAN CLASS="MATH"><I>u</I></SPAN> is  not in the interval <SPAN CLASS="MATH">(0, 1)</SPAN>
    * 
    *    @exception ArithmeticException if the inverse cannot be computed,
    *      for example if it would give infinity in a theoritical context
    * 
    * 
    */
   public int inverseFInt (double u) {
      int i, j, k;

      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u is not in [0,1]");
      if (u <= 0.0)
         return supportA;
      if (u >= 1.0)
         return supportB;


      // Remember: the upper part of cdf contains the complementary distribu-
      // tion for xmed < s <= xmax, and the lower part of cdf the
      // distribution for xmin <= x <= xmed 

      if (u <= cdf[xmed - xmin]) {
         // In the lower part of cdf
         if (u <= cdf[0])
            return xmin;
         i = 0;
         j = xmed - xmin;
         while (i < j) {
            k = (i + j) / 2;
            if (u > cdf[k])
               i = k + 1;
            else
               j = k;
         }
      }
      else {
         // In the upper part of cdf
         u = 1 - u;
         if (u < cdf[xmax - xmin])
            return xmax;

         i = xmed - xmin + 1;
         j = xmax - xmin;
         while (i < j) {
            k = (i + j) / 2;
            if (u < cdf[k])
               i = k + 1;
            else
               j = k;
         }
         i--;
      }

      return i + xmin;
   }

}
