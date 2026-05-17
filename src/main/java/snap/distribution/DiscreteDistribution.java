
package snap.distribution;

import java.util.Formatter;
import java.util.Locale;



/**
 * Classes implementing discrete distributions over a 
 * <SPAN  CLASS="textit">finite set of real numbers</SPAN> should inherit from this class. 
 * For discrete distributions over integers, 
 * see {@link snap.distribution.DiscreteDistributionInt DiscreteDistributionInt}.
 * 
 * <P>
 * We assume that the random variable <SPAN CLASS="MATH"><I>X</I></SPAN> of interest can take one of the
 * <SPAN CLASS="MATH"><I>n</I></SPAN> values 
 * <SPAN CLASS="MATH"><I>x</I><SUB>0</SUB> &lt; <SUP> ... </SUP> &lt; <I>x</I><SUB>n-1</SUB></SPAN> (which are <SPAN  CLASS="textit">sorted</SPAN> by increasing order).
 * It takes the value <SPAN CLASS="MATH"><I>x</I><SUB>k</SUB></SPAN> with probability 
 * <SPAN CLASS="MATH"><I>p</I><SUB>k</SUB> = <I>P</I>[<I>X</I> = <I>x</I><SUB>k</SUB>]</SPAN>.
 * In addition to the methods specified in the interface 
 * {@link snap.distribution.Distribution Distribution},
 * a method that returns the probability <SPAN CLASS="MATH"><I>p</I><SUB>k</SUB></SPAN> is supplied. 
 * 
 * <P>
 * Note that the default implementation of the complementary distribution function
 * returns <TT>1.0 - cdf(x - 1)</TT>, which is not accurate when <SPAN CLASS="MATH"><I>F</I>(<I>x</I>)</SPAN> is near 1.
 * 
 */
public class DiscreteDistribution implements Distribution {
  /*
     For better precision in the tails, we keep the cumulative probabilities
     (F) in cdf[x] for x <= xmed (i.e. cdf[x] is the sum off all the probabi-
     lities pr[i] for i <= x),
     and the complementary cumulative probabilities (1 - F) in cdf[x] for
     x > xmed (i.e. cdf[x] is the sum off all the probabilities pr[i]
     for i >= x).
  */ 
  
   protected double cdf[] = null;    // cumulative probabilities
   protected double pr[] = null;     // probability terms or mass distribution
   protected int xmin = 0;           // pr[x] = 0 for x < xmin
   protected int xmax = 0;           // pr[x] = 0 for x > xmax
   protected int xmed = 0;           // cdf[x] = F (x) for x <= xmed, and 
                                     // cdf[x] = bar_F (x) for x > xmed
   protected int    nObs;            // number of different observations
   protected int    nObsTot;         // total number of observations
   protected double obs[];
   protected double sortedObs[];
   protected double supportA = Double.NEGATIVE_INFINITY;
   protected double supportB = Double.POSITIVE_INFINITY;



   protected DiscreteDistribution () {} 
   // Default constructor called by subclasses such as 'EmpiricalDist'


   /**
    * Constructs a discrete distribution over the <SPAN CLASS="MATH"><I>n</I></SPAN> values
    *  contained in array <TT>obs</TT>, with probabilities given in array <TT>prob</TT>.
    *  Both arrays must have at least <SPAN CLASS="MATH"><I>n</I></SPAN> elements, the probabilities must
    *  sum to 1, and the observations are assumed to be sorted by increasing order.
    * 
    */
   public DiscreteDistribution (double[] obs, double[] prob, int n) {
      init(n, obs, prob);
   }


   /**
    * Constructs a discrete distribution whose parameters are given 
    *    in a single ordered array: <TT>params[0]</TT> contains <SPAN CLASS="MATH"><I>n</I></SPAN>, the number of 
    *    values to consider. Then the next <SPAN CLASS="MATH"><I>n</I></SPAN> values of <TT>params</TT> are the
    *     observation values, and the last <SPAN CLASS="MATH"><I>n</I></SPAN> values of <TT>params</TT>
    *    are the probabilities values.
    * 
    * 
    */
   public DiscreteDistribution (double[] params) {

      if (params.length != 1+params[0]*2)
         throw new IllegalArgumentException("Wrong parameter size");

      int n =  (int)params[0];
      double[] obs = new double[n];
      double[] prob = new double[n];

      //int indice = 1;
      System.arraycopy (params, 1, obs, 0, n);
      System.arraycopy (params, n+1, prob, 0, n);
      /*for(int i=0; i<no; i++) {
         obs[i] = params[indice];
         indice++;
      }
      
      for(int i=0; i<np; i++) {
         prob[i] = params[indice];
         indice++;
      }*/
      init(n, obs, prob);
    }


   private void init(int n, double[] obs, double[] prob) {

     int no = obs.length;
     int np = prob.length;
     if (n <= 0) throw new IllegalArgumentException ("n <= 0");
     if (no < (n-1) || np < (n-1))
        throw new IllegalArgumentException 
                  ("Size of arrays 'obs' or 'prob' less than 'n'");

     nObs     = n;
     this.obs = obs;
     pr      = prob;

     // cdf
     sortedObs = new double[nObs];
     System.arraycopy (obs, 0, sortedObs, 0, nObs);

     supportA = sortedObs[0];
     supportB = sortedObs[nObs-1];

     cdf = new double[nObs];
     cdf[0] = pr[0];
     cdf[nObs - 1] = pr[nObs - 1];
     xmin = 0;
     xmax = nObs-1;

    /* Compute the cumulative probabilities until F >= 0.5, and keep them in
       the lower part of cdf */
    int i = 0;
    while (i < xmax && cdf[i] < 0.5) {
       i++;
       cdf[i] = pr[i] + cdf[i - 1];
    }
    // This is the boundary between F and barF in the CDF
    xmed = i;

    /* Compute the cumulative probabilities of the complementary distribution
       and keep them in the upper part of cdf. */
    i = nObs - 2;
    while (i > xmed) {
       cdf[i] = pr[i] + cdf[i + 1];
       i--;
    }
}

   /**
    * Returns <SPAN CLASS="MATH"><I>p</I><SUB>k</SUB></SPAN>, the probability of 
    *   the <SPAN CLASS="MATH"><I>k</I></SPAN>-th observation, for <SPAN CLASS="MATH">0&nbsp;&lt;=&nbsp;<I>k</I> &lt; <I>n</I></SPAN>.
    *   The result should be a real number in the interval <SPAN CLASS="MATH">[0, 1]</SPAN>.
    * 
    * @param k observation number, 
    * <SPAN CLASS="MATH">0&nbsp;&lt;=&nbsp;<I>k</I> &lt; <I>n</I></SPAN>
    * 
    *    @return the probability of observation <TT>k</TT>
    * 
    * 
    */
   public double prob (int k) {
      return pr[k];
   }


   /**
    * @param x value at which the distribution function must be evaluated
    * 
    *    @return the distribution function evaluated at <TT>x</TT>
    * 
    */
   public double cdf (double x) {
      if (x < sortedObs[0])
         return 0.0;
      if (x >= sortedObs[nObs-1])
         return 1.0;
      if ((xmax == xmed) || (x < sortedObs[xmed+1])) {
         for (int i = 0; i <= xmed; i++)
            if (x >= sortedObs[i] && x < sortedObs[i+1])
               return cdf[i];
      } else {
         for (int i = xmed + 1; i < nObs-1; i++)
            if (x >= sortedObs[i] && x < sortedObs[i+1])
               return 1.0 - cdf[i+1];
      }
      throw new IllegalStateException();

   }


   /**
    * @param x value at which the complementary distribution function must be evaluated
    * 
    *    @return the complementary distribution function evaluated at <TT>x</TT>
    * 
    */
   public double barF (double x) {
      if (x <= sortedObs[0])
         return 1.0;
      if (x > sortedObs[nObs-1])
         return 0.0;
      if ((xmax == xmed) || (x <= sortedObs[xmed+1])) {
         for (int i = 0; i <= xmed; i++)
            if (x > sortedObs[i] && x <= sortedObs[i+1])
               return 1.0 - cdf[i];
      } else {
         for (int i = xmed + 1; i < nObs-1; i++)
            if (x > sortedObs[i] && x <= sortedObs[i+1])
               return cdf[i + 1];
      }
      throw new IllegalStateException();
   }


   public double inverseF (double u) {

      int i, j, k;

      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u not in [0,1]");
      if (u <= 0.0)
         return supportA;
      if (u >= 1.0)
         return supportB;


      // Remember: the upper part of cdf contains the complementary distribu-
      // tion for xmed < s <= xmax, and the lower part of cdf the
      // distribution for xmin <= s <= xmed 

      if (u <= cdf[xmed - xmin]) {
         // In the lower part of cdf
         if (u <= cdf[0])
            return sortedObs[xmin];
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
            return sortedObs[xmax];

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

      return sortedObs[i + xmin];
   }


   /**
    * Computes the mean 
    * <SPAN CLASS="MATH"><I>E</I>[<I>X</I>] = &sum;<SUB>i=1</SUB><SUP>n</SUP><I>p</I><SUB>i</SUB><I>x</I><SUB>i</SUB></SPAN> of the distribution.
    * 
    */
   public double getMean() {
      double mean = 0.0;

      for (int i = 0; i < nObs; i++)
         mean += obs[i] * pr[i];

      return mean;
   }


   /**
    * Computes the variance 
    * <SPAN CLASS="MATH">Var[<I>X</I>] = &sum;<SUB>i=1</SUB><SUP>n</SUP><I>p</I><SUB>i</SUB>(<I>x</I><SUB>i</SUB> - <I>E</I>[<I>X</I>])<SUP>2</SUP></SPAN>
    *    of the distribution.
    * 
    */
   public double getVariance() {
      double variance = 0.0;
      double mean = getMean();

      for (int i = 0; i < nObs; i++)
         variance += (obs[i] - mean) * (obs[i] - mean) * pr[i];

      return (variance / (double) nObs);
   }


   /**
    * Computes the standard deviation of the distribution.
    * 
    */
   public double getStandardDeviation() {
      return Math.sqrt (getVariance());
   }


   /**
    * Returns a table containing the parameters of the current distribution.
    *    This table is built in regular order, according to constructor
    *    <TT>DiscreteDistribution(double[] params)</TT> order.
    * 
    */
   public double[] getParams() {
      double[] retour = new double[1+nObs*2];
      double sum = 0;
      retour[0] = nObs;
      System.arraycopy (obs, 0, retour, 1, nObs);
      for(int i = 0; i<nObs-1; i++) {
         retour[nObs+1+i] = cdf[i] - sum;
         sum = cdf[i];
      }
      retour[2*nObs] = 1.0 - sum;

      return retour;
   }


   /**
    * Returns the lower limit <SPAN CLASS="MATH"><I>x</I><SUB>a</SUB></SPAN> of the support of the distribution.
    *  The probability is 0 for all <SPAN CLASS="MATH"><I>x</I> &lt; <I>x</I><SUB>a</SUB></SPAN>. 
    * 
    * @return <SPAN CLASS="MATH"><I>x</I></SPAN> lower limit of support
    * 
    */
   public double getXinf() {
      return supportA;
   }


   /**
    * Returns the upper limit <SPAN CLASS="MATH"><I>x</I><SUB>b</SUB></SPAN> of the support of the distribution.
    *  The probability is 0 for all <SPAN CLASS="MATH"><I>x</I> &gt; <I>x</I><SUB>b</SUB></SPAN>.
    * 
    * @return <SPAN CLASS="MATH"><I>x</I></SPAN> upper limit of support
    * 
    */
   public double getXsup() {
      return supportB;
   }


   /**
    * Returns a <TT>String</TT> containing information about the current distribution.
    * 
    */
   public String toString() {
      System.out.println(cdf + " : " + obs);
      StringBuilder sb = new StringBuilder();
      Formatter formatter = new Formatter(sb, Locale.US);
      formatter.format("%s%n", getClass().getSimpleName());
      formatter.format("%s :  %s%n", "value", "cdf");
      int i;
      for(i = 0; i<nObs-1; i++) {
         formatter.format("%f : %f%n", obs[i], cdf[i]);
      }
         formatter.format("%f : %f%n", obs[i], 1.0);
      return sb.toString();
   }

}
