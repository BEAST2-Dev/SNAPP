
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
 * the <EM>exponential</EM> distribution
 * with mean <SPAN CLASS="MATH">1/<I>&#955;</I></SPAN> where 
 * <SPAN CLASS="MATH"><I>&#955;</I> &gt; 0</SPAN>.
 * Its density is
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay">
 * <I>f</I> (<I>x</I>) = <I>&#955;e</I><SUP>-<I>&#955;</I>x</SUP>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;for <I>x</I>&nbsp;&gt;=&nbsp;0,
 * </DIV><P></P>
 * its distribution function is
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay">
 * <I>F</I>(<I>x</I>) = 1 - <I>e</I><SUP>-<I>&#955;</I>x</SUP>,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;for <I>x</I>&nbsp;&gt;=&nbsp;0,
 * </DIV><P></P>
 * and its inverse distribution function is
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay">
 * <I>F</I><SUP>-1</SUP>(<I>u</I>) = - ln(1 - <I>u</I>)/<I>&#955;</I>,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;for 0 &lt; <I>u</I> &lt; 1.
 * </DIV><P></P>
 * 
 */
public class ExponentialDist extends ContinuousDistribution {
   private double lambda;



   /**
    * Constructs an <TT>ExponentialDist</TT> object with parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN> = 1.
    * 
    */
   public ExponentialDist() {
      setLambda (1.0);
   }


   /**
    * Constructs an <TT>ExponentialDist</TT> object with parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN> =
    *   <TT>lambda</TT>.
    * 
    */
   public ExponentialDist (double lambda) {
      setLambda (lambda);
  }


   public double density (double x) {
      return density (lambda, x);
   }

   public double cdf (double x) {
      return cdf (lambda, x);
   }

   public double barF (double x) {
      return barF (lambda, x);
   }

   public double inverseF (double u) {
      return inverseF (lambda, u);
   }

   public double getMean() {
      return ExponentialDist.getMean (lambda);
   }

   public double getVariance() {
      return ExponentialDist.getVariance (lambda);
   }

   public double getStandardDeviation() {
      return ExponentialDist.getStandardDeviation (lambda);
   }

   /**
    * Computes the density function.
    * 
    */
   public static double density (double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      return x < 0 ? 0 : lambda*Math.exp (-lambda*x);
   }


   /**
    * Computes the  distribution function.
    * 
    */
   public static double cdf (double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (x <= 0.0)
         return 0.0;
      double y = lambda * x;
      if (y >= XBIG)
         return 1.0;
      return -Math.expm1 (-y);
   }


   /**
    * Computes the complementary distribution function.
    * 
    */
   public static double barF (double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (x <= 0.0)
         return 1.0;
      if (lambda*x >= XBIGM)
         return 0.0;
         return Math.exp (-lambda*x);
   }


   /**
    * Computes the inverse distribution function.
    * 
    */
   public static double inverseF (double lambda, double u) {
        if (lambda <= 0)
           throw new IllegalArgumentException ("lambda <= 0");
        if (u < 0.0 || u > 1.0)
            throw new IllegalArgumentException ("u not in [0,1]");
        if (u >= 1.0)
            return Double.POSITIVE_INFINITY;
        if (u <= 0.0)
            return 0.0;
        return -Math.log1p (-u)/lambda;
   }


   /**
    * Estimates the parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN> of the exponential distribution
    *   using the maximum likelihood method, from the <SPAN CLASS="MATH"><I>n</I></SPAN> observations
    *    <SPAN CLASS="MATH"><I>x</I>[<I>i</I>]</SPAN>, 
    * <SPAN CLASS="MATH"><I>i</I> = 0, 1,&#8230;, <I>n</I> - 1</SPAN>. The estimate is returned in a one-element
    *     array, as element 0.
    * 
    * @param x the list of observations used to evaluate parameters
    * 
    *    @param n the number of observations used to evaluate parameters
    * 
    *    @return returns the parameter [
    * <SPAN CLASS="MATH">hat(&lambda;)</SPAN>]
    * 
    */
   public static double[] getMLE (double[] x, int n)
   {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[];
      double sum = 0.0;
      parameters = new double[1];
      for (int i = 0; i < n; i++)
         sum+= x[i];
      parameters[0] = (double) n / sum;

      return parameters;
   }


   /**
    * Creates a new instance of an exponential distribution with parameter
    *    <SPAN CLASS="MATH"><I>&#955;</I></SPAN> estimated using
    *    the maximum likelihood method based on the <SPAN CLASS="MATH"><I>n</I></SPAN> observations <SPAN CLASS="MATH"><I>x</I>[<I>i</I>]</SPAN>,
    *    
    * <SPAN CLASS="MATH"><I>i</I> = 0, 1,&#8230;, <I>n</I> - 1</SPAN>.
    * 
    * @param x the list of observations to use to evaluate parameters
    * 
    *    @param n the number of observations to use to evaluate parameters
    * 
    * 
    */
   public static ExponentialDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new ExponentialDist (parameters[0]);
   }


   /**
    * Computes and returns the mean, 
    * <SPAN CLASS="MATH"><I>E</I>[<I>X</I>] = 1/<I>&#955;</I></SPAN>,
    *    of the exponential distribution with parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN>.
    * 
    * @return the mean of the exponential distribution 
    * <SPAN CLASS="MATH"><I>E</I>[<I>X</I>] = 1/<I>&#955;</I></SPAN>
    * 
    */
   public static double getMean (double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (1 / lambda);
   }


   /**
    * Computes and returns the variance, 
    * <SPAN CLASS="MATH">Var[<I>X</I>] = 1/<I>&#955;</I><SUP>2</SUP></SPAN>,
    *    of the exponential distribution with parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN>.
    * 
    * @return the variance of the Exponential distribution 
    * <SPAN CLASS="MATH">Var[<I>X</I>] = 1/<I>&#955;</I><SUP>2</SUP></SPAN>
    * 
    */
   public static double getVariance (double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (1 / (lambda * lambda));
   }


   /**
    * Computes and returns the standard deviation of the
    *    exponential distribution with parameter <SPAN CLASS="MATH"><I>&#955;</I></SPAN>.
    * 
    * @return the standard deviation of the exponential distribution
    * 
    */
   public static double getStandardDeviation (double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (1 / lambda);
   }


   /**
    * Returns the value of <SPAN CLASS="MATH"><I>&#955;</I></SPAN> for this object.
    * 
    */
   public double getLambda() {
      return lambda;
   }



   /**
    * Sets the value of <SPAN CLASS="MATH"><I>&#955;</I></SPAN> for this object.
    * 
    */
   public void setLambda (double lambda) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      this.lambda = lambda;
      supportA = 0.0;
   }


   /**
    * Return a table containing the parameters of the current distribution.
    * 
    * 
    */
   public double[] getParams () {
      double[] retour = {lambda};
      return retour;
   }


   public String toString () {
      return getClass().getSimpleName() + " : lambda = " + lambda;
   }

}
