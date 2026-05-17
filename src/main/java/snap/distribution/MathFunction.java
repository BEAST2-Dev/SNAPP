
package snap.distribution;

/**
 * This interface should be implemented by classes which represent univariate 
 * mathematical functions. It is used to pass an arbitrary function of one 
 * variable as argument to another function. For example, it is used
 * in {@link snap.distribution.RootFinder RootFinder} to find
 *  the zeros of a function.
 * 
 */
public interface MathFunction {


   /**
    * Returns the value of the function evaluated at <SPAN CLASS="MATH"><I>x</I></SPAN>.
    * 
    * @param x value at which the distribution function is evaluated
    * 
    *    @return function evaluated at <TT>x</TT>
    * 
    */
   public double evaluate (double x);

}
