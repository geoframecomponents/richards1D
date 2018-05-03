package richards_classes;
/**
 * Since soil moisture is a nonnegative function with bounded variations, it is almost
 * everywhere differentiable, admit only discontinuities of the first kind, and can be
 * expressed as the difference of two nonnegative, nondecreasing, and bounded functions
 * (the Jordan decomposition [8]), say p(suction) and q(suction), so that c(suction) = p(suction)-q(suction) >= 0
 * and 0 <= q(suction) <= p(suction) for all suction values.
 * Look at: A NESTED NEWTON-TYPE ALGORITHM FOR FINITE VOLUME METHODS SOLVING RICHARDS’ EQUATION IN MIXED FORM, Casulli V., 2010
 *  
 * Here:
 * p(suction) is called dTheta1 and is computed by dWaterContent1
 * q(suction) is called dTheta2 and is computed by dWaterContent2
 * 
 */
public abstract class JordanDecomposition {
	protected double f1;
	protected double f2;
	protected double df1;
	protected double df2;
	protected Object myFunction;
	
	//public JordanDecomposition(Object myFunction){
	//	this.myFunction = myFunction;
	//}
	
	public abstract void setSoilParametrization (double par1SWRC, double par2SWRC, double thetaR, double thetaS);
	
	/**
	 * @param variable
	 * @return integral of p 
	 */
	public abstract double pIntegral(double variable);
	
	
	
	/**
	 * @param variable
	 * @return integral of q
	 */
	public abstract double qIntegral(double suction);
	
	/**
	 * @param variable
	 * @return p
	 */
	public abstract double p(double suction);
	
	/**
	 * @param variable
	 * @return q
	 */
	public abstract double q(double suction);

}