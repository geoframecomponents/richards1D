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
public class JordanDecompositionTotalDepth extends JordanDecomposition {

	protected TotalDepth totalDepth;
	
	public JordanDecompositionTotalDepth(Object myFunction){
		totalDepth = (TotalDepth) myFunction;
	}
	
	
	@Override
	public void setSoilParametrization(double par1SWRC, double par2SWRC, double thetaR, double thetaS) {
		// TODO Auto-generated method stub
		
	}
	
	/**
	 * @param suction
	 * @return H1 
	 */
	public double pIntegral(double suction){
		this.f1 = totalDepth.totalDepth(suction);

		return this.f1;
	}
	
	/**
	 * @param suction
	 * @return H2
	 */
	public double qIntegral(double suction){
		this.f2 = 0;

		return this.f2;
	}
	
	/**
	 * @param suction
	 * @return dH1
	 */
	public double p(double suction){
		this.df1 = totalDepth.dTotalDepth(suction);

		return this.df1;
	}
	
	/**
	 * @param suction
	 * @return dH2
	 */
	public double q(double suction){
		this.df2 = 0;

		return this.df2;
	}



}