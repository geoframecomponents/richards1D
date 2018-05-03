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
public class JordanDecompositionSoilMoisture extends JordanDecomposition {

	protected SoilParametrization soilPar;
	
	public JordanDecompositionSoilMoisture(Object myFunction){
		soilPar = (SoilParametrization) myFunction;
	}
	
	public void setSoilParametrization (double par1SWRC, double par2SWRC, double thetaR, double thetaS) {
		this.soilPar.set(par1SWRC, par2SWRC, thetaR, thetaS, -999);
	}
	
	/**
	 * @param suction
	 * @return theta1 
	 */
	public double pIntegral(double suction){
		if(suction <= this.soilPar.getPsiStar()) {
			this.f1 = this.soilPar.waterContent(suction);
		} else {
			this.f1 = this.soilPar.waterContent(this.soilPar.getPsiStar()) + this.soilPar.dWaterContent(this.soilPar.getPsiStar())*(suction - this.soilPar.getPsiStar());
		}

		return this.f1;
	}
	
	/**
	 * @param suction
	 * @return theta2
	 */
	public double qIntegral(double suction){
		this.f2 = pIntegral(suction) - this.soilPar.waterContent(suction);

		return this.f2;
	}
	
	/**
	 * @param suction
	 * @return dtheta1
	 */
	public double p(double suction){
		if (suction <= this.soilPar.getPsiStar()) {
		    // left of critical value, take the original derivative
		    this.df1 = this.soilPar.dWaterContent(suction);
		}
		else {
		    // on the right of the critical value, keep the maximum derivative
		    this.df1 = this.soilPar.dWaterContent(this.soilPar.getPsiStar());
		}

		return this.df1;
	}
	
	/**
	 * @param suction
	 * @return dtheta2
	 */
	public double q(double suction){
		this.df2 = p(suction) - this.soilPar.dWaterContent(suction);

		return this.df2;
	}

}