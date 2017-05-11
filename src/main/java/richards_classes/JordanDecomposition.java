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
public class JordanDecomposition {
	protected double theta1;
	protected double theta2;
	protected double dTheta1;
	protected double dTheta2;
	protected SoilParametrization soilPar;
	
	public JordanDecomposition(SoilParametrization sp){
		soilPar = sp;
	}
	
	/**
	 * @param suction
	 * @return theta1 
	 */
	public double waterContent1(double suction){
		if(suction <= soilPar.psiStar) {
			this.theta1 = soilPar.waterContent(suction);
		} else {
			this.theta1 = soilPar.waterContent(soilPar.psiStar) + soilPar.dWaterContent(soilPar.psiStar)*(suction - soilPar.psiStar);
		}

		return this.theta1;
	}
	
	/**
	 * @param suction
	 * @return theta2
	 */
	public double waterContent2(double suction){
		this.theta2 = waterContent1(suction) - soilPar.waterContent(suction);

		return this.theta2;
	}
	
	/**
	 * @param suction
	 * @return dtheta1
	 */
	public double dWaterContent1(double suction){
		if (suction <= soilPar.psiStar) {
		    // left of critical value, take the original derivative
		    this.dTheta1 = soilPar.dWaterContent(suction);
		}
		else {
		    // on the right of the critical value, keep the maximum derivative
		    this.dTheta1 = soilPar.dWaterContent(soilPar.psiStar);
		}

		return this.dTheta1;
	}
	
	/**
	 * @param suction
	 * @return dtheta2
	 */
	public double dWaterContent2(double suction){
		this.dTheta2 = dWaterContent1(suction) - soilPar.dWaterContent(suction);

		return this.dTheta2;
	}

}