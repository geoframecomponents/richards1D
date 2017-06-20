package richards_classes;

public class Gardner extends SoilParametrization {
	
	// Gardner's parameter
	private double alpha;
	
	/**
	 * 
	 * @param alpha Gardner's parameter      >0
	 * @param thetaR residual water content	     ]0,1[
	 * @param thetaS water content at saturation ]0,1[
	 * @param kappaSaturation hydraulic conductivity at saturation >0
	 * @exception if alpha<0 error check the value
	 * @exception if thetaR>1 or thetaR<0 error check the value
	 * @exception if thetaS>1 or thetaR<0 error check the value
	 * @exception if thetaR > thetaS check the value: thetaR must be less then thetaS
	 * @exception if kappaSaturation <0 check the value
	 */
	
	public Gardner(double alpha, double thetaR, double thetaS, double kappaSaturation){
		
		this.alpha = alpha; 
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.kappaSaturation = kappaSaturation;

		if(this.alpha < 0){
			throw new IllegalArgumentException( "ERROR: Check the value of the Gardner's parameter  alpha \n");
		}
		if(this.thetaR < 0 | this.thetaR > 1 ){
			throw new IllegalArgumentException( "ERROR: Check the value of residual water content \n");
		}
		if(this.thetaS < 0 | this.thetaS > 1 ){
			throw new IllegalArgumentException( "ERROR: Check the value of water content at saturation \n");
		}
		if(this.thetaR >this.thetaS){
			throw new IllegalArgumentException( "ERROR: Check the value of residual water content or the value of water content at saturation \n");
		}
		if(this.kappaSaturation <0){
			throw new IllegalArgumentException( "ERROR: Check the value of hydraulic conductivity at saturation \n");
		}
		
		this.psiStar = 0;
	}
	
	/**
	 * @param suction 
	 * @return theta water content at suction value
	 */
	public double waterContent(double suction){
				
		if(suction <= 0) {
		    this.theta = this.thetaR + (this.thetaS - this.thetaR)*Math.exp(this.alpha*suction);
		} else {
		    this.theta = this.thetaS;
		}

		return this.theta;
	}
	
	/**
	 * @param suction
	 * @return dTheta the value of the moisture capacity 
	 */
	public double dWaterContent(double suction){
		
		if (suction <= 0) {
		    this.dTheta = this.alpha*(this.thetaS - this.thetaR)*Math.exp(this.alpha*suction);
		} else {
		    this.dTheta = 0;
		}
		
		return this.dTheta;
	}
	
	/**
	 * @param suction
	 * @return kappa hydraulic conductivity at suction value
	 */
	public double hydraulicConductivity(double suction){
		
		if (suction <= 0) {
			this.kappa = this.kappaSaturation * Math.exp(this.alpha*suction);
		} else {
			this.kappa = this.kappaSaturation;
		}
		
		return this.kappa;
	}
}
