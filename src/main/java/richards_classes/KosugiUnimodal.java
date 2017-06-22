package richards_classes;

import org.apache.commons.math3.special.Erf;

public class KosugiUnimodal extends SoilParametrization {
	private double rMedian; // radius median of pore-size distribution
	private double sigma;   // standard deviation of pore-size distribution
	private double psiMedian; // suction value related to rMedian by Young-Laplace equation
		
	public KosugiUnimodal(double rMedian, double sigma, double thetaR, double thetaS, double kappaSaturation){
		this.rMedian = rMedian;
		this.sigma = sigma; 
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.kappaSaturation = kappaSaturation;
		
		if(this.rMedian < 0){
			throw new IllegalArgumentException( "ERROR: Check the value of the median radius of the pore-size distribution rMedian \n");
		}
		if(this.sigma < 0){
			throw new IllegalArgumentException( "ERROR: Check the value of the variance of the pore-size distribution sigma \n");
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
		
		this.psiStar = -1.49*Math.pow(10, -5)/this.rMedian/Math.exp(Math.pow(this.sigma,2));  // see Brutsaert, 1996
		this.psiMedian = -1.49*Math.pow(10, -5)/this.rMedian;
	}
	
	/**
	 * @param suction 
	 * @return theta water content at suction value
	 */
	public double waterContent(double suction){
				
		if(suction <= 0) {
		    this.theta = this.thetaR + (this.thetaS - this.thetaR)*0.5*( 1-Erf.erf(Math.log(suction/this.psiMedian)*2*this.sigma) ) ;
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
		
		if (suction < 0) {
		    this.dTheta = (this.thetaS-this.thetaR)/(Math.sqrt(2*Math.PI)*this.sigma*(-suction)) * Math.exp(-Math.pow(Math.log(suction/this.psiMedian),2)/(2*Math.pow(this.sigma,2)));
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
		final double l = 1;     // (Kosugi, 2002)
		final double gamma = 2; // Mualem model (Kosugi, 2002)
 		final double eta = 1;   // Mualem model (Kosugi, 2002)
		this.saturationDegree = (waterContent(suction) - thetaR) / (thetaS - thetaR); 
		this.kappa = Math.pow(this.saturationDegree, l)*Math.pow( ( 0.5*(1-Erf.erf( Math.pow((1-Erf.erf(this.saturationDegree)),-1) + eta*this.sigma ) ) ),gamma );
		
		return this.kappa;
	}
}
