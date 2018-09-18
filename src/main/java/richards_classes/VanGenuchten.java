/*
 * GNU GPL v3 License
 *
 * Copyright 2017 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package richards_classes;

/**
 * Van Genuchten's SWRC model
 * @author Niccolo' Tubini
 */

public class VanGenuchten extends SoilParametrization {
	
	// Van Genuchten's parameters
	private double n;
	private double m;
	private double alpha;
	
	
	/**
	 * General constructor to be used when there are several soil layers
	 * each one with its own parameters
	 */
	
	public VanGenuchten() {};
	
	
	/**
	 * 
	 * @param n		Van Genuchten's parameter      >1
	 * @param alpha Van Genuchten's parameter      >0
	 * @param thetaR residual water content	     ]0,1[
	 * @param thetaS water content at saturation ]0,1[
	 * @param kappaSaturation hydraulic conductivity at saturation >0
	 * @exception if n<1 error check the value
	 * @exception if alpha<0 error check the value
	 * @exception if thetaR>1 or thetaR<0 error check the value
	 * @exception if thetaS>1 or thetaR<0 error check the value
	 * @exception if thetaR > thetaS check the value: thetaR must be less then thetaS
	 * @exception if kappaSaturation <0 check the value
	 */
	
	public VanGenuchten(double n, double alpha, double thetaR, double thetaS, double kappaSaturation){
		
		this.n = n;
		this.alpha = alpha; 
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.kappaSaturation = kappaSaturation;
		
		if(this.n < 1){
			throw new IllegalArgumentException( "ERROR: Check the value of the Van Genuchten's parameter  n \n");
		}
		if(this.alpha < 0){
			throw new IllegalArgumentException( "ERROR: Check the value of the Van Genuchten's parameter  alpha \n");
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
		
		this.m = 1-1/this.n;
		this.psiStar1 = (-1.0/this.alpha)*Math.pow((this.n-1.0)/this.n,1.0/this.n);
	}
	
	
	
	public void set(double n, double alpha, double par3, double par4, double par5, double alphaSpecificStorage, double betaSpecificStorage, double thetaR, double thetaS, double kappaSaturation) {
		
		this.n = n;
		this.alpha = alpha; 
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.kappaSaturation = kappaSaturation;
		this.m = 1-1/this.n;
		super.alphaSpecificStorage = alphaSpecificStorage;
		super.betaSpecificStorage = betaSpecificStorage;
		
		this.psiStar1 = (-1.0/this.alpha)*Math.pow((this.n-1.0)/this.n,1.0/this.n);
	}
	

	
	/**
	 * @param suction 
	 * @return theta water content at suction value
	 */
	public double waterContent(double suction){
				
		if(suction <= 0) {
		    this.theta = this.thetaR + (this.thetaS - this.thetaR) / Math.pow(1.0 + Math.pow(Math.abs(this.alpha*suction), this.n), this.m);
		} else {
		    this.theta = this.thetaS + 9.81*( this.alphaSpecificStorage + this.thetaS*this.betaSpecificStorage)*suction;
		}

		return this.theta;
	}
	
	
	
	/**
	 * @param suction 
	 * @return dTheta the value of the moisture capacity 
	 */
	public double dWaterContent(double suction){
		
		if (suction <= 0) {
		    this.dTheta = this.alpha*this.n*this.m*(this.thetaS - this.thetaR) / Math.pow(1.0 + Math.pow(Math.abs(this.alpha*suction), this.n), this.m + 1.0)*Math.pow(Math.abs(this.alpha*suction), this.n - 1.0);
		} else {
		    this.dTheta =  9.81*( this.alphaSpecificStorage + this.thetaS*this.betaSpecificStorage);
		}
		
		return this.dTheta;
	}
	
	
		
	/**
	 * @param suction
	 * @return kappa hydraulic conductivity at suction value
	 */
	public double hydraulicConductivity(double suction){
		
		this.saturationDegree = (waterContent(suction) - thetaR) / (thetaS - thetaR); 
		if(this.saturationDegree<1) {
			this.kappa = this.kappaSaturation * Math.pow(this.saturationDegree, 0.5 ) * Math.pow(1.0 - Math.pow(1.0 - Math.pow(this.saturationDegree, 1.0/this.m), this.m), 2.0);
		} else {
			this.kappa = this.kappaSaturation;
		}
		return this.kappa;
	}

	
	
	/**
	 * @param suction
	 * @return theta1 
	 */
	public double pIntegral(double suction){
		if(suction <= this.psiStar1) {
			super.f1 = this.waterContent(suction);
		} else {
			this.f1 = this.waterContent(this.psiStar1) + this.dWaterContent(this.psiStar1)*(suction - this.psiStar1);
		}

		return this.f1;
	}
	
	
	
	/**
	 * @param suction
	 * @return theta2
	 */
	public double qIntegral(double suction){
		super.f2 = pIntegral(suction) - this.waterContent(suction);

		return super.f2;
	}
	
	/**
	 * @param suction
	 * @return dtheta1
	 */
	public double p(double suction){
		if (suction <= this.psiStar1) {
		    // left of critical value, take the original derivative
			super.df1 = this.dWaterContent(suction);
		}
		else {
		    // on the right of the critical value, keep the maximum derivative
			super.df1 = this.dWaterContent(this.psiStar1);
		}

		return super.df1;
	}
	
	
	
	/**
	 * @param suction
	 * @return dtheta2
	 */
	public double q(double suction){
		super.df2 = p(suction) - this.dWaterContent(suction);

		return super.df2;
	}



}
