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

import org.apache.commons.math3.special.Erf;

/**
 * Romano et al. bimodal retention curve
 * @author Niccolo' Tubini
 *
 */

public class Romano extends SoilParametrization {
	
	// Romano parameters
	private double w;
	private double psiM1;
	private double psiM2;
	private double sigma1;
	private double sigma2;
	
	private double gamma1;
	private double gamma2;
	private double aa;
	private double bb;
	private double r;
	/**
	 * General constructor to be used when there are several soil layers
	 * each one with its own parameters
	 */
	
	public Romano() {};
	
	
	
	/**
	 * 
	 * @param w		      >1
	 * @param psiM1  Brooks-Corey's parameter      >0
	 * @param psiM2 residual water content	     ]0,1[
	 * @param sigma1 water content at saturation ]0,1[
	 * @param sigma2 h
	 * @param kappaSaturation hydraulic conductivity at saturation >0
	 */
	public Romano(double w, double psiM1, double psiM2, double sigma1, double sigma2, double thetaR, double thetaS, double kappaSaturation){
		this.w = w;
		this.psiM1 = psiM1; 
		this.psiM2 = psiM2;
		this.sigma1 = sigma1;
		this.sigma2 = sigma2; 
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.kappaSaturation = kappaSaturation;
		/*
		if(this.n < 1){
			throw new IllegalArgumentException( "ERROR: Check the value of the Brooks-Corey's parameter  n \n");
		}
		if(this.psiD > 0){
			throw new IllegalArgumentException( "ERROR: Check the value of the Brooks-Corey's parameter  psiD \n");
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
		*/
		
	}
	
	
	
	public void set(double w, double psiM1, double psiM2, double sigma1, double sigma2, double alphaSpecificStorage, double betaSpecificStorage, double thetaR, double thetaS, double kappaSaturation) {
		
		this.w = w;
		this.psiM1 = psiM1; 
		this.psiM2 = psiM2;
		this.sigma1 = sigma1;
		this.sigma2 = sigma2;
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.kappaSaturation = kappaSaturation;
		super.alphaSpecificStorage = alphaSpecificStorage;
		super.betaSpecificStorage = betaSpecificStorage;
		
		//this.psiStar = this.psiD;
	}
	
	
	
	
	/**
	 * @param suction 
	 * @return theta water content at suction value
	 */
	public double waterContent(double suction){
				
		if(suction <= 0) {
			this.theta = this.thetaR + (this.thetaS-this.thetaR) * (this.w/2 * Erf.erfc(Math.log(suction/this.psiM1)/(this.sigma1*Math.sqrt(2)) ) + (1-w)/2* Erf.erfc(Math.log(suction/this.psiM2)/(this.sigma2*Math.sqrt(2)) ) );
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
			this.gamma1 = Math.exp( -( Math.log(suction/this.psiM1)/(this.sigma1*Math.sqrt(2))) );
			this.gamma2 = Math.exp( -( Math.log(suction/this.psiM2)/(this.sigma2*Math.sqrt(2))) );
		    this.dTheta = 1/(Math.sqrt(2*Math.PI)*suction/this.psiM1)*(this.thetaS - this.thetaR) * ( this.w/this.sigma1*this.gamma1 + (1-w)/this.sigma2*this.gamma2 );
		} else {
		    this.dTheta = + 9.81*( this.alphaSpecificStorage + this.thetaS*this.betaSpecificStorage);
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
			this.aa = ( Math.pow(this.sigma1, 2) + Math.log(suction/psiM1)) / ( this.sigma1*Math.sqrt(2) );
			this.aa = ( Math.pow(this.sigma2, 2) + Math.log(suction/psiM2)) / ( this.sigma2*Math.sqrt(2) );
			this.r = this.psiM1/this.psiM2 * (1-w)/w * Math.exp(0.5*(Math.pow(this.sigma1,2)-Math.pow(this.sigma2,2)))
			this.kappa = this.kappaSaturation*Math.sqrt(this.saturationDegree)*Math.pow( 0.5*Erf.erfc(this.aa)/(1+this.r) + 0.5*Erf.erfc(this.bb)/(1+this.r),2);
		} else {
			this.kappa = this.kappaSaturation;
		}
		
		return this.kappa;
	}
}