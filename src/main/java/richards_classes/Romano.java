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
	private double[] w;
	private double[] psiM1;
	private double[] psiM2;
	private double[] sigma1;
	private double[] sigma2;
	
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
	/*
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
		
		
	}
	*/
	
	
	
	public void set(double[] w, double[] sigma1, double[] sigma2, double[] psiM1, double[] psiM2, double[] psiStar1, double[] psiStar2, double[] psiStar3, double[] alphaSpecificStorage, double[] betaSpecificStorage, double[] thetaR, double[] thetaS, double[] kappaSaturation) {
		
		this.w = w;
		this.psiM1 = psiM1; 
		this.psiM2 = psiM2;
		this.sigma1 = sigma1;
		this.sigma2 = sigma2;
		super.psiStar1 = psiStar1;
		super.psiStar2 = psiStar2;
		super.psiStar3 = psiStar3;
		super.thetaR = thetaR;
		super.thetaS = thetaS;
		super.kappaSaturation = kappaSaturation;
		super.alphaSpecificStorage = alphaSpecificStorage;
		super.betaSpecificStorage = betaSpecificStorage;
		
	}
	
	
	
	
	/**
	 * @param suction 
	 * @return theta water content at suction value
	 */
	public double waterContent(double suction,int i){
				
		if(suction <= 0) {
			this.theta = this.thetaR[i] + (this.thetaS[i]-this.thetaR[i]) * ( this.w[i]/2 * Erf.erfc( Math.log(suction/this.psiM1[i])/(this.sigma1[i]*Math.sqrt(2)) ) + (1-w[i])/2* Erf.erfc(Math.log(suction/this.psiM2[i])/(this.sigma2[i]*Math.sqrt(2)) ) );
		} else {
			this.theta = this.thetaS[i] + 9.81*( this.alphaSpecificStorage[i] + this.thetaS[i]*this.betaSpecificStorage[i])*suction;
		}
		return this.theta;
	}
	
	
	
	/**
	 * @param suction
	 * @return dTheta the value of the moisture capacity
	 */
	public double dWaterContent(double suction,int i){
		
		if (suction <= 0) {
			this.gamma1 = Math.exp( -Math.pow(( Math.log(suction/this.psiM1[i])/(this.sigma1[i]*Math.sqrt(2))),2) );
			this.gamma2 = Math.exp( -Math.pow(( Math.log(suction/this.psiM2[i])/(this.sigma2[i]*Math.sqrt(2))),2) );
		    this.dTheta = 1/(Math.sqrt(2*Math.PI)*suction/this.psiM1[i])*(this.thetaS[i] - this.thetaR[i]) * ( this.w[i]/this.sigma1[i]*this.gamma1 + (1-this.w[i])/this.sigma2[i]*this.gamma2 );
		} else {
		    this.dTheta = + 9.81*( this.alphaSpecificStorage[i] + this.thetaS[i]*this.betaSpecificStorage[i]);
		}
		
		return this.dTheta;
	}
	
	
	
	/**
	 * @param suction 
	 * @return kappa hydraulic conductivity at suction value
	 */
	public double hydraulicConductivity(double suction,int i){
		
		this.saturationDegree = (waterContent(suction,i) - super.thetaR[i]) / (super.thetaS[i] - super.thetaR[i]);
		if(this.saturationDegree<1) {
			this.aa = ( Math.pow(this.sigma1[i], 2) + Math.log(suction/this.psiM1[i]) ) / ( this.sigma1[i]*Math.sqrt(2) );
			this.bb = ( Math.pow(this.sigma2[i], 2) + Math.log(suction/this.psiM2[i]) ) / ( this.sigma2[i]*Math.sqrt(2) );
			this.r = this.psiM1[i]/this.psiM2[i] * (1-this.w[i])/this.w[i] * Math.exp(0.5*(Math.pow(this.sigma1[i],2)-Math.pow(this.sigma2[i],2)));
			this.kappa = this.kappaSaturation[i]*Math.sqrt(this.saturationDegree)*Math.pow( 0.5*Erf.erfc(this.aa)/(1+this.r) + 0.5*Erf.erfc(this.bb)/(1+1/this.r),2);
		} else {
			this.kappa = this.kappaSaturation[i];
		}
		
		return this.kappa;
	}
	
	
	
	/**
	 * @param suction
	 * @return theta1 
	 */
	public double pIntegral(double suction,int i){
		if(suction <= this.psiStar1[i]) {
			super.f1 = this.waterContent(suction,i);
		} else if(super.psiStar1[i]<suction && suction<=super.psiStar3[i]) {
			this.f1 = this.waterContent(super.psiStar1[i],i) + this.dWaterContent(this.psiStar1[i],i)*(suction - this.psiStar1[i]);
		} else if(super.psiStar3[i]<suction && suction<=super.psiStar2[i]) {
			//qui
			this.f1 = this.waterContent(suction, i) - this.waterContent(super.psiStar3[i], i) + this.waterContent(super.psiStar1[i], i) + this.dWaterContent(super.psiStar1[i], i)*(suction-super.psiStar1[i]);
		} else if(super.psiStar2[i]<suction && suction<=0) {
			this.f1 = this.waterContent(super.psiStar1[i], i) + this.waterContent(super.psiStar2[i], i) - this.waterContent(super.psiStar3[i], i) + this.dWaterContent(super.psiStar2[i], i)*(suction-super.psiStar2[i]) + this.dWaterContent(super.psiStar1[i], i)*(suction-super.psiStar1[i]);
		} else {
			this.f1 = this.waterContent(super.psiStar1[i], i) + this.waterContent(super.psiStar2[i], i) - this.waterContent(super.psiStar3[i], i) + this.dWaterContent(super.psiStar2[i], i)*(suction-super.psiStar2[i]) + this.dWaterContent(super.psiStar1[i], i)*(suction-super.psiStar1[i]) + this.dWaterContent(suction, i)*(suction-0);
		}
		return this.f1;
	}
	
	
	
	/**
	 * @param suction
	 * @return theta2
	 */
	public double qIntegral(double suction, int i){
		super.f2 = pIntegral(suction,i) - this.waterContent(suction,i);

		return super.f2;
	}
	
	/**
	 * @param suction
	 * @return dtheta1
	 */
	public double p(double suction,int i){
		if (suction <= this.psiStar1[i]) {
			super.df1 = this.dWaterContent(suction,i);
		} else if(super.psiStar1[i]<suction && suction<=super.psiStar3[i]) {
			super.df1 = this.dWaterContent(this.psiStar1[i],i);
		} else if (super.psiStar3[i]<suction && suction <=super.psiStar2[i]) {
			super.df1 = this.dWaterContent(suction, i) + this.dWaterContent(super.psiStar1[i], i);
		} else if (super.psiStar2[i]<suction && suction <=0){
			super.df1 = this.dWaterContent(super.psiStar2[i], i) + this.dWaterContent(super.psiStar1[i], i);
		} else {
			super.df1 = this.dWaterContent(super.psiStar2[i], i) + this.dWaterContent(super.psiStar1[i], i) + this.dWaterContent(suction, i);
		}
			
		return super.df1;
	}
	
	
	
	/**
	 * @param suction
	 * @return dtheta2
	 */
	public double q(double suction,int i){
		super.df2 = p(suction,i) - this.dWaterContent(suction,i);

		return super.df2;
	}

}
