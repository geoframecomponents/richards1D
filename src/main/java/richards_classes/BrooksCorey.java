/*
 * GNU GPL v3 License
 *
 * Copyright 2017  Niccolo` Tubini
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
 * Brooks and Corey' SWRC model
 * @author Niccolo' Tubini
 *
 */

public class BrooksCorey extends SoilParametrization {
	
	// Brooks - Corey's parameters
	private double[] n;
	private double[] psiD;
	
	
	/**
	 * General constructor to be used when there are several soil layers
	 * each one with its own parameters
	 */
	public BrooksCorey() {};
	
	
		
	/**
	 * Set parameters
	 */
	public void set(double[] n, double[] psiD, double[] par3, double[] par4, double[] par5, double[] psiStar1, double[] psiStar2, double[] psiStar3, double[] alphaSpecificStorage, double[] betaSpecificStorage, double[] thetaR, double[] thetaS, double[] kappaSaturation) {
		
		this.n = n;
		this.psiD = psiD; 
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.kappaSaturation = kappaSaturation;
		super.alphaSpecificStorage = alphaSpecificStorage;
		super.betaSpecificStorage = betaSpecificStorage;
		
		this.psiStar1 = this.psiD;
	}
	
	
	
	/**
	 * Compute the water for a given water suction value
	 * @param suction 
	 * @return theta water content at suction value
	 */
	public double waterContent(double suction, int i){
				
		if(suction <= this.psiD[i]) {
			this.theta = this.thetaR[i] + (this.thetaS[i]-this.thetaR[i]) * Math.pow((this.psiD[i]/suction), this.n[i]);
		} else if (this.psiD[i]<suction && suction<=0) {
			this.theta = this.thetaS[i];
		} else {
			this.theta = this.thetaS[i] + 9.81*( this.alphaSpecificStorage[i] + this.thetaS[i]*this.betaSpecificStorage[i])*suction;
		}
		return this.theta;
	}
	
	
	
	/**
	 * Compute the derivative of the water content for a given water suction value
	 * @param suction 
	 * @return dTheta the value of the moisture capacity 
	 */
	public double dWaterContent(double suction,int i){
		
		if (suction <= this.psiD[i]) {
		    this.dTheta = this.n[i]*(this.thetaS[i] - this.thetaR[i])/Math.abs(this.psiD[i]) * Math.pow(this.psiD[i]/suction,this.n[i]+1);
		} else if (this.psiD[i]<suction && suction<=0) {
		    this.dTheta = 0;
		} else {
			this.dTheta = + 9.81*( this.alphaSpecificStorage[i] + this.thetaS[i]*this.betaSpecificStorage[i]);
		}
		
		return this.dTheta;
	}
	
	
	
	/**
	 * Compute the unsaturated hydraulic conductivity for a given water suction value
	 * @param suction
	 * @return kappa hydraulic conductivity at suction value
	 */
	public double hydraulicConductivity(double suction,int i){
		
		this.saturationDegree = (waterContent(suction,i) - thetaR[i]) / (thetaS[i] - thetaR[i]);
		if(this.saturationDegree<1) {
			this.kappa = this.kappaSaturation[i] * Math.pow(this.saturationDegree, 3+2/this.n[i]);
		} else {
			this.kappa = this.kappaSaturation[i];
		}
		
		return this.kappa;
	}

	
	
	/**
	 * Evaluate the function theta_1 of the Jordan decomposition (Casulli and Zanolli, 2010) for a
	 * given water suction value
	 * @param suction
	 * @return theta1 
	 */
	public double pIntegral(double suction,int i){
		if(suction <= this.psiStar1[i]) {
			super.f1 = this.waterContent(suction,i);
		} else if (this.psiD[i]<suction && suction<=0) {
			this.f1 = this.waterContent(this.psiStar1[i],i) + this.dWaterContent(this.psiStar1[i],i)*(suction - this.psiStar1[i]);
		} else {
			this.f1 = this.waterContent(this.psiStar1[i],i) + this.dWaterContent(this.psiStar1[i],i)*(suction - this.psiStar1[i]) + this.dWaterContent(suction, i)*(suction-0);
		}

		return this.f1;
	}
	
	
	
	/**
	 * Evaluate the function theta_2 of the Jordan decomposition (Casulli and Zanolli, 2010) for a
	 * given water suction value
	 * @param suction
	 * @return theta2
	 */
	public double qIntegral(double suction,int i){
		super.f2 = pIntegral(suction,i) - this.waterContent(suction,i);

		return super.f2;
	}
	
	
	
	/**
	 * Evaluate the function p of the Jordan decomposition (Casulli and Zanolli, 2010) for a
	 * given water suction value
	 * @param suction
	 * @return dtheta1
	 */
	public double p(double suction,int i){
		if (suction <= this.psiStar1[i]) {
			super.df1 = this.dWaterContent(suction,i);
		} else if (this.psiD[i]<suction && suction<=0) {
			super.df1 = this.dWaterContent(this.psiStar1[i],i);
		} else {
			super.df1 = this.dWaterContent(this.psiStar1[i],i) + this.dWaterContent(suction, i);
		}

		return super.df1;
	}
	
	
	
	/**
	 * Evaluate the function q of the Jordan decomposition (Casulli and Zanolli, 2010) for a
	 * given water suction value
	 * @param suction
	 * @return dtheta2
	 */
	public double q(double suction,int i){
		super.df2 = p(suction,i) - this.dWaterContent(suction,i);

		return super.df2;
	}

}