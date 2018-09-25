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
	 * 
	 * @param n		Brooks-Corey's parameter      >1
	 * @param psiD  Brooks-Corey's parameter      >0
	 * @param thetaR residual water content	     ]0,1[
	 * @param thetaS water content at saturation ]0,1[
	 * @param kappaSaturation hydraulic conductivity at saturation >0
	 * @exception if n<1 error check the value
	 * @exception if psiD>0 error check the value
	 * @exception if thetaR>1 or thetaR<0 error check the value
	 * @exception if thetaS>1 or thetaR<0 error check the value
	 * @exception if thetaR > thetaS check the value: thetaR must be less then thetaS
	 * @exception if kappaSaturation <0 check the value
	 */
	/*
	public BrooksCorey(double n, double psiD, double thetaR, double thetaS, double kappaSaturation){
		this.n = n;
		this.psiD = psiD; 
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.kappaSaturation = kappaSaturation;
		
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
	
	
	
	public void set(double[] n, double[] psiD, double[] par3, double[] par4, double[] par5, double[] psiStar1, double[] psiStar2, double[] psiStar3, double[] alphaSpecificStorage, double[] betaSpecificStorage, double[] thetaR, double[] thetaS, double[] kappaSaturation) {
		
		this.n = n;
		this.psiD = psiD; 
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.kappaSaturation = kappaSaturation;
		super.alphaSpecificStorage = alphaSpecificStorage;
		super.betaSpecificStorage = betaSpecificStorage;
		
		//super.psiStar1 = this.psiD;
		this.psiStar1 = this.psiD;
	}
	
	
	
	
	/**
	 * @param suction 
	 * @return theta water content at suction value
	 */
	public double waterContent(double suction, int i){
				
		if(suction <= this.psiD[i]) {
			this.theta = this.thetaR[i] + (this.thetaS[i]-this.thetaR[i]) * Math.pow((this.psiD[i]/suction), this.n[i]);
		}//else if (this.psiD[i]<suction && suction<0){
		//	this.theta = this.thetaS[i];
		else {
			this.theta = this.thetaS[i];// + 9.81*( this.alphaSpecificStorage[i] + this.thetaS[i]*this.betaSpecificStorage[i])*suction;
		}
		return this.theta;
	}
	
	
	
	/**
	 * @param suction
	 * @return dTheta the value of the moisture capacity
	 */
	public double dWaterContent(double suction,int i){
		
		if (suction <= this.psiD[i]) {
		    this.dTheta = this.n[i]*(this.thetaS[i] - this.thetaR[i])/Math.abs(this.psiD[i]) * Math.pow(this.psiD[i]/suction,this.n[i]+1);
		} else {
		    this.dTheta = 0;
		}// else {
		//	this.dTheta = + 9.81*( this.alphaSpecificStorage[i] + this.thetaS[i]*this.betaSpecificStorage[i]);
		//}
		
		return this.dTheta;
	}
	
	
	
	/**
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
	 * @param suction
	 * @return theta1 
	 */
	public double pIntegral(double suction,int i){
		if(suction <= this.psiStar1[i]) {
			super.f1 = this.waterContent(suction,i);
		} else {
			this.f1 = this.waterContent(this.psiStar1[i],i) + this.dWaterContent(this.psiStar1[i],i)*(suction - this.psiStar1[i]);
		}

		return this.f1;
	}
	
	
	
	/**
	 * @param suction
	 * @return theta2
	 */
	public double qIntegral(double suction,int i){
		super.f2 = pIntegral(suction,i) - this.waterContent(suction,i);

		return super.f2;
	}
	
	/**
	 * @param suction
	 * @return dtheta1
	 */
	public double p(double suction,int i){
		if (suction <= this.psiStar1[i]) {
		    // left of critical value, take the original derivative
			super.df1 = this.dWaterContent(suction,i);
		}
		else {
		    // on the right of the critical value, keep the maximum derivative
			super.df1 = this.dWaterContent(this.psiStar1[i],i);
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