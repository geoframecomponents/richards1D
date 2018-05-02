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
 * @author Niccolò Tubini
 *
 */

public class BrooksCorey extends SoilParametrization {
	
	// Brooks - Corey's parameters
	private double n;
	private double psiD;
	
	
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
	
	
	
	public void set(double n, double psiD, double thetaR, double thetaS, double kappaSaturation) {
		
		this.n = n;
		this.psiD = psiD; 
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.kappaSaturation = kappaSaturation;
		// DA CALCOLARE!!!
		//this.psiStar = (-1.0/this.alpha)*Math.pow((this.n-1.0)/this.n,1.0/this.n);
	}
	
	
	
	
	/**
	 * @param suction 
	 * @return theta water content at suction value
	 */
	public double waterContent(double suction){
				
		if(suction <= this.psiD) {
			this.theta = this.thetaR + (this.thetaS-this.thetaR) * Math.pow((this.psiD/suction), this.n);
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
		
		if (suction <= this.psiD) {
		    this.dTheta = this.n*(this.thetaS - this.thetaR) /Math.abs(this.psiD) * Math.pow(this.psiD/suction,this.n+1);
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
		
		this.saturationDegree = (waterContent(suction) - thetaR) / (thetaS - thetaR);
		this.kappa = this.kappaSaturation * Math.pow(this.saturationDegree, 3+2/this.n);
		return this.kappa;
	}
}