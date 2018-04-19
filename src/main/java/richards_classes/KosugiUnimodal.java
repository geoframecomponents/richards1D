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
 * Kosugi's two-parameter lognormal distribution SWRC model
 * @author Niccolò Tubini
 */

import org.apache.commons.math3.special.Erf;

public class KosugiUnimodal extends SoilParametrization {
	
	private double rMedian; // radius median of pore-size distribution
	private double sigma;   // standard deviation of pore-size distribution
	private double psiMedian; // suction value related to rMedian by Young-Laplace equation

	
	
	public KosugiUnimodal(HydraulicConductivityTemperature kappaWithTemperature, double rMedian, double sigma, double thetaR, double thetaS, double kappaSaturation){
		this.rMedian = rMedian;
		this.sigma = sigma; 
		this.thetaR = thetaR;
		this.thetaS = thetaS;
		this.kappaSaturation = kappaSaturation;
		super.kappaTemperature = kappaWithTemperature;
		
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
		    this.theta = this.thetaR + (this.thetaS - this.thetaR)*0.5*( 1-Erf.erf(Math.log(suction/this.psiMedian)/(this.sigma*Math.pow(2,0.5))) ) ; // ho medificato la parte relativa al sigma dentro il logaritmo
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
		final double l = 0.5;   // (Kosugi, 1996)
		final double gamma = 2; // Mualem model (Kosugi, 1996)
 		final double eta = 1;   // Mualem model (Kosugi, 1996)
		this.saturationDegree = (waterContent(suction) - thetaR) / (thetaS - thetaR); 
		this.kappa = this.kappaSaturation * Math.pow(this.saturationDegree, l)*Math.pow( ( 0.5*Erf.erfc( Erf.erfcInv(2*this.saturationDegree) + eta*this.sigma/Math.sqrt(2)  ) ),gamma );
		
		return this.kappa;
	}
}
