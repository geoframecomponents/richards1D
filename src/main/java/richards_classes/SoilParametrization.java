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
 * The soil prametrization abstract class.
 * @author Niccolò Tubini
 *
 */

public abstract class SoilParametrization {
	
	protected double psiStar; // suction value at which the derivative of hydraulic capacity is null
	protected double thetaR;  // residual water content 
	protected double thetaS;  // water content at saturation
	protected double kappaSaturation; // hydraulic conductivity at saturation
	
	protected double theta; // water content
	protected double dTheta;// first derivative of theta with respect of suction
	protected double saturationDegree; 
	protected double kappa; // hydraulic conductivity
	
	protected HydraulicConductivityTemperature kappaTemperature;
	
	
	/**
	 * This method return the value of suction at which the derivative
	 * of theta with respect of suction has its maximum.
	 * @return
	 */
	public double getPsiStar(){
		return this.psiStar;
	}
	
	
	/**
	 * This method compute the water content given the suction value 
	 * @param suction
	 * @return
	 */
	public abstract double waterContent(double suction);
	
	
	/**
	 * This method compute the value of the derivative of theta given the suction value
	 * @param suction
	 * @return
	 */
	public abstract double dWaterContent(double suction);
	
	
	/**
	 * This method compute the hydraulic conductivity accordingly with Mualem's assumption.
	 * @param suction
	 * @return
	 */
	public abstract double hydraulicConductivity(double suction);
	
	
	
	/**
	 * This method compute the hydraulic conductivity accordingly with Mualem's assumption and considering the effect of temperature.
	 * @param suction
	 * @param temperature
	 * @return
	 */
	public double hydraulicConductivity(double suction, double temperature) {
		
		kappa = hydraulicConductivity(suction)*kappaTemperature.temperatureCorrection(temperature);
		
		return kappa;
				
	}
	
	
	/**
	 * This method creates a data set to plot the hydraulic properties of the soil
	 * SWRC(psi), hydraulic conductivity(Se), moisture capacity(psi)
	 */
	public double[][] hydraulicModelCurves(){
		
		double[][] result  = new double[200][6];

		for(int i=0; i<result.length; i++){
			result[i][0] = (double)(-i);
			result[i][1] = (double)(i)/200;
		}
				
		for(int i=0; i<result.length; i++){
			result[i][2] = waterContent(result[i][0]);
			result[i][3] = dWaterContent(result[i][0]);
			result[i][4] = hydraulicConductivity(result[i][0]);
			
		}
		
	return result;
	} 
}