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

import java.util.ArrayList;
import java.util.LinkedHashMap;

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
	 * This method creates a data set to plot the hydraulic properties of the soil
	 * SWRC(psi), hydraulic conductivity(Se), moisture capacity(psi)
	 */
	public double[][] hydraulicModelCurves(){
		
		double[][] result  = new double[200][5];

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
	
	/**
	 * This method creates a data set to plot the hydraulic properties of the soil
	 * SWRC(psi), hydraulic conductivity(Se), moisture capacity(psi)
	 */
	public LinkedHashMap<String, double[]> hydraulicModelCurves1(){
		
		LinkedHashMap<String, double[]> result  = new LinkedHashMap<String, double[]>();
		double[] tempValues = new double[200];
		//System.out.println("\n\nPsi");
		for(int i=0; i<tempValues.length; i++){
			tempValues[i] = (double)(-i);
			//System.out.println(tempValues[i]);
		}
		result.put("Psi [m]", tempValues.clone());
		
		//System.out.println("\n\nTheta");
		for(int i=0; i<result.get("Psi [m]").length; i++){
			tempValues[i] = waterContent(result.get("Psi [m]")[i]);
			//System.out.println(tempValues[i]);

		}
		result.put("Theta [-]", tempValues.clone());
		//System.out.println("\n\ndTheta");
		for(int i=0; i<result.get("Psi [m]").length; i++){
			//System.out.println(result.get("Psi [m]")[i]);
			tempValues[i] = dWaterContent(result.get("Psi [m]")[i]);
			//System.out.println("	"+tempValues[i]);
		}
		result.put("dTheta [1/m]", tempValues.clone());
		
		//System.out.println("\n\nK");
		for(int i=0; i<result.get("Psi [m]").length; i++){
			tempValues[i] = hydraulicConductivity(result.get("Psi [m]")[i]);
			//System.out.println(tempValues[i]);
		}
		result.put("K [m/s]", tempValues.clone());
		
	return result;
	} 
}