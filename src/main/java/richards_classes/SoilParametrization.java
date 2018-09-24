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
 * @author Niccolo' Tubini
 *
 */

public abstract class SoilParametrization {
	
	protected double[] psiStar1; // suction value at which the derivative of hydraulic capacity is null
	protected double[] psiStar2; // suction value at which the derivative of hydraulic capacity is null
	protected double[] psiStar3; // suction value at which the derivative of hydraulic capacity is null
	protected double[] thetaR;  // residual water content 
	protected double[] thetaS;  // water content at saturation
	protected double[] kappaSaturation; // hydraulic conductivity at saturation
	protected double[] alphaSpecificStorage;
	protected double[] betaSpecificStorage;
	
	protected double theta; // water content
	protected double dTheta;// first derivative of theta with respect of suction
	protected double saturationDegree; 
	protected double kappa; // hydraulic conductivity
	protected double f1;
	protected double f2;
	protected double df1;
	protected double df2;
	
	
	/**
	 * This method set SWRC parameters
	 * @param par1
	 * @param par2
	 * @param par3
	 * @param par4
	 * @param par5
	 * @param thetaR adimensional residual water content
	 * @param thetaS adimensional water content at saturation
	 * @param kappaSaturation hydraulic conductivity at saturation
	 * 
	 */
	//public abstract void set(double par1, double par2, double par3, double par4, double par5, double alphaSpecificStorage, double betaSpecificStorage, double thetaR, double thetaS, double kappaSaturation);
	
	
	
	/**
	 * This method set SWRC parameters
	 * @param par1
	 * @param par2
	 * @param par3
	 * @param par4
	 * @param par5
	 * @param thetaR adimensional residual water content
	 * @param thetaS adimensional water content at saturation
	 * @param kappaSaturation hydraulic conductivity at saturation
	 * 
	 */
	public abstract void set(double[] par1, double[] par2, double[] par3, double[] par4, double[] par5, double[] psiStar1, double[] psiStar2, double[] psiStar3, double[] alphaSpecificStorage, double[] betaSpecificStorage, double[] thetaR, double[] thetaS, double[] kappaSaturation);
	
	
	
	/**
	 * This method return the value of suction at which the derivative
	 * of moisture capacity with respect of suction has its maximum.
	 * @return
	 */
	public double getPsiStar1(int i){
		return this.psiStar1[i];
	}
	
	
	/**
	 * This method return the value of suction at which the derivative
	 * of moisture capacity with respect of suction has its maximum.
	 * @return
	 */
	//public double getPsiStar2(){
	//	return this.psiStar2;
	//}
	
	
	/**
	 * This method compute the water content given the suction value 
	 * @param suction
	 * @return
	 */
	public abstract double waterContent(double suction,int i);
	
	
	/**
	 * This method compute the value of the derivative of theta given the suction value
	 * @param suction
	 * @return
	 */
	public abstract double dWaterContent(double suction,int i);
	
	
	/**
	 * This method compute the hydraulic conductivity accordingly with Mualem's assumption.
	 * @param suction
	 * @return
	 */
	public abstract double hydraulicConductivity(double suction,int i);
	
	
	/**
	 * This method compute the integral of p function.
	 * @param suction
	 * @return
	 */
	public abstract double pIntegral(double suction, int i);
	
	
	/**
	 * This method compute the integral of q function.
	 * @param suction
	 * @return
	 */
	public abstract double qIntegral(double suction, int i);
	
	
	/**
	 * This method compute the p function.
	 * @param suction
	 * @return
	 */
	public abstract double p(double suction, int i);
	
	
	/**
	 * This method compute the q function.
	 * @param suction
	 * @return
	 */
	public abstract double q(double suction,int i);
	
	
	// questo metodo e` da mettere nello script python della creazione del dominio di calcolo
	/**
	 * This method creates a data set to plot the hydraulic properties of the soil
	 * SWRC(psi), hydraulic conductivity(Se), moisture capacity(psi)
	 */
	// DEPRECATED
	/*
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
	*/
	
	
	// To be reviewed since this method has to be applied to each layer
	/**
	 * This method creates a data set to plot the hydraulic properties of the soil
	 * SWRC(psi), hydraulic conductivity(Se), moisture capacity(psi)
	 */
	/*
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
	*/
}