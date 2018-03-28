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
package thermalParameterization;

import richards_classes.SimpleSoilParametrizationFactory;
import richards_classes.SoilParametrization;

/**
 * Integration of the energy conservation equation within the soil requires the definition of 
 * parametrizations to evaluate the heat capacity coefficient and the heat diffusion coefficient. Indeed, it depends on soil 
 * mineral composition, organic matter fraction, and of course the water content.
 * Since in literature different parametrization are proposed it is usefull defining an abstract class
 * to switch from a parametrization to another one with much less effort.
 * 
 * Note: extending to freezing-thawing processes one has to take into account of the ice content but also 
 * of the so called apparent heat capacity due to the water phase change.
 * 
 * @author Niccolò Tubini
 *
 */

public abstract class ThermalSoilParameterization {
	
	final double waterSpecificHeatCapacity = 4188; //[J/(Kg K)]
	final double waterDensity = 1000; //[Kg/m3]
	
	protected double sandFraction;
	protected double clayFraction;
	protected double thetaS;
	protected double cT;
	protected SoilParametrization soilPar;
	
	
	public ThermalSoilParameterization(double thetaS, double sandFraction, double clayFraction, SoilParametrization soilPar){
		this.thetaS = thetaS;
		this.sandFraction = sandFraction;
		this.clayFraction = clayFraction;
		this.soilPar = soilPar;
	}
	
	public abstract double thermalConductivity(double suction);
	
	/**
	 * This method computes the heat capacity as the sum of the soil particles heat capacity per unit volume
	 * and the that due to the water content.
	 * It is worthwhile that the SWRC is the same used to solve Richards' equation: this is garanteed by the
	 * constructor.
	 * @param suction
	 * @return heatCapacity
	 */
	public double heatCapacity(double suction){
		
		this.cT = waterDensity*waterSpecificHeatCapacity*soilPar.waterContent(suction) + (1-thetaS)*(2.128*sandFraction+2.385*clayFraction)/(sandFraction+clayFraction)*Math.pow(10, 6);
		
		return cT;
		
	}

	
	/**
	 * This method creates a data set to plot the thermal properties of the soil
	 * C_T(psi)
	 */
	public double[][] thermalSoilPropertiesCurves(){
		
		double[][] result  = new double[200][5];

		for(int i=0; i<result.length; i++){
			result[i][0] = (double)(-i);
		}
				
		for(int i=0; i<result.length; i++){
			result[i][1] = heatCapacity(result[i][0]);
			result[i][2] = thermalConductivity(result[i][0]);
			result[i][3] = -999;
			result[i][4] = -999;
		}
		
	return result;
	} 
	
	
	/**
	 * Just to test this class
	 * @param args
	 */
	public static void main(String[] args) {
		double thetaS = 0.41;
		double thetaR = 0.095;
		double alpha = 1.9;
		double n = 1.1;
		double psiE = -999;
		double lambda = -999;
		double rMedian = -999;
		double sigma = -999;
		double ks = 1;
		double[][] thermalParametrization;
		String soilHydraulicModel = "Van Genuchten";
		String thermalModel = "Johansen";
		
		SoilParametrization soilPar;
		SimpleSoilParametrizationFactory soilParFactory = new SimpleSoilParametrizationFactory();
		soilPar = soilParFactory.createSoilParametrization(soilHydraulicModel,alpha,n,psiE,lambda,rMedian,sigma,thetaR,thetaS,ks);
		
		ThermalSoilParameterization thermalPar;
		SimpleThermalSoilParameterizationFactory thermalParFactory = new SimpleThermalSoilParameterizationFactory();
		thermalPar = thermalParFactory.createThermalSoilParameterization(thermalModel, thetaS, 0.35, 0.4,
				soilPar, "mineral soil", 2.5, 0.5);
		
		//PrintTXT print = new PrintTXT();
		
				
		double myHeatCapacity = thermalPar.heatCapacity(-10);
		System.out.println("myHeatCapacity: " +myHeatCapacity);
		double myThermalConductivity = thermalPar.thermalConductivity(-10);
		System.out.println("myThermalConductivity: " +myThermalConductivity);
		
		//thermalParametrization = thermalPar.thermalSoilPropertiesCurves();
		//print.setValueMatrix(thermalParametrization);
		//print.PrintMatrix("resources/Output", "heatCapacityParametrization"+".csv", soilHydraulicModel+"&"+thermalModel, "Psi, CT, lambda, na, na");
	} 
	

}