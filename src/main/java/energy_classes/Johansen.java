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
package energy_classes;

import richards_classes.SoilParametrization;

/**
 * 
 * This class contains the Johansen's model to compute the thermal conductivity (Johansen, 1975).
 * @author Niccolò Tubini
 *
 */
public class Johansen extends ThermalSoilParametrization {
	
	private double soilBulkDensity; //rho_dry
	private double lambdaDry;
	private double lambdaSaturated;
	private double lambdaG;
	private double lambda0;
	private double quartzFraction;
	private double lambda;
	private String drySoil;
	final double lambdaW = 0.6;
	
	public Johansen(double thetaS, double sandFraction, double clayFraction,
			SoilParametrization soilPar, double soilBulkDensity, String drySoil, double lambda0, double quartzFraction) {
		super(thetaS, sandFraction, clayFraction, soilPar);
		this.soilBulkDensity = soilBulkDensity;
		this.lambda0 = lambda0;
		this.drySoil = drySoil;
		this.quartzFraction = quartzFraction;
		
		if(this.lambda0 > 3 || this.lambda0 < 2){
			throw new IllegalArgumentException( "ERROR: Check the value of the lambda0 parameter \n");
		}
		if(this.quartzFraction > 1 || this.quartzFraction <= 0){
			throw new IllegalArgumentException( "ERROR: Check the value of the quartz fraction \n");
		}
		/*if(this.drySoil.equalsIgnoreCase("unfrozen peat") || this.drySoil.equalsIgnoreCase("mineral soil") || this.drySoil.equalsIgnoreCase("crushed rock")){
			throw new IllegalArgumentException( "ERROR: Check the kind of dry soil you selected \n");
		}
		*/
	}
	
	
	
	private double dryThermalConductivity(String drySoil){
		
		if(drySoil.equalsIgnoreCase("unfrozen peat")){
			lambdaDry = 0.55;
		}
		else if(drySoil.equalsIgnoreCase("mineral soil")){
			lambdaDry = (0.135*soilBulkDensity + 64.7)/(2700-0.947*soilBulkDensity);
		}
		else if(drySoil.equalsIgnoreCase("crushed rock")){
			lambdaDry = 0.39*Math.pow(super.thetaS,-2.2);
		}
		
		return lambdaDry;
	}
	
	
	private double saturatedThermalConductivity(){
		lambdaSaturated = Math.pow(quartzThermalConductivity(), 1-super.thetaS)*Math.pow(lambdaW, super.thetaS);
		
		return lambdaSaturated;
	}
	
	
	private double quartzThermalConductivity(){
		
		if(quartzFraction < 0.2){
			lambdaG = Math.pow(7.7, quartzFraction) * Math.pow(lambda0, 1-quartzFraction);
		}
		else{
			lambdaG = (8.80*super.sandFraction + 2.92*super.clayFraction)/(super.sandFraction + super.clayFraction);
		}
		return lambdaG;
	}

	
	@Override
	public double thermalConductivity(double suction) {
		lambda = saturatedThermalConductivity()*super.soilPar.waterContent(suction)/super.thetaS + (1-super.soilPar.waterContent(suction)/super.thetaS)*dryThermalConductivity(this.drySoil); 
		return lambda;
	}

}
