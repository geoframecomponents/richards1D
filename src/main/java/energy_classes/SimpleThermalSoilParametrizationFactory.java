package energy_classes;

import richards_classes.SoilParametrization;

public class SimpleThermalSoilParametrizationFactory {
	
	public ThermalSoilParametrization createThermalSoilParametrization (String type, double thetaS, double sandFraction, double clayFraction,
			SoilParametrization soilPar, double soilBulkDensity, String drySoil, double lambda0, double quartzFraction) {

		ThermalSoilParametrization thermalPar = null;
		if(type.equalsIgnoreCase("Johansen")){
			thermalPar = new Johansen( thetaS, sandFraction, clayFraction,
					 soilPar, soilBulkDensity,  drySoil, lambda0, quartzFraction);
		}
		
		return thermalPar;
	}	
}
