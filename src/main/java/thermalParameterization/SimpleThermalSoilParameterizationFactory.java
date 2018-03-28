package thermalParameterization;

import richards_classes.SoilParametrization;

public class SimpleThermalSoilParameterizationFactory {
	
	public ThermalSoilParameterization createThermalSoilParameterization (String type, double thetaS, double sandFraction, double clayFraction,
			SoilParametrization soilPar, String drySoil, double lambda0, double quartzFraction) {

		ThermalSoilParameterization thermalPar = null;
		if(type.equalsIgnoreCase("Johansen")){
			thermalPar = new Johansen( thetaS, sandFraction, clayFraction,
					 soilPar, drySoil, lambda0, quartzFraction);
		}
		
		return thermalPar;
	}	
}