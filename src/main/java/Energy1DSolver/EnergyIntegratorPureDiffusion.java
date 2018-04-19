package Energy1DSolver;

import richards_classes.BoundaryCondition;
import richards_classes.SoilParametrization;
import thermalParameterization.ThermalSoilParameterization;

public class EnergyIntegratorPureDiffusion extends EnergyIntegrator {
	
	public EnergyIntegratorPureDiffusion (ThermalSoilParameterization thermalSoilPar, int NUM_CONTROL_VOLUMES, double[] dx, double[] spaceDelta,
			BoundaryCondition tempTopBoundaryCondition, BoundaryCondition tempBottomBoundaryCondition,
			BoundaryCondition tempConvectionTopBoundaryCondition, BoundaryCondition tempConvectionBottomBoundaryCondition) {
	
		super( thermalSoilPar, NUM_CONTROL_VOLUMES, dx, spaceDelta, tempTopBoundaryCondition, tempBottomBoundaryCondition,
				tempConvectionTopBoundaryCondition, tempConvectionBottomBoundaryCondition);
		System.out.println("Energy equation with pure conduction.");
	}

	
	
	
	@Override
	public double[] solver1D() {

		internalEnergy = super.computeTotalInternalEnergy(psis, temperatures);
				
		for(int i=0; i<NUM_CONTROL_VOLUMES; i++) {
			if(i==0) {
			
				lambdaP = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i+1]) );
				lambdaM = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i]) );
				
				upperDiagonal[i] = tempBottomBoundaryCondition.upperDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999);;
				mainDiagonal[i] = super.thermalSoilPar.heatCapacity(super.psis[i])*dx[i] + tempBottomBoundaryCondition.mainDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999);;
				lowerDiagonal[i] = tempBottomBoundaryCondition.lowerDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999);
				rhss[i] = super.thermalSoilPar.heatCapacity(super.psis[i])*temperatures[i]*dx[i] + tempBottomBoundaryCondition.rightHandSide(temperatureBottomBC, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999);
			
				
			} else if(i==NUM_CONTROL_VOLUMES-1) {
				
				lambdaP = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i]) );
				lambdaM = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i-1]) );
				
				upperDiagonal[i] = tempTopBoundaryCondition.upperDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999);;
				mainDiagonal[i] = super.thermalSoilPar.heatCapacity(super.psis[i])*dx[i] + tempTopBoundaryCondition.mainDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999);;
				lowerDiagonal[i] = tempTopBoundaryCondition.lowerDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999);
				rhss[i] = super.thermalSoilPar.heatCapacity(super.psis[i])*temperatures[i]*dx[i] + tempTopBoundaryCondition.rightHandSide(temperatureTopBC, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i],timeDelta,-999);
			
			} else {
							
				lambdaP = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i+1]) );
				lambdaM = 0.5*( super.thermalSoilPar.thermalConductivity(super.psis[i]) + super.thermalSoilPar.thermalConductivity(super.psis[i-1]) );
				
				upperDiagonal[i] = -timeDelta/spaceDelta[i+1]*lambdaP;
				mainDiagonal[i] = super.thermalSoilPar.heatCapacity(super.psis[i])*dx[i] + timeDelta/spaceDelta[i]*lambdaM + timeDelta/spaceDelta[i+1]*lambdaP;
				lowerDiagonal[i] = -timeDelta/spaceDelta[i]*lambdaM;
				rhss[i] = super.thermalSoilPar.heatCapacity(super.psis[i])*temperatures[i]*dx[i];
				
			}
			
		}
		
		/*for(int i=0; i<NUM_CONTROL_VOLUMES; i++) {
			System.out.println(super.thermalSoilPar.heatCapacity(super.psis[i]));
		}*/
				
		tempThomasAlg.set(upperDiagonal, mainDiagonal, lowerDiagonal, rhss);
		super.temperatures = tempThomasAlg.solver();
		
		internalEnergyNew = super.computeTotalInternalEnergy(psis, temperatures);
		
		errorEnergy = internalEnergyNew - internalEnergy + timeDelta*( -super.thermalSoilPar.thermalConductivity(super.psis[NUM_CONTROL_VOLUMES-1])*(temperatureTopBC-temperatures[NUM_CONTROL_VOLUMES-1])/spaceDelta[NUM_CONTROL_VOLUMES] 
																				+ super.thermalSoilPar.thermalConductivity(super.psis[0])*(temperatures[0]-temperatureBottomBC)/spaceDelta[0]);
		
		
		System.out.println("     errorEnergy: "+errorEnergy);
		return super.temperatures;
	}

}
