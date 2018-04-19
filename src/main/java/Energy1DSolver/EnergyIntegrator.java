package Energy1DSolver;

import richards_classes.BoundaryCondition;
import richards_classes.JordanDecomposition;
import richards_classes.PrintTXT;
import richards_classes.SoilParametrization;
import richards_classes.Thomas;
import richards_classes.TotalDepth;
import thermalParameterization.SimpleThermalSoilParameterizationFactory;
import thermalParameterization.ThermalSoilParameterization;

public abstract class EnergyIntegrator {
	
	protected int NUM_CONTROL_VOLUMES;
	protected double timeDelta;

	protected double temperatureTopBC;
	protected double temperatureBottomBC;
	protected double lambdaP;
	protected double lambdaM;
	protected double internalEnergy;
	protected double internalEnergyNew;
	protected double totalEnergy;
	protected double errorEnergy;
	protected double totalDiffusionFlux;
	protected double totalConvectionFlux;
	protected double Peclet;
	
	protected double[] diffusionFluxes;
	protected double[] convectionFluxes;
	protected double[] Peclets;
	
	protected double[] temperatures;
	protected double[] psis;
	protected double[] psisNew;
	protected double[] velocities;
	protected double[] spaceDelta;
	protected double[] dx;
	
	protected double[] mainDiagonal;
	protected double[] upperDiagonal;
	protected double[] lowerDiagonal;
	protected double[] rhss;
	
	protected ThermalSoilParameterization thermalSoilPar;
	protected BoundaryCondition tempTopBoundaryCondition;
	protected BoundaryCondition tempBottomBoundaryCondition;
	protected BoundaryCondition tempConvectionTopBoundaryCondition;
	protected BoundaryCondition tempConvectionBottomBoundaryCondition;
	protected Thomas tempThomasAlg;

	
	public EnergyIntegrator() {};
	
	
	
	public EnergyIntegrator(ThermalSoilParameterization thermalSoilPar, int NUM_CONTROL_VOLUMES, double[] dx, double[] spaceDelta,
			BoundaryCondition tempTopBoundaryCondition, BoundaryCondition tempBottomBoundaryCondition,
			BoundaryCondition tempConvectionTopBoundaryCondition, BoundaryCondition tempConvectionBottomBoundaryCondition) {
		
		this.thermalSoilPar = thermalSoilPar;
		this.NUM_CONTROL_VOLUMES = NUM_CONTROL_VOLUMES-1;
		this.dx = dx;
		this.spaceDelta = spaceDelta;
		// Forse e meglio creali direttamente qui come viene fatto per thermalSoilPar
		this.tempTopBoundaryCondition = tempTopBoundaryCondition;
		this.tempBottomBoundaryCondition = tempBottomBoundaryCondition;
		this.tempConvectionTopBoundaryCondition = tempConvectionTopBoundaryCondition;
		this.tempConvectionBottomBoundaryCondition = tempConvectionBottomBoundaryCondition;
		
		tempThomasAlg = new Thomas();
				
		upperDiagonal = new double[this.NUM_CONTROL_VOLUMES];
		mainDiagonal = new double[this.NUM_CONTROL_VOLUMES];
		lowerDiagonal = new double[this.NUM_CONTROL_VOLUMES];
		rhss = new double[this.NUM_CONTROL_VOLUMES];
		
		diffusionFluxes = new double[this.NUM_CONTROL_VOLUMES+1];
		convectionFluxes = new double[this.NUM_CONTROL_VOLUMES+1];
		Peclets = new double[this.NUM_CONTROL_VOLUMES+1];
		
	}
	
	
	
	
	public void set(double[] psis, double[] psisNew, double[] velocities, double[] temperatures, double timeDelta, double temperatureBottomBC, double temperatureTopBC) {
		
		this.psis = psis;
		this.psisNew = psisNew;
		this.temperatures = temperatures;
		this.velocities = velocities;
		this.timeDelta = timeDelta;
		this.temperatureTopBC = temperatureTopBC;
		this.temperatureBottomBC = temperatureBottomBC;

		
	}
	
	
	
	
	public double computeTotalInternalEnergy(double[] psis, double[] temperatures) {
		totalEnergy = 0.0;
		
		for(int i=0; i< NUM_CONTROL_VOLUMES; i++ ) {
			
			totalEnergy += this.thermalSoilPar.heatCapacity(psis[i])*temperatures[i]*this.dx[i];
			
		}
		
		return totalEnergy;
		
	}
	
	
	
	public double[] computeDiffusionFluxes(double[] psis, double[] temperatures, double temperatureBottomBC, double temperatureTopBC) {
		
		for(int i=0; i< NUM_CONTROL_VOLUMES+1; i++ ) {
			
			if(i==0) {
				diffusionFluxes[i] = thermalSoilPar.thermalConductivity(psis[i])*(temperatures[i]-temperatureBottomBC)/this.spaceDelta[i];

				
			} else if(i==NUM_CONTROL_VOLUMES) {
				diffusionFluxes[i] = thermalSoilPar.thermalConductivity(psis[i-1])*(temperatureTopBC-temperatures[i-1])/this.spaceDelta[i];

				
			} else {
				
				diffusionFluxes[i] = 0.5*(thermalSoilPar.thermalConductivity(psis[i-1])+thermalSoilPar.thermalConductivity(psis[i]))*(temperatures[i]-temperatures[i-1])/this.spaceDelta[i];

			}
		}
				
		return diffusionFluxes;
	}
	
	
	
	public double[] computeConvectionFluxes(double[] velocities, double[] temperatures, double temperatureBottomBC, double temperatureTopBC) {
		
		for(int i=0; i< NUM_CONTROL_VOLUMES+1; i++ ) {
			
			if(i==0) {
				convectionFluxes[i] = 4188000*0.5*( velocities[i]*(temperatures[i]+temperatureBottomBC) - Math.abs(velocities[i])*(temperatures[i]-temperatureBottomBC) );
				
			} else if(i==NUM_CONTROL_VOLUMES) {
				convectionFluxes[i] = 4188000*0.5*( velocities[i]*(temperatureTopBC+temperatures[i-1]) - Math.abs(velocities[i])*(temperatureTopBC-temperatures[i-1]) );

				
			} else {
				
				convectionFluxes[i] = 4188000*0.5*( velocities[i]*(temperatureTopBC+temperatures[i-1]) - Math.abs(velocities[i])*(temperatureTopBC-temperatures[i-1]) );

			}
		}
				
		return convectionFluxes;
	}
	
	
	
	
	public double[] computePeclets() {
		
		for(int i=0; i< NUM_CONTROL_VOLUMES+1; i++ ) {
			//Peclets[i] = Math.abs(this.convectionFluxes[i])/Math.abs(this.diffusionFluxes[i]);
			if(i==0) {
				Peclets[i] = ( Math.abs(this.velocities[i])*this.spaceDelta[i])/( thermalSoilPar.thermalConductivity(this.psis[i])/4188000);
				
			} else if(i==NUM_CONTROL_VOLUMES) {
				Peclets[i] = ( Math.abs(this.velocities[i])*this.spaceDelta[i])/( thermalSoilPar.thermalConductivity(this.psis[i])/4188000);
			} else {
				Peclets[i] = ( Math.abs(this.velocities[i])*this.spaceDelta[i])/( 0.5*( thermalSoilPar.thermalConductivity(this.psis[i]) +thermalSoilPar.thermalConductivity(this.psis[i+1]) )/4188000);
			}
		}
		
		return Peclets;
	}
	
	
	public abstract double[] solver1D();

}
