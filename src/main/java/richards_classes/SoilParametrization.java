package richards_classes;

public abstract class SoilParametrization {
	protected double psiStar; // suction value at which the derivative of hydraulic capacity is null
	protected double thetaR;  // residual water content 
	protected double thetaS;  // water content at saturation
	protected double kappaSaturation; // hydraulic conductivity at saturation
	
	protected double theta; // water content
	protected double dTheta;
	protected double saturationDegree; // hydraulic conductivity at saturation
	protected double kappa; // hydraulic conductivity
	
	public double getPsiStar(){
		return this.psiStar;
	}
	
	public abstract double waterContent(double suction);
	
	public abstract double dWaterContent(double suction);
	
	public abstract double hydraulicConductivity(double suction);
	
	
}
