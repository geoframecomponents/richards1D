package richards_classes;


public abstract class BoundaryCondition {
	
	protected double bC;
	protected double kP;
	protected double kM;
	protected double spaceDeltaP;
	protected double spaceDeltaM;
	protected double tTimestep;
	protected double delta;
	protected double term;
	
	public abstract double upperDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta);
	
	public abstract double mainDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta);
	
	public abstract double lowerDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta);
	
	public abstract double rightHandSide(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta);
}
