package richards_classes;


public abstract class BoundaryCondition {
	
	protected double bC;
	protected double kP;
	protected double kM;
	protected double gridvarP;
	protected double gridvarM;
	protected double gridvarsqP;
	protected double gridvarsqM;
	protected double term;
	
	public abstract double upperDiagonal(double bC, double kP, double kM, double gridvarsqP, double gridvarsqM, double gridvarP, double gridvarM);
	
	public abstract double mainDiagonal(double bC, double kP, double kM, double gridvarsqP, double gridvarsqM, double gridvarP, double gridvarM);
	
	public abstract double lowerDiagonal(double bC, double kP, double kM, double gridvarsqP, double gridvarsqM, double gridvarP, double gridvarM);
	
	public abstract double rightHandSide(double bC, double kP, double kM, double gridvarsqP, double gridvarsqM, double gridvarP, double gridvarM);
}
