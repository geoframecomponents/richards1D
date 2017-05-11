package richards_classes;

public abstract class SuctionBoundaryCondition {
	
	protected double suction;
	
	public abstract double BC(double time);

}
