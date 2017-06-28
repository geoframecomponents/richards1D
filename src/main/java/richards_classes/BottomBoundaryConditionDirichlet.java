package richards_classes;

public class BottomBoundaryConditionDirichlet extends BoundaryCondition {

	public double upperDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.tTimestep = tTimestep;
		this.delta = delta;
		
		this.term =-this.kP * this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)*1/this.spaceDeltaP*1/Math.pow(Math.cos(this.delta),2);
		
		return term;
	}
	
	public double mainDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.tTimestep = tTimestep;
		this.delta = delta;
		
		this.term = this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)*1/this.spaceDeltaM*1/Math.pow(Math.cos(this.delta),2)* this.kM + this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)*1/this.spaceDeltaP*1/Math.pow(Math.cos(this.delta),2)*this.kP;

		return term;

	}
	
	public double lowerDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.tTimestep = tTimestep;
		this.delta = delta;
		
		this.term = 0;

		return term;

	}

	public double rightHandSide(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double tTimestep, double delta) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.spaceDeltaP = spaceDeltaP;
		this.spaceDeltaM = spaceDeltaM;
		this.tTimestep = tTimestep;
		this.delta = delta;
		
		this.term = this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)* this.kP - this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)*this.kM + this.kM * this.tTimestep/(this.spaceDeltaM+this.spaceDeltaP/2)*1/this.spaceDeltaM*1/Math.pow(Math.cos(this.delta),2)* this.bC;

		return term;

	}
}
