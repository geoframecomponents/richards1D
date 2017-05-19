package richards_classes;

public class TopBoundaryConditionNeumann extends BoundaryCondition{

	public double upperDiagonal(double bC, double kP, double kM, double gridvarsqP, double gridvarsqM, double gridvarP, double gridvarM) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.gridvarP = gridvarP;
		this.gridvarM = gridvarM;
		this.gridvarsqP = gridvarsqP;
		this.gridvarsqM = gridvarsqM;
		this.term = 0;

		return term;
	}
	
	public double mainDiagonal(double bC, double kP, double kM, double gridvarsqP, double gridvarsqM, double gridvarP, double gridvarM) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.gridvarP = gridvarP;
		this.gridvarM = gridvarM;
		this.gridvarsqP = gridvarsqP;
		this.gridvarsqM = gridvarsqM;
		
		this.term = this.gridvarsqM * (this.kM);

		return term;

	}
	
	public double lowerDiagonal(double bC, double kP, double kM, double gridvarsqP, double gridvarsqM, double gridvarP, double gridvarM) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.gridvarP = gridvarP;
		this.gridvarM = gridvarM;
		this.gridvarsqP = gridvarsqP;
		this.gridvarsqM = gridvarsqM;
		
		this.term = -this.kM * this.gridvarsqM;

		return term;

	}

	public double rightHandSide(double bC, double kP, double kM, double gridvarsqP, double gridvarsqM, double gridvarP, double gridvarM) {
		this.bC = bC;
		this.kP = kP;
		this.kM = kM;
		this.gridvarP = gridvarP;
		this.gridvarM = gridvarM;
		this.gridvarsqP = gridvarsqP;
		this.gridvarsqM = gridvarsqM;
		
		this.term = this.gridvarM * (-this.kM) + this.gridvarP * this.bC;

		return term;

	}
}
