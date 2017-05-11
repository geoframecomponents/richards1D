package richards_classes;

public class SuctionBCConstant extends SuctionBoundaryCondition {
	/**
	 * This class set the boundary condition for suction as constant in time
	 * @param time
	 * @return suction 
	 */
	public double BC(double time){
		this.suction = 0.0; // set the value for the boundary condition
		return this.suction;
	}

}
