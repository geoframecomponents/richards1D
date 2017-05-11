package richards_classes;

public class HydrostaticSuction extends SuctionInitialCondition{
	/**
	 * This class set the initial condition for suction as hydrostatic
	 * @param space is the space coordinate stored as the center of the control volume
	 * @return suction 
	 */
	
	public double IC(double space){
		
		this.suction = space-2;
		// Initial condition like Dumbser (Matlab code)
		//this.suction = -space;
		return this.suction;
	}

}
