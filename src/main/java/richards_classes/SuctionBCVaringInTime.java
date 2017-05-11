package richards_classes;

public class SuctionBCVaringInTime extends SuctionBoundaryCondition {
	/**
	 * This class set the boundary condition for suction as a function of time
	 * @param time
	 * @return suction 
	 */
	public double BC(double time){
		
		if(time <= Math.pow(10,5)) {
	        this.suction = -0.05 + 0.03 * Math.sin(2 * Math.PI * time/Math.pow(10,5));
	    } else if(time > Math.pow(10,5) && time <= 1.8 * Math.pow(10,5)) {
	    	this.suction = +0.1;
	    }
	    else {
	    	this.suction = -0.05 + 2952.45 * Math.exp(-time / 18204.8);
	    }
		
		return this.suction;
	}		
}
