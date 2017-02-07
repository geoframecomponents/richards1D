package richards_utils;

public class DomainRectangularGridRegularDiscretization {
	/**
	 * Representation by its centres of a regular gridded (linear, quadratic or cuboid) domain
	 */	
	
	private double[] lenghts 	= new double[3];
	private int[] numpoints		= new int[3];
	private double[][][] centres = new double[numpoints[0]][numpoints[1]][numpoints[2]];

	/**
	 * This class creates a matrix of centres 
	 * @param min value of the sequence
	 * @param max value of the sequence
	 * @param points, number of values in the sequence
	 * @return sequence[] values in the sequence
	 */
	public double[][][] makeGrid(double[] lengths, int[] points, double[][][] centres) {
		
		this.lengths = lengths;
		this.points = points;
		this.centres = centres;
		
		for(int i = 0; i<3;i++) {
			Sequencer sequencer = new Sequencer();
			centres[i] = sequencer.seq(0,lengths[i],points[i])
		}
		return centres;
	}
	
	
}
