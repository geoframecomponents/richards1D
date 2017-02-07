package richards_utils;

public class DomainRectangularGridRegularDiscretization {
	/**
	 * Representation by its centres of a regular gridded (linear, quadratic or cuboid) domain
	 */	
	
	private double[] lengths 	= new double[3];
	private int[] numpoints		= new int[3];

	/**
	 * This class creates a matrix of centres 
	 * @param min value of the sequence
	 * @param max value of the sequence
	 * @param points, number of values in the sequence
	 * @return sequence[] values in the sequence
	 */
	public double[][][][] makeGrid(double[] lengths, int[] numpoints) {
		
		this.lengths = lengths;
		this.numpoints = numpoints;
		double[][][][] centres = new double[numpoints[0]][numpoints[1]][numpoints[2]][3];
		
		double[] sequence_i = new double[numpoints[0]];
		double[] sequence_j = new double[numpoints[1]];
		double[] sequence_k = new double[numpoints[2]];
		
		Sequencer sequencer = new Sequencer();
		
		sequence_i = sequencer.seq(0, lengths[0], numpoints[0]);
		sequence_j = sequencer.seq(0, lengths[1], numpoints[1]);
		sequence_k = sequencer.seq(0, lengths[2], numpoints[2]);
		
		for(int i=0;i<numpoints[0];i++) {
			for(int j=0;j<numpoints[1];j++) {
				for(int k=0;k<numpoints[2];k++) {
					centres[i][j][k][0] = sequence_i[i];
					centres[i][j][k][1] = sequence_j[j];
					centres[i][j][k][2] = sequence_k[k];
				}
			}
		}
		
		return centres;
	}
	
	
}
