package richards_utils;

public class Sequencer {
	/**
	 * This class creates a sequence of values 
	 * @param min value of the sequence
	 * @param max value of the sequence
	 * @param points, number of values in the sequence
	 * @return sequence[] values in the sequence
	 */
	public static double[] seq(double min, double max, int points) {  
	    double[] sequence = new double[points];  
	    for (int i = 0; i < points; i++){  
	        sequence[i] = min + i * (max - min) / (points - 1);  
	    }  
	    return sequence;  
	}	
}