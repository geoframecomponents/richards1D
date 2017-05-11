package richards_classes;

public class DomainDiscretization {
	/**
	 * Returns a vector (array of doubles) of equally spaced numbers between a given 
	 * lower and upper bound
	 * <p>
	 * The lower and upper bounds need not be in growing order: in 
	 * case min>max, the method generates a sequence of decreasing 
	 * real numbers
	 *
	 * @param  max  	a real number, the value of the lower bound of the sequence
	 * @param  max  	a real number, the value of the upper bound of the sequence
	 * @param  points  	an integer, the number of equally spaced points in the sequence 
	 * @return 			a vector of equally spaced real numbers
	 */	
	public static double[] seq(double min, double max, int points) {  
	    double[] sequence = new double[points];  
	    for (int i = 0; i < points; i++){  
	        sequence[i] = min + i * (max - min) / (points - 1);  
	    }  
	    return sequence;  
	}	
}