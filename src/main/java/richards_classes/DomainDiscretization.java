package richards_classes;

public class DomainDiscretization {
	/**
	 * This class discretizes the 1D domain 
	 * @param min 
	 * @param max
	 * @param points number of finite volumes
	 * @return
	 */
	public static double[] seq(double min, double max, int points) {  
	    double[] sequence = new double[points];  
	    for (int i = 0; i < points; i++){  
	        sequence[i] = min + i * (max - min) / (points - 1);  
	    }  
	    return sequence;  
	}	
}
