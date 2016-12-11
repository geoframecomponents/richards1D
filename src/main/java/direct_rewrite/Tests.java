import java.lang.Math;
import java.util.Arrays;



// Used to test some dubious things in Richards code...
public class Tests {
	static double number	= 5.0;
	static double exponent 	= 2.0;	
	static double 	space_bottom		= 0.0;
	static double 	space_top			= 2.0;
	static int 	NUM_CONTROL_VOLUMES	= 100; 
	static double 	space_delta			= (space_top - space_bottom) / NUM_CONTROL_VOLUMES; // delta

	public static void main(String[] args) {
		double exponent_doubled;
		double result;

		result = Math.pow(number,exponent);
		print("Number " + number + " elevated at " + exponent + " gives " + result);

		exponent = 0.5;
		number   = 4.0;
		result = Math.pow(number,exponent);
		print("Number " + number + " elevated at " + exponent + " gives " + result);

		number = -4.0;
		result = Math.abs(number);
		print("Number " + number + " absolute value is " + result);		

		print(seq(space_bottom + space_delta / 2,space_top - space_delta / 2,NUM_CONTROL_VOLUMES));
		print(seq(space_bottom + space_delta / 2,space_top - space_delta / 2,NUM_CONTROL_VOLUMES).length);


	}
	private static void print(String arg) {
		System.out.println(arg);
	}
	private static void print(double arg) {
		System.out.println(arg);
	}
	private static void print(double[] arg) {
		System.out.println(java.util.Arrays.toString(arg));
	}

	public static double[] seq(double min, double max, int points) {  
	    double[] sequence = new double[points];
	    for (int i = 0; i < points; i++){  
	        sequence[i] = min + i * (max - min) / (points - 1);  
	    }  
	    return sequence;  
	}	


}