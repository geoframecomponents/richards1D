package richards;

/*
 * GNU GPL v3 License
 *
 * Copyright 2015 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
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