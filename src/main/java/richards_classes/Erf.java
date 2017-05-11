package richards_classes;

public class Erf {
	//http://introcs.cs.princeton.edu/java/21function/ErrorFunction.java.html
    public double erf(double z) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(z));

        // use Horner's method
        double ans = 1 - t * Math.exp( -z*z   -   1.26551223 +
                                            t * ( 1.00002368 +
                                            t * ( 0.37409196 + 
                                            t * ( 0.09678418 + 
                                            t * (-0.18628806 + 
                                            t * ( 0.27886807 + 
                                            t * (-1.13520398 + 
                                            t * ( 1.48851587 + 
                                            t * (-0.82215223 + 
                                            t * ( 0.17087277))))))))));
        if (z >= 0) return  ans;
        else        return -ans;
    }
}
//public static double erf(double x) {
//// constants
//final double a1 =  0.254829592;
//final double a2 = -0.284496736;
//final double a3 =  1.421413741;
//final double a4 = -1.453152027;
//final double a5 =  1.061405429;
//final double p  =  0.3275911;
//
//// Save the sign of x
//double sign = 1;
//if (x < 0) {
//  sign = -1;
//}
//x = Math.abs(x);
//
//// A&S formula 7.1.26
//double t = 1.0/(1.0 + p*x);
//double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*Math.exp(-x*x);
//
//return sign*y;
//}