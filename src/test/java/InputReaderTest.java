import richards_utils.Expr;
import richards_utils.TextIO;

import org.junit.Test;

public class InputReaderTest {
	private static boolean interactive = true;
	private static String predefinedFunc = "sin(x)^2+4";
	private static int[] values = {1,2,3,4};
	
	@Test
	public void test(){
	    String line;
	    Expr function = null;
	    double val;

		System.out.println("This program will evaluate a specified function, f(x), at");
		System.out.println("specified values of the variable x.  The definition of f(x)");
		System.out.println("can use the operators +, -, *, /, and ^ as well as mathematical");
		System.out.println("functions such as sin, abs, and ln.");

		if(interactive==true) {

			line = TextIO.getln().trim();

			try {
				function = new Expr(line);
			}
			catch (IllegalArgumentException e) {
				System.out.println("Error! The definition of f(x) is not valid.");
				System.out.println(e.getMessage());
			}
		} else {
			function = new Expr(predefinedFunc);			
		}
		
		for(int i=0;i<values.length;i++) {
			val = function.value(values[i]);
			System.out.println(val);
		}				
		
	} 
}
