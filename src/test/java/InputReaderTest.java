import richards_utils.Expr;
import richards_utils.TextIO;
import oms3.annotations.*;
import org.junit.Test;
import ic.*;

@Author(
		name="Aaron Iemma",
		org="DICAM - Departement of Environmental and Civil Engineering - Trento, UNITN",
		contact="iemma.ron -at- gmail.com"
	)
	@Description("Testing the reading of IC and boundary for the RIchards1D problem")
	@Keywords("R,Renjin")

public class InputReaderTest {
	@In public boolean interactive = true;
	@In public String givenfunction = "sin(x)^2+4";
	private static int[] values = {1,2,3,4};
	
	@Test
	public void testfunctionparsing(){
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
			function = new Expr(givenfunction);			
		}
		
		for(int i=0;i<values.length;i++) {
			val = function.value(values[i]);
			System.out.println(val);
		}				
		
	}
	
	@Test
	public void testICReading() {
		Domain1D mydomain = new Domain1D();

		mydomain.performReadDomainFile("/home/cifciaf/Desktop/testRead/ictest");
		mydomain.readicfromfile.parse();
		mydomain.readicfromfile.show();		
	}
}
