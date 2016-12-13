import org.junit.Test;
import richards_classes.RichardsSolver;
import richards_classes.Richards;

public class RichardsTest {

	@Test
	public void testTest() {

		RichardsSolver richardssolver = new RichardsSolver();
 
		Richards richards = richardssolver.solveRichards("1D");
		
		System.out.println("My GOD this was hard!");

	}
}
