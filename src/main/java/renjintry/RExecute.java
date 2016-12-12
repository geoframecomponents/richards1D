package renjintry;

// R script
import javax.script.*;
// OMS3 packages
import oms3.annotations.*;
// Others
import static java.lang.Math.*;

@Status(Status.CERTIFIED)

@Author(
	name="Aaron Iemma",
	org="DICAM - Departement of Environmental and Civil Engineering - Trento, UNITN",
	contact="iemma.ron@gmail.com"
)
@Description("A class to execute a given R script")
@Keywords("R,Renjin")

public class RExecute {
	/**
	 * Executes a given R script
	 *
	 * @param script name
	 * @see script output
	 */
    @In public String scriptName;

    @Execute
    public void execute() throws Exception {
        ScriptEngineManager manager = new ScriptEngineManager();
        ScriptEngine renjin = manager.getEngineByName("Renjin");
        // check if the engine has loaded correctly:
        if(renjin == null) {
            throw new RuntimeException("Renjin Script Engine not found on the classpath.");
        }

		renjin.eval(new java.io.FileReader(scriptName));          	
    }
}

