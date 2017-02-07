package renjintry;

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

// R script
import javax.script.*;
// OMS3 packages
import oms3.annotations.*;
// Others
import static java.lang.Math.*;

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
    @In private String scriptName;

    @Execute
    public void execute(String scriptName) throws Exception {
        ScriptEngineManager manager = new ScriptEngineManager();
        ScriptEngine renjin = manager.getEngineByName("Renjin");
        // check if the engine has loaded correctly:
        if(renjin == null) {
            throw new RuntimeException("Renjin Script Engine not found on the classpath.");
        }

		renjin.eval(new java.io.FileReader(this.scriptName));          	
    }
}

