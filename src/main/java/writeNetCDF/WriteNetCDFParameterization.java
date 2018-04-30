/*
 * GNU GPL v3 License
 *
 * Copyright 2016 Marialaura Bancheri
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

package writeNetCDF;

import java.io.IOException;
import java.util.LinkedHashMap;
import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Keywords;
import oms3.annotations.License;
import ucar.ma2.ArrayDouble;
import ucar.ma2.DataType;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFileWriter;
import ucar.nc2.Variable;

@Description("This class writes in a NetCDF hydraulic parameterizations: theta(psi), k(psi), moisture capacity(psi).")
@Documentation("")
@Author(name = "Niccolo' Tubini, Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Hydrology, Richards, Infiltration")
//@Label(JGTConstants.HYDROGEOMORPHOLOGY)
//@Name("shortradbal")
//@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")


public class WriteNetCDFParameterization {

	int NLVL;
	double[] myTempVariable; 

	// Create the file.
	String filename;
	NetcdfFileWriter dataFile;

	
	
	public void writeNetCDF(LinkedHashMap<String,double[]> myVariables,String fileName,String model) {

		final int NLVL = myVariables.get("Psi [m]").length;
		//System.out.println(NLVL);
		// human readable date will be converted in unix format, the format will be an input and it has to be consistent with that used in OMS        

		// Create the file.
		//String filename = fileName;
		NetcdfFileWriter dataFile = null;

		try {
			// Create new netcdf-3 file with the given filename
			dataFile = NetcdfFileWriter.createNew(NetcdfFileWriter.Version.netcdf3, fileName);
			// add a general attribute describing the problem and containing other relevant information for the user
			dataFile.addGroupAttribute(null, new Attribute("Description of the problem",model));

			//add dimensions  where time dimension is unlimit
			// in 1D case dimension are time and the depth
			Dimension lvlDim = dataFile.addDimension(null, "psi", NLVL);

			// Define the coordinate variables.
			Variable psiVar = dataFile.addVariable(null, "psi", DataType.DOUBLE, "psi");

			// Define units attributes for data variables.
			dataFile.addVariableAttribute(psiVar, new Attribute("units", "m"));
			dataFile.addVariableAttribute(psiVar, new Attribute("long_name", "Water suction"));

			// Define the netCDF variables for the psi and theta data.
			String dims = "psi";

			Variable thetaVar = dataFile.addVariable(null, "theta", DataType.DOUBLE, dims);
			Variable dThetaVar = dataFile.addVariable(null, "dTheta", DataType.DOUBLE, dims);
			Variable hydraulicConductivityVar = dataFile.addVariable(null, "hydraulic conductivity", DataType.DOUBLE, dims);

			// Define units attributes for data variables.
			dataFile.addVariableAttribute(thetaVar, new Attribute("units", "-"));
			dataFile.addVariableAttribute(thetaVar, new Attribute("long_name", "Adimensional water content"));
			dataFile.addVariableAttribute(dThetaVar, new Attribute("units", "1/m"));
			dataFile.addVariableAttribute(dThetaVar, new Attribute("long_name", "Soil moisture capacity"));
			dataFile.addVariableAttribute(hydraulicConductivityVar, new Attribute("units", "m/s"));
			dataFile.addVariableAttribute(hydraulicConductivityVar, new Attribute("long_name", "Hydraulic conductivity (Mualem model)"));

			// These data are those created by bufferWriter class. If this wasn't an example program, we
			// would have some real data to write for example, model output.
			// times variable is filled later
			ArrayDouble.D1 psis = new ArrayDouble.D1(lvlDim.getLength());

			// These data are those created by bufferWriter class. This will write our hydraulic head (psi) and
			// adimensional water content (theta) data
			ArrayDouble.D1 dataTheta = new ArrayDouble.D1(lvlDim.getLength());
			ArrayDouble.D1 dataDTheta = new ArrayDouble.D1(lvlDim.getLength());
			ArrayDouble.D1 dataHydraulicConductivity =  new ArrayDouble.D1(lvlDim.getLength());



			myTempVariable =  myVariables.get("Psi [m]");
			for (int lvl = 0; lvl < NLVL; lvl++) {
				//System.out.println(myTempVariable[lvl]);
				psis.set(lvl, myTempVariable[lvl]);

			}

			myTempVariable =  myVariables.get("Theta [-]");
			for (int lvl = 0; lvl < NLVL; lvl++) {
				//System.out.println(myTempVariable[lvl]);
				dataTheta.set(lvl, myTempVariable[lvl]);

			}


			myTempVariable =  myVariables.get("dTheta [1/m]");
			for (int lvl = 0; lvl < NLVL; lvl++) {

				dataDTheta.set(lvl, myTempVariable[lvl]);

			}

			myTempVariable =  myVariables.get("K [m/s]");
			for (int lvl = 0; lvl < NLVL; lvl++) {

				dataHydraulicConductivity.set(lvl, myTempVariable[lvl]);

			}


			//Create the file. At this point the (empty) file will be written to disk
			dataFile.create();

			// A newly created Java integer array to be initialized to zeros.
			int[] origin = new int[1];

			dataFile.write(psiVar, psis);
			dataFile.write(thetaVar,dataTheta);
			dataFile.write(dThetaVar, dataDTheta);
			dataFile.write(hydraulicConductivityVar, dataHydraulicConductivity);

		} catch (IOException e) {
			e.printStackTrace(System.err);

		} catch (InvalidRangeException e) {
			e.printStackTrace(System.err);

		} finally {
			if (dataFile != null)
				try {
					dataFile.close();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
		}

		System.out.println("*** SUCCESS writing output file, " + fileName);

	}


}
