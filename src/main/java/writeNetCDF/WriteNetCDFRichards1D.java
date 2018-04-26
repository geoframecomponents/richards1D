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
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Execute;
import oms3.annotations.Finalize;
import oms3.annotations.In;
import oms3.annotations.Initialize;
import oms3.annotations.Keywords;
import oms3.annotations.License;
import oms3.annotations.Unit;
import ucar.ma2.Array;
import ucar.ma2.ArrayDouble;
import ucar.ma2.ArrayFloat;
import ucar.ma2.DataType;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFileWriter;
import ucar.nc2.Variable;

@Description("Buffer for 1D simulation in which there is one spatial coordinate and time.")
@Documentation("")
@Author(name = "Niccolo' Tubini, Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Hydrology, Richards, Infiltration")
//@Label(JGTConstants.HYDROGEOMORPHOLOGY)
//@Name("shortradbal")
//@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")


public class WriteNetCDFRichards1D {
	
	@Description()
	@In
	@Unit ()
	public LinkedHashMap<String,ArrayList<double[]>> myVariables; // consider the opportunity to save varibale as float instead of double
	
	@Description()
	@In
	@Unit ()
	public double[] mySpatialCoordinate;
	
	@Description()
	@In
	@Unit ()
	public String fileName;
	
	@Description("Brief descritpion of the problem")
	@In
	@Unit ()
	public String briefDescritpion;
	
	int NLVL;
	int NREC;
	// human readable date will be converted in unix format, the format will be an input and it has to be consistent with that used in OMS
	DateFormat dateFormat;
	Date date = null;
	long unixTime;
	double[] myTempVariable; 
	Iterator it;
    

	// Create the file.
	String filename;
	NetcdfFileWriter dataFile;

	@Execute
	public void writeNetCDF() {
		
		final int NLVL = mySpatialCoordinate.length;
		final int NREC = myVariables.keySet().size();
		// human readable date will be converted in unix format, the format will be an input and it has to be consistent with that used in OMS
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm");
		Date date = null;
		long unixTime;
		double[] myTempVariable; 
		//Iterator it;
        

		// Create the file.
		//String filename = fileName;
		NetcdfFileWriter dataFile = null;

		try {
			// Create new netcdf-3 file with the given filename
			dataFile = NetcdfFileWriter.createNew(NetcdfFileWriter.Version.netcdf3, fileName);
			// add a general attribute describing the problem and containing other relevant information for the user
			dataFile.addGroupAttribute(null, new Attribute("Description of the problem",briefDescritpion));
			
			//add dimensions  where time dimension is unlimit
			// in 1D case dimension are time and the depth
			Dimension lvlDim = dataFile.addDimension(null, "depth", NLVL);
			Dimension timeDim = dataFile.addUnlimitedDimension("time");

			// Define the coordinate variables.
			Variable depthVar = dataFile.addVariable(null, "depth", DataType.DOUBLE, "depth");
			Variable timeVar = dataFile.addVariable(null, "time", DataType.INT, "time");

			// Define units attributes for data variables.
			dataFile.addVariableAttribute(depthVar, new Attribute("units", "m"));
			dataFile.addVariableAttribute(depthVar, new Attribute("long_name", "Soil depth"));
			dataFile.addVariableAttribute(timeVar, new Attribute("units","unix convention"));

			// Define the netCDF variables for the psi and theta data.
			String dims = "time depth";

			Variable psiVar = dataFile.addVariable(null, "psi", DataType.DOUBLE, dims);
			Variable waterHeightVar = dataFile.addVariable(null, "water heigth", DataType.DOUBLE, dims);
			Variable errorVar = dataFile.addVariable(null, "error", DataType.DOUBLE, "time");
			
			// Define units attributes for data variables.
			dataFile.addVariableAttribute(psiVar, new Attribute("units", "m"));
			dataFile.addVariableAttribute(psiVar, new Attribute("long_name", "Hydraulic head"));
			dataFile.addVariableAttribute(waterHeightVar, new Attribute("units", "m"));
			dataFile.addVariableAttribute(waterHeightVar, new Attribute("long_name", "water height"));
			dataFile.addVariableAttribute(errorVar, new Attribute("units", "m"));
			dataFile.addVariableAttribute(errorVar, new Attribute("long_name", "volume error at each time step"));
			
			// These data are those created by bufferWriter class. If this wasn't an example program, we
			// would have some real data to write for example, model output.
			// times variable is filled later
			ArrayDouble.D1 depths = new ArrayDouble.D1(lvlDim.getLength());
			Array times = Array.factory(DataType.LONG, new int[] {NREC});

			int z;

			for (z = 0; z < lvlDim.getLength(); z++) {
				depths.set(z, mySpatialCoordinate[z]);
			}

			// These data are those created by bufferWriter class. This will write our hydraulic head (psi) and
			// adimensional water content (theta) data
			ArrayDouble.D2 dataWaterHeight = new ArrayDouble.D2(NREC, lvlDim.getLength());
			ArrayDouble.D2 dataPsi = new ArrayDouble.D2(NREC, lvlDim.getLength());
			ArrayDouble.D1 dataError =  new ArrayDouble.D1(NREC);
			
			int i=0;
			it = myVariables.entrySet().iterator();
			while (it.hasNext()) {
				
				Entry<String, ArrayList<double[]>> entry = (Entry<String, ArrayList<double[]>>) it.next();
				
				try {
					date = dateFormat.parse(entry.getKey());
				} catch (ParseException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				unixTime = (long) date.getTime()/1000;
				// think if there is a better way instead of using i
				times.setLong(i, unixTime);
				
						
				myTempVariable =  entry.getValue().get(0);
				for (int lvl = 0; lvl < NLVL; lvl++) {
					
					dataPsi.set(i, lvl, myTempVariable[lvl]);

				}
				
				myTempVariable =  entry.getValue().get(1);
				for (int lvl = 0; lvl < NLVL; lvl++) {
					
					dataWaterHeight.set(i, lvl, myTempVariable[lvl]);

				}

				dataError.set(i, entry.getValue().get(2)[0]);
				
				i++;
			}
		

			//Create the file. At this point the (empty) file will be written to disk
			dataFile.create();

			// A newly created Java integer array to be initialized to zeros.
			int[] origin = new int[2];

			dataFile.write(depthVar, depths);
			dataFile.write(timeVar, origin, times);
			dataFile.write(psiVar, origin, dataPsi);
			dataFile.write(waterHeightVar, origin, dataWaterHeight);
			dataFile.write(errorVar, origin, dataError);

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
