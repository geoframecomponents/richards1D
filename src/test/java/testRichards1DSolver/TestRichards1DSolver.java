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

package testRichards1DSolver;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;

import Richards1DSolver.Richards1DSolver;

import org.junit.Test;

/**
 * Test the {@link TestRichards1DSolver} module.
 * 
 * @author Niccolò Tubini 
 */
public class TestRichards1DSolver {

	@Test
	public void Test() throws Exception {


		String startDate = "1991-09-18 00:00" ;
		String endDate = "1992-07-14 00:00";
		int timeStepMinutes = 60*24;
		String fId = "ID";


		String pathTopBC ="resources/Input/TopBoundaryCondition.csv";
		String pathBottomBC ="resources/Input/BottomBoundaryCondition.csv";
		String pathIC = "resources/Input/InitialCondition.txt";

		OmsTimeSeriesIteratorReader topBCReader = getTimeseriesReader(pathTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);

		double[] iC = ReadAndStoreDouble(pathIC);

		Richards1DSolver R1DSolver = new Richards1DSolver();

		R1DSolver.ks = 0.0000007 ;
		R1DSolver.thetaS =0.41 ;
		R1DSolver.thetaR = 0.095;
		R1DSolver.n = 1.31;
		R1DSolver.alpha =1.9 ;
		R1DSolver.lambda =1.9 ;
		R1DSolver.psiE = 1.9;
		R1DSolver.rMedian =1.9 ;
		R1DSolver.sigma =1.9 ;
		R1DSolver.soilHydraulicModel = "VanGenuchten";
		R1DSolver.spaceBottom = 2.0;
		R1DSolver.spaceTop = 0.0;
		R1DSolver.tTimestep = 1000; // !! timeStepMinutes*60
		R1DSolver.newtonTolerance = Math.pow(10,-12);
		R1DSolver.iC = iC;
		R1DSolver.dir = "resources/Output/";
		while( topBCReader.doProcess  ) {

			topBCReader.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = topBCReader.outData;
			R1DSolver.inTopBC= bCValueMap;


			bottomBCReader.nextRecord();
			bCValueMap = bottomBCReader.outData;
			R1DSolver.inBottomBC = bCValueMap;

			R1DSolver.inCurrentDate = topBCReader.tCurrent;
			
			R1DSolver.solve();		

		}

		topBCReader.close();
		bottomBCReader.close();
	}

	private OmsTimeSeriesIteratorReader getTimeseriesReader( String inPath, String id, String startDate, String endDate,
			int timeStepMinutes ) throws URISyntaxException {
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = inPath;
		reader.idfield = "ID";
		reader.tStart = startDate;
		reader.tTimestep = timeStepMinutes;
		reader.tEnd = endDate;
		reader.fileNovalue = "-9999";
		reader.initProcess();
		return reader;
	}

	private double[] ReadAndStoreDouble(String filePath) {

		double[] readVector;
		List<Double> list = new ArrayList<>();	

		try (Stream<String> stream = Files.lines(Paths.get(filePath))) {
			list = stream.
					filter(line -> !line.startsWith("!")).
					flatMap(line->Arrays.stream(line.split(" "))).
					map(Double::valueOf).
					collect(Collectors.toList());
		} catch (IOException e) { 
			e.printStackTrace();
		}

		System.out.println("Reading completed");


		readVector = new double[list.size()];
		for(int i = 0;i < readVector.length;i++){
			readVector[i] = list.get(i);
		}

		return readVector;

	}
}