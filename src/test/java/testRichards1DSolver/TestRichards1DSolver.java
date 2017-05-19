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
import java.io.*;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;

import Richards1DSolver.*;

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
		String endDate = "1992-09-18 00:00";
		int timeStepMinutes = 60*24;
		String fId = "ID";


		String pathTopBC ="resources/Input/D_TopBoundaryConditionREAL.csv";
		String pathBottomBC ="resources/Input/D_BottomBoundaryCondition.csv";
		String pathIC = "resources/Input/InitialConditionTEST.txt";

		OmsTimeSeriesIteratorReader topBCReader = getTimeseriesReader(pathTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);

		double[] iC = ReadAndStoreIC(pathIC);
		double[] depth = ReadAndStoreDepth(pathIC);

		Richards1DSolver R1DSolver = new Richards1DSolver();

		R1DSolver.ks = 0.0000743;
		R1DSolver.thetaS =0.375 ;
		R1DSolver.thetaR = 0.058;
		R1DSolver.n = 3.177;
		R1DSolver.alpha =3.52 ;
		R1DSolver.lambda =1.9 ;
		R1DSolver.psiE = 1.9;
		R1DSolver.rMedian =1.9 ;
		R1DSolver.sigma =1.9 ;
		R1DSolver.soilHydraulicModel = "VanGenuchten";
		R1DSolver.topBCType = "Top Dirichlet";
		R1DSolver.bottomBCType = "Bottom Free Drainage";
		R1DSolver.spaceBottom = 2.0;
		R1DSolver.tTimestep = timeStepMinutes*60;
		R1DSolver.newtonTolerance = Math.pow(10,-3);
		R1DSolver.iC = iC;
		R1DSolver.depth = depth;
		R1DSolver.dir = "resources/Output";
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
	/* Per la lettura della condizione iniziale con la vecchia formattazione
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
	*/
	private double[] ReadAndStoreIC(String filePath) throws IOException{
		List<Double> l1 = new ArrayList<>();
		List<Double> l2 = new ArrayList<>();
		double[] readVector;
		
		Scanner s = new Scanner(new FileReader(filePath));

		while (s.hasNext()) {
		    l1.add(s.nextDouble());
		    l2.add(s.nextDouble());
		}

		s.close();
		
		readVector = new double[l1.size()];
		
		for(int i = 0;i < l1.size();i++){
			readVector[i] = l2.get(i);
		}
		return readVector;

	}
	private double[] ReadAndStoreDepth(String filePath) throws IOException{
		List<Double> l1 = new ArrayList<>();
		List<Double> l2 = new ArrayList<>();
		double[] readVector;
		
		Scanner s = new Scanner(new FileReader(filePath));

		while (s.hasNext()) {
		    l1.add(s.nextDouble());
		    l2.add(s.nextDouble());
		}

		s.close();
		
		readVector = new double[l1.size()];
		
		for(int i = 0;i < l1.size();i++){
			readVector[i] = l1.get(i);
		}
		return readVector;

	}
	
}