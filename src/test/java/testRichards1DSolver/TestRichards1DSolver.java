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
import java.net.URISyntaxException;
import java.util.*;

import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;

import Richards1DSolver.*;

import org.junit.Test;
import richards_classes.ReadCsvTwoColumns;

/**
 * Test the {@link TestRichards1DSolver} module.
 * 
 * @author Niccolo' Tubini, Francesco Serafin
 */
public class TestRichards1DSolver {

	@Test
	public void Test() throws Exception {


		String startDate = "2017-01-01 18:00" ;
		String endDate = "2017-01-02 23:55";
		int timeStepMinutes = 5;
		String fId = "ID";


		String pathTopBC ="resources/Input/Trento0.csv";
		String pathBottomBC ="resources/Input/TrentoBottom - Copy.csv";
		String pathIC = "resources/Input/InitialConditionHydrostatic.csv";
		String pathSourceSink = "resources/Input/SourceSink0.csv";

		String pathTempTopBC ="resources/Input/Temperature_2Gradini.csv";
		String pathTempBottomBC ="resources/Input/Temperature_TrentoBottom.csv";
		String pathTempIC = "resources/Input/InitialConditionTemperature.csv";
	
		OmsTimeSeriesIteratorReader topBCReader = getTimeseriesReader(pathTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);

		OmsTimeSeriesIteratorReader tempTopBCReader = getTimeseriesReader(pathTempTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader tempBottomBCReader = getTimeseriesReader(pathTempBottomBC, fId, startDate, endDate, timeStepMinutes);
		
		ReadCsvTwoColumns readIC = new ReadCsvTwoColumns();
		readIC.setFilePath(pathIC);
		readIC.process();
		double[] iC = readIC.getVariable();
		double[] depth = readIC.getDepth();
		
		ReadCsvTwoColumns readSourceSink = new ReadCsvTwoColumns();
		readSourceSink.setFilePath(pathSourceSink);
		readSourceSink.process();
		double[] sourceSink = readSourceSink.getVariable();

		ReadCsvTwoColumns readTempIC = new ReadCsvTwoColumns();
		readIC.setFilePath(pathTempIC);
		readIC.process();
		double[] temperatureIC = readIC.getVariable();
		double[] temperatureDepth = readIC.getDepth();
		
		Richards1DSolver R1DSolver = new Richards1DSolver();

		R1DSolver.ks = 0.00000007;
		R1DSolver.thetaS =0.35;
		R1DSolver.thetaR = 0.02;
		R1DSolver.n = 1.964;
		R1DSolver.alpha = 4.1;
		R1DSolver.lambda =1.5 ;
		R1DSolver.psiE = -1/(0.0286/0.01);
		R1DSolver.rMedian =0.0000020781 ;
		R1DSolver.sigma =0.6 ;
		R1DSolver.soilHydraulicModel = "VanGenuchten";
		R1DSolver.topBCType = "Psi Top Neumann";
		R1DSolver.bottomBCType = "Psi Bottom Dirichlet";
		R1DSolver.tempTopBCType = "Temperature Diffusion Top Dirichlet";
		R1DSolver.tempBottomBCType = "Temperature Diffusion Bottom Dirichlet";
		R1DSolver.tempConvectionTopBCType = "Temperature Convection Top Dirichlet";
		R1DSolver.tempConvectionBottomBCType = "Temperature Convection Bottom Dirichlet";
		R1DSolver.energySolver = "diffusion convection";
		R1DSolver.soilType = "mineral soil"; // "unfrozen peat" "mineral soil" "crushed rock"
		R1DSolver.sandFraction = 0.5;
		R1DSolver.clayFraction = 0.2;
		R1DSolver.lambda0 = 2.5; // 2<lambda0<3
		R1DSolver.quartzFraction = 0.1;
		R1DSolver.soilThermalModel = "Johansen";
		R1DSolver.kappaWithTemperatureModel = "Isothermal flow";
		R1DSolver.delta = 0;
		R1DSolver.spaceBottom = 2.0;
		R1DSolver.tTimestep = 300;
		R1DSolver.timeDelta =300;
		R1DSolver.newtonTolerance = Math.pow(10,-14);
		R1DSolver.iC = iC;
		R1DSolver.depth = depth;
		R1DSolver.sourceSink = sourceSink;
		R1DSolver.temperatureIC = temperatureIC;
		R1DSolver.temperatureDepth = temperatureDepth;
		R1DSolver.dir = "resources/Output clay 2";
		R1DSolver.nestedNewton =1;
		

		while( topBCReader.doProcess  ) {

			topBCReader.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = topBCReader.outData;
			R1DSolver.inTopBC= bCValueMap;


			bottomBCReader.nextRecord();
			bCValueMap = bottomBCReader.outData;
			R1DSolver.inBottomBC = bCValueMap;
			
			tempTopBCReader.nextRecord();	
			bCValueMap = tempTopBCReader.outData;
			R1DSolver.inTemperatureTopBC= bCValueMap;
			
			tempBottomBCReader.nextRecord();
			bCValueMap = tempBottomBCReader.outData;
			R1DSolver.inTemperatureBottomBC = bCValueMap;

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
}