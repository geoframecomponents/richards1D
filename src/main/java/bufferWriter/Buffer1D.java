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

package bufferWriter;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import oms3.annotations.*;

@Description("Buffer for 1D simulation in which there is one spatial coordinate and time.")
@Documentation("")
@Author(name = "Niccolo' Tubini, Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Hydrology, Richards, Infiltration")
//@Label(JGTConstants.HYDROGEOMORPHOLOGY)
//@Name("shortradbal")
//@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")

public class Buffer1D {
	
	@Description("Varible to store")
	@In 
	@Unit ("-")
	public ArrayList<double[]> inputVariable;
	
	@Description("Date at which the varible is computed")
	@In 
	@Unit ("YYYY-MM-DD HH:mm")
	public String inputDate;
	
	@Description("Spatial coordinate: is the position of the centroids ")
	@In 
	@Unit ("m")
	public double[] inputSpatialCoordinate;
	
	
	@Description()
	@Out
	@Unit ()
	public LinkedHashMap<String,ArrayList<double[]>> myVariable; // consider the opportunity to save varibale as float instead of double
	
	@Description()
	@Out
	@Unit ()
	public double[] mySpatialCoordinate;
	
	@Description("")
	int step=0;
	
	ArrayList<double[]> tempVariable;
	
	public Buffer1D() {
		
		myVariable = new LinkedHashMap<String,ArrayList<double[]>>();
		
		
	}
	
	@Execute
	public void solve() {
		//System.out.println("Buffer1D step:" + step);
		if(step==0){
			
		mySpatialCoordinate = inputSpatialCoordinate;
		tempVariable = new ArrayList<double[]>();
		//System.out.println(mySpatialCoordinate.toString());
		
		}
		
		
		tempVariable.add(inputVariable.get(0).clone());
		tempVariable.add(inputVariable.get(1).clone());
		tempVariable.add(inputVariable.get(2).clone());
		tempVariable.add(inputVariable.get(3).clone());
		tempVariable.add(inputVariable.get(4).clone());
		tempVariable.add(inputVariable.get(5).clone());
		
		myVariable.put(inputDate,(ArrayList<double[]>) tempVariable.clone());
		//System.out.println(myVariable.size() +"       "+ myVariable.keySet());
		//System.out.println(myVariable.toString());
		tempVariable.clear();
		step++;
		
	}
	

}
