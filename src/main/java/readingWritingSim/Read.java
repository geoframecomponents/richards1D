package readingWritingSim;

import oms3.annotations.*;
import oms3.io.CSTable;
import oms3.io.DataIO;
//import richards_utils.TextIO;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Read {
	//@Description("Path of the input file")
	//@In
	public String filePath = "C:/Users/Nico/Documents/DOC U N I V E R S I T A/Summer School OMS/oms3.prj.provaLetturaScrittura/data/Input.txt";
	
	//@Description("Input value read")
	//@Out
	public double value;
	
	//@Execute
	//public void readTxt() throws IOException {
		try{
			TextIO.readFile(filePath);
		}
		catch (IllegalArgumentException e) {
			System.out.println("Can't open file " + filePath + " for reading!");
			System.exit(1);  
		}	
	    
	}
}
