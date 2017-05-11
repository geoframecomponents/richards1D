
package readingWritingSim;

import oms3.annotations.*;
import oms3.io.CSTable;
import oms3.io.DataIO;

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


public class ReadAndStoreInteger2 {
		

		public String filePath = "C:/Users/Nico/Documents/DOC U N I V E R S I T A/Summer School OMS/oms3.prj.provaLetturaScrittura/data/Input.txt";
		

		public double value;
		
		
		public void readTxt() throws IOException {
		   
		    File file = new File(filePath);
		    Scanner scanner = new Scanner(file);
		
		    while (scanner.hasNext()) {
		        value = scanner.nextDouble();
		        System.out.println("the value is: "+value);
		    	//System.out.println(scanner.next());
		    }
		    
		}

}
