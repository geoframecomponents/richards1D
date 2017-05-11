package readingWritingSim;

// Modified http://www.mkyong.com/java8/java-8-stream-read-a-file-line-by-line/
// Completed by couple of others: 
// http://stackoverflow.com/questions/41699688/how-to-read-a-file-containing-text-and-integers-and-save-the-integers-in-an-arra
// For further comments see ReadAFile.java
// @Author Riccardo Rigon

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import Operazione.*;

public class ReadAndStoreInteger3 {
	
	public static void main(String[] args) {
		
		String fileName =  "C:/Users/Nico/Documents/DOC U N I V E R S I T A/Summer School OMS/oms3.prj.provaLetturaScrittura/data/Input.txt";
		List<Integer> list = new ArrayList<>();	

		try (Stream<String> stream = Files.lines(Paths.get(fileName))) {
			list = stream.
					filter(line -> !line.startsWith("!")).
					flatMap(line->Arrays.stream(line.split(" "))).
					map(Integer::valueOf).
					collect(Collectors.toList());
		} catch (IOException e) { 
			e.printStackTrace();
		}
		
		System.out.println(list.toString());

		//int lunghezza = list.size();
		//Integer[] vettore = new Integer[lunghezza];
		
		//vettore = (Integer[]) list.toArray();
		
		int[] vettore = new int[list.size()];
		for(int i = 0;i < vettore.length;i++){
			vettore[i] = list.get(i);
		  
		}
		
		somma2 s = new somma2();
		
		int[] risultato = s.sum(vettore);
		System.out.println("la somma è:\n");
		for(int i = 0;i < risultato.length;i++){
			System.out.println(+ risultato[i]);
		}
	}

}