package ngsep.sequences;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Scanner;

import ngsep.sequences.io.FastqFileReader;


public class RawReadTrimmer {
	
	private String FASTQ_FILE_EXTENSION = ".fastq.gz";
	private String ILUMINA_ADAPTER = "AGATCG";
	
	private String inputDirectory;
	private String outDirectory;
	private String fileNamesFile;
	private boolean pairedEnd = false;
	private ArrayList<String> fileNames = new ArrayList<>();
	
	public static void main(String[] args) throws Exception {
		RawReadTrimmer instance = new RawReadTrimmer();
		instance.inputDirectory = args[0];
		instance.fileNamesFile = args[1];
		instance.outDirectory = args[2];
		if(args[3] == "-p") {
			instance.pairedEnd = true;
		}
		
		instance.run();
	}

	private void run() throws FileNotFoundException, IOException {
		loadFileNames();
	}
	
	private void loadFileNames() throws FileNotFoundException, IOException {
		Scanner scanner = new Scanner(new File(fileNamesFile));
		while (scanner.hasNextLine()) {
		   String line = scanner.nextLine();
//		   line = line.strip();
//		   fileNames.add(line);
		}
	}
	
	private void walkFiles() throws IOException {
		if(pairedEnd) {
			for(String file: fileNames) {
				String file_1 = file + "_1" + FASTQ_FILE_EXTENSION;
				String file_2 = file + "_2" + FASTQ_FILE_EXTENSION;
				try(FastqFileReader reader_1 = new FastqFileReader(file_1);
						FastqFileReader reader_2 = new FastqFileReader(file_2)) {
					trimAdapters(reader_1, file_1, reader_2, file_2);
				}
			}
		} else {
			for(String file: fileNames) {
				try(FastqFileReader reader = new FastqFileReader(file + FASTQ_FILE_EXTENSION);) {
					trimAdapters(reader);
				}
			}
		}
		
	}

	private void trimAdapters(FastqFileReader reader) {
		// TODO Auto-generated method stub
		
	}

	private void trimAdapters(FastqFileReader reader_1, String file_1, FastqFileReader reader_2, String file_2) throws FileNotFoundException {
		try(PrintStream out_1 = new PrintStream(outDirectory + file_1);
				PrintStream out_2 = new PrintStream(outDirectory + file_2)) {
			Iterator<RawRead> it1 = reader_1.iterator();
			Iterator<RawRead> it2 = reader_2.iterator();
			while(it1.hasNext() && it2.hasNext()) {
				int cutIndex = 0;
				RawRead forward = it1.next();
				RawRead reverse = it2.next();
				cutIndex = forward.getSequenceString().indexOf(ILUMINA_ADAPTER);
				
			}
		}
 	}
	
}
