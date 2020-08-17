package ngsep.sequencing;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.Scanner;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import ngsep.sequences.DegenerateSequence;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;


public class RawReadsTrimmer {
	
	private String adapter = "AGATCG";
	
	private String inputDirectory;
	private String outDirectory;
	private String sampleNamesFile;
	
	public static void main(String[] args) throws Exception {
		RawReadsTrimmer instance = new RawReadsTrimmer();
		instance.inputDirectory = args[0];
		instance.sampleNamesFile = args[1];
		instance.outDirectory = args[2];
		instance.adapter = args[3];
		instance.run();
	}

	private void run() throws FileNotFoundException, IOException {
		loadFileNames();
	}
	
	private void loadFileNames() throws IOException {
		try (Scanner scanner = new Scanner(new File(sampleNamesFile))) {
			while (scanner.hasNextLine()) {
				String line = scanner.nextLine();
				String [] items = line.split("\t");
				try(FastqFileReader reader1 = new FastqFileReader(inputDirectory+File.separator+items[1]);
					FastqFileReader reader2 = new FastqFileReader(inputDirectory+File.separator+items[2])) {
					trimAdapters(reader1, items[1], reader2, items[2]);
				}
			}
		}
	}

	private void trimAdapters(FastqFileReader reader1, String file1, FastqFileReader reader2, String file2) throws IOException {
		Pattern pattern = Pattern.compile(DegenerateSequence.makeRegularExpression(adapter));
		try(OutputStream os1 = new GZIPOutputStream(new FileOutputStream(outDirectory + File.separator+file1));
			PrintStream out1 = new PrintStream(os1);
			OutputStream os2 = new GZIPOutputStream(new FileOutputStream(outDirectory + File.separator+file2));
			PrintStream out2 = new PrintStream(os2)) {
			Iterator<RawRead> it1 = reader1.iterator();
			Iterator<RawRead> it2 = reader2.iterator();
			while(it1.hasNext() && it2.hasNext()) {
				RawRead read1 = it1.next();
				RawRead read2 = it2.next();
				int l1 = read1.getLength();
				read1.trimFromSequence(pattern);

				if(read1.getLength()!=l1) {
					l1 = read1.getLength();
					read2.trimToLength(l1);
				}
				if(l1>=40 ) {
					read1.save(out1);
					read2.save(out2);
				}
				
			}
		}
 	}
	
}
