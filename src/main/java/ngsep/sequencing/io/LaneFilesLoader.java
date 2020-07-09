package ngsep.sequencing.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.sequencing.SequencingLane;

public class LaneFilesLoader {

	public void loadFiles(String laneFilesDescriptor, List<SequencingLane> lanes) throws IOException {
		Map<String, SequencingLane> lanesMap = new HashMap<>();
		for(SequencingLane lane:lanes) {
			lanesMap.put(lane.getUniqueKey(), lane);
		}
		try (FileInputStream fis = new FileInputStream(laneFilesDescriptor);
			 BufferedReader in = new BufferedReader(new InputStreamReader(fis));) {
			String line = in.readLine();
			while (line != null) {
				String[] items = line.split("\t");
				SequencingLane lane = lanesMap.get((new SequencingLane(items[0], items[1]).getUniqueKey()));
				if(lane==null) continue;
				if(items.length<4) lane.addSingleEndFile(items[2]);
				else lane.addPairedEndFiles(items[2],items[3]);
				line = in.readLine();
			}
		}
	}
}
