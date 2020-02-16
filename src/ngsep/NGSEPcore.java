/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep;

import java.lang.reflect.Method;
import java.util.Arrays;


import ngsep.main.Command;
import ngsep.main.CommandsDescriptor;

public class NGSEPcore {

		
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		CommandsDescriptor descriptor = CommandsDescriptor.getInstance();
		if(args.length == 0 || args[0].equals("help") || args[0].equals("-h") || args[0].equals("--help")){
			descriptor.printUsage();
			return;
		} else if(args[0].equals("version") || args[0].equals("-v") || args[0].equals("--version")){
			descriptor.printVersion();
			return;
		} else if(args[0].equals("citing") || args[0].equals("-c") || args[0].equals("--citing")){
			descriptor.printCiting();
			return;
		}
		
		Command command = descriptor.getCommand(args[0]);
		if(command == null) {
			String currentCommandId = descriptor.getCurrentCommandId(args[0]);
			if(currentCommandId!=null) {
				System.err.println();
				System.err.println("ERROR: Command "+args[0]+" was replaced by command "+currentCommandId+ ". See new usage below.");
				System.err.println();
				command = descriptor.getCommand(currentCommandId);
				descriptor.printHelp(command);
			} else {
				System.err.println("ERROR: Unrecognized command "+args[0]);
			}
			System.exit(1);
		}
		Class<?> program = command.getProgram();
	    Class<?>[] argTypes = new Class[] { String[].class };
	    Method main = program.getDeclaredMethod("main", argTypes);
  	    String[] mainArgs = Arrays.copyOfRange(args, 1, args.length);
	    main.invoke(null, (Object)mainArgs);
	}

}
