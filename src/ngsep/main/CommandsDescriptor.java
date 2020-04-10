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
package ngsep.main;

import java.io.InputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

/**
 * Class to handle the different commands available in NGSEP
 * @author Jorge Duitama
 *
 */
public class CommandsDescriptor {
	public static final String ATTRIBUTE_VERSION="version";
	public static final String ATTRIBUTE_DATE="date";
	public static final String ATTRIBUTE_ID="id";
	public static final String ATTRIBUTE_CLASSNAME="class";
	public static final String ATTRIBUTE_TYPE="type";
	public static final String ATTRIBUTE_DEFAULT_VALUE="default";
	public static final String ATTRIBUTE_DEFAULT_CONSTANT="defaultConstant";
	public static final String ATTRIBUTE_ATTRIBUTE="attribute";
	public static final String ATTRIBUTE_GROUPID="groupId";
	public static final String ATTRIBUTE_MULTIPLE="multiple";
	public static final String ATTRIBUTE_FORMERID="formerId";
	public static final String ATTRIBUTE_DEPRECATED="deprecated";
	public static final String ATTRIBUTE_PRINTHELP="printHelp";
	public static final String ELEMENT_COMMAND="command";
	public static final String ELEMENT_COMMANDGROUP="commandgroup";
	public static final String ELEMENT_TITLE="title";
	public static final String ELEMENT_INTRO="intro";
	public static final String ELEMENT_DESCRIPTION="description";
	public static final String ELEMENT_USAGE="usage";
	public static final String ELEMENT_ARGUMENT="argument";
	public static final String ELEMENT_OPTION="option";
	
	private String resource = "/ngsep/main/CommandsDescriptor.xml";
	private String swVersion;
	private String releaseDate;
	private String swTitle;
	private Map<String,List<Command>> commandsByGroup = new HashMap<String,List<Command>>();
	private Map<String,Command> commandsByClass = new HashMap<String,Command>();
	private Map<String,Command> commandsById = new HashMap<String,Command>();
	private Map<String, String> commandGroupNames = new LinkedHashMap<String, String>();
	private Map<String,String> formerCommandIds = new HashMap<String, String>();
	public static CommandsDescriptor instance = new CommandsDescriptor();
	/**
	 * Private constructor to implement the singleton pattern
	 */
	private CommandsDescriptor () {
		load();
	}
	public static CommandsDescriptor getInstance() {
		return instance;
	}
	/**
	 * Loads the commands descriptor XML
	 */
	private void load() {
		
		Document doc;
		try (InputStream is = this.getClass().getResourceAsStream(resource)) { 
			DocumentBuilder documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
			doc = documentBuilder.parse(new InputSource(is));
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		Element rootElement = doc.getDocumentElement();
		swVersion = rootElement.getAttribute(ATTRIBUTE_VERSION);
		releaseDate = rootElement.getAttribute(ATTRIBUTE_DATE);
		loadSoftwareDescription(rootElement);
		
	}
	/**
	 * Loads the description of the software
	 * @param parent
	 */
	private void loadSoftwareDescription(Element parent) {
		NodeList offspring = parent.getChildNodes();
		for(int i=0;i<offspring.getLength();i++){  
			Node node = offspring.item(i);
			if (node instanceof Element){ 
				Element elem = (Element)node;
				if(ELEMENT_COMMANDGROUP.equals(elem.getNodeName())) {
					loadCommandGroup(elem);
				}
				if(ELEMENT_COMMAND.equals(elem.getNodeName())) {
					Command c;
					try {
						c = loadCommand(elem);
					} catch (RuntimeException e) {
						throw new RuntimeException("Can not load command with id "+elem.getAttribute(ATTRIBUTE_ID),e);
					}
					List<Command> commandsList = commandsByGroup.get(c.getGroupId());
					if(commandsList==null)  throw new RuntimeException("Group "+c.getGroupId()+" not registered for command with id: "+c.getId());
					if(commandsById.containsKey(c.getId())) throw new RuntimeException("Duplicated command id: "+c.getId());
					commandsById.put(c.getId(),c);
					commandsList.add(c);
					commandsByClass.put(c.getProgram().getName(), c);
				} else if(ELEMENT_TITLE.equals(elem.getNodeName())) {
					swTitle = loadText(elem);
				} 
			}
		}
		
	}
	/**
	 * Loads the basic information of a command group specified by the given element
	 * @param elem XML element with the command group description
	 */
	private void loadCommandGroup(Element elem) {
		String id = elem.getAttribute(ATTRIBUTE_ID);
		if(id==null) throw new RuntimeException("Every command group must have an id");
		String name = loadText(elem);
		commandGroupNames.put(id,name);
		commandsByGroup.put(id,new ArrayList<Command>());
	}
	/**
	 * Loads the information of a specific command
	 * @param cmdElem XML element with the command description
	 * @return Command
	 */
	private Command loadCommand(Element cmdElem) {
		String id = cmdElem.getAttribute(ATTRIBUTE_ID);
		if(id==null) throw new RuntimeException("Every command must have an id");
		String className = cmdElem.getAttribute(ATTRIBUTE_CLASSNAME);
		if(className==null) throw new RuntimeException("Every command must have a class name");
		Class<?> program;
		try {
			program = Class.forName(className);
			Class<?>[] argTypes = new Class[] { String[].class };
			program.getDeclaredMethod("main",argTypes);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException("Can not load class for command: "+id,e);
		} catch (NoSuchMethodException e) {
			throw new RuntimeException("Class "+className+" for command: "+id+" does not have a main method",e);
		} catch (SecurityException e) {
			throw new RuntimeException("Class "+className+" for command: "+id+" can not be called",e);
		}
		Command cmd = new Command(id,program);
		String printHelpStr = cmdElem.getAttribute(ATTRIBUTE_PRINTHELP);
		cmd.setPrintHelp(!"false".equals(printHelpStr));
		String groupId = cmdElem.getAttribute(ATTRIBUTE_GROUPID);
		cmd.setGroupId(groupId);
		String formerId = cmdElem.getAttribute(ATTRIBUTE_FORMERID);
		if(formerId!=null) {
			formerCommandIds.put(formerId, cmd.getId());
		}
		NodeList offspring = cmdElem.getChildNodes(); 
		for(int i=0;i<offspring.getLength();i++){  
			Node node = offspring.item(i);
			if (node instanceof Element){ 
				Element elem = (Element)node;
				if(ELEMENT_TITLE.equals(elem.getNodeName())) {
					cmd.setTitle(loadText(elem));
				} else if(ELEMENT_INTRO.equals(elem.getNodeName())) {
					cmd.setIntro(loadText(elem));
				} else if(ELEMENT_DESCRIPTION.equals(elem.getNodeName())) {
					cmd.setDescription(loadText(elem));
				} else if(ELEMENT_ARGUMENT.equals(elem.getNodeName())) {
					String optMultiple = elem.getAttribute(ATTRIBUTE_MULTIPLE);
					boolean multiple = false;
					if(optMultiple!=null) multiple = "true".equals(optMultiple.trim().toLowerCase());
					cmd.addArgument(loadText(elem),multiple);
				} else if(ELEMENT_OPTION.equals(elem.getNodeName())) {
					String optId = elem.getAttribute(ATTRIBUTE_ID);
					if(optId==null || optId.length()==0) throw new RuntimeException("Every option must have an id");
					CommandOption opt = new CommandOption(optId);
					String optType = elem.getAttribute(ATTRIBUTE_TYPE);
					if(optType!=null && optType.length()>0) opt.setType(optType);
					String optDefault = elem.getAttribute(ATTRIBUTE_DEFAULT_VALUE);
					if(optDefault!=null && optDefault.trim().length()>0) opt.setDefaultValue(optDefault);
					String optDefaultConstant = elem.getAttribute(ATTRIBUTE_DEFAULT_CONSTANT);
					if(optDefaultConstant!=null && optDefaultConstant.trim().length()>0) opt.setDefaultValue(loadValue(program,optDefaultConstant));
					String optAttribute = elem.getAttribute(ATTRIBUTE_ATTRIBUTE);
					if(optAttribute!=null && optAttribute.trim().length()>0) opt.setAttribute(optAttribute);
					String optDeprecated = elem.getAttribute(ATTRIBUTE_DEPRECATED);
					if(optDeprecated!=null) opt.setDeprecated("true".equals(optDeprecated.trim().toLowerCase()));
					String description = loadText(elem);
					if(description==null || description.length()==0) throw new RuntimeException("Option "+optId+" does not have a description");
					opt.setDescription(description);
					cmd.addOption(opt);
				}	
			}
		}
		return cmd;
	}
	private String loadValue(Class<?> program, String constantName) {
		try {
			return ""+program.getDeclaredField(constantName).get(null);
		} catch (IllegalArgumentException | IllegalAccessException | NoSuchFieldException | SecurityException e) {
			throw new RuntimeException(e);
		}
		
	}
	/**
	 * Loads a text node as a String
	 * @param elem Text element
	 * @return String loaded text
	 */
	private String loadText(Element elem) {
		NodeList offspring = elem.getChildNodes();
		for (int i=0; i < offspring.getLength(); i++) {
	        Node subnode = offspring.item(i);
	        if (subnode.getNodeType() == Node.TEXT_NODE) {
	            String desc = subnode.getNodeValue();
	            if(desc!=null) {
	            	desc = desc.trim();
	            	//return desc;
	            	return desc.replaceAll("\\s", " ");
	            }
	        }
		}
		return null;
	}
	public String getSwVersion() {
		return swVersion;
	}
	public String getReleaseDate() {
		return releaseDate;
	}
	public Command getCommand(String name) {
		return commandsById.get(name);
	}
	public Command getCommandByClass(String classname) {
		return commandsByClass.get(classname);
	}
	/**
	 * Prints the general usage including the command names and intro information
	 */
	public void printUsage(){
		System.err.println();
		printVersionHeader();
		System.err.println("=============================================================================");
		System.err.println();
		System.err.println("USAGE: java -jar NGSEPcore_"+swVersion+".jar <COMMAND> <OPTIONS> <ARGUMENTS>");
		System.err.println();
		for(String commandGroupId:commandGroupNames.keySet()) {
			System.err.println("Commands for "+commandGroupNames.get(commandGroupId));
			List<Command> commandsGroup = commandsByGroup.get(commandGroupId);
			System.err.println();
			for(Command c:commandsGroup) {
				if(c.isPrintHelp()) {
					System.err.println("  > " + c.getId());
					System.err.println("          "+c.getIntro());
				}
			}
			System.err.println();
		}
		
		
		System.err.println();
		System.err.println("See http://sourceforge.net/projects/ngsep/files/Library/ for more details.");
		System.err.println();		
	}
	/**
	 * Prints the software version
	 */
	private void printVersionHeader() {
		System.err.println(" NGSEP - "+swTitle);
		System.err.println(" Version " + swVersion + " ("+releaseDate+")");
	}
	/**
	 * Prints the help for a specific program. Locates the command from the Class where the program is implemented
	 * @param program
	 */
	public void printHelp(Class<?> program) {
		Command c = commandsByClass.get(program.getName());
		printHelp(c);
	}
	/**
	 * Prints the help for a specific program. Locates the command from the Class where the program is implemented
	 * @param program
	 */
	public void printHelp(Command c) {
		
		int titleLength = c.getTitle().length();
		for(int i=0;i<titleLength;i++)System.err.print("-");
		System.err.println();
		System.err.println(c.getTitle());
		for(int i=0;i<titleLength;i++)System.err.print("-");
		System.err.println();
		System.err.println();
		System.err.println(c.getDescription());
		System.err.println();
		System.err.println("USAGE:");
		System.err.println();
		System.err.print("java -jar NGSEPcore_"+swVersion+".jar "+ c.getId()+" <OPTIONS>");
		for(String arg:c.getArguments()) {
			
			System.err.print(" <"+arg+">");
			if(c.isMultiple(arg)) System.err.print("*");
		}
		System.err.println();
		System.err.println();
		System.err.println("OPTIONS:");
		System.err.println();
		List<CommandOption> options = c.getOptionsList(); 
		int longestOpt = getLongestOption(options);
		for(CommandOption option:options) {
			if(option.isDeprecated()) continue;
			System.err.print("        -"+option.getId());
			if(option.printType() ) System.err.print(" "+option.getType());
			int optLength = option.getPrintLength();
			int diff = longestOpt-optLength;
			for(int i=0;i<diff+1;i++)System.err.print(" ");
			System.err.print(": ");
			String desc = option.getDescription();
			if(option.getDefaultValue()!=null) desc+=" Default: "+option.getDefaultValue();
			printDescription(desc,longestOpt+3);
		}
		System.err.println();
	}
	
	private int getLongestOption(List<CommandOption> options) {
		int max = 0;
		for(CommandOption opt:options) {
			if(opt.isDeprecated()) continue;
			int l = opt.getPrintLength();
			if(max<l) max = l;
		}
		return max;
	}
	/**
	 * Prints a general description
	 * @param desc Description to be printed in standard error
	 * @param startColumn for formatting
	 */
	private void printDescription(String desc, int startColumn) {
		//TODO: Print in a command line friendly format
		System.err.println(desc);
		
	}
	/**
	 * Prints the version plus general usage
	 */
	public void printVersion() {
		System.err.println();
		printVersionHeader();
		System.err.println();
		System.err.println(" For usage type     java -jar NGSEPcore_"+swVersion+".jar --help");
		System.err.println();
		System.err.println(" For citing type     java -jar NGSEPcore_"+swVersion+".jar --citing");
		System.err.println();
	}
	/**
	 * Prints the citing information
	 */
	public void printCiting(){
		System.err.println("------");
		System.err.println("Citing");
		System.err.println("------");
		System.err.println();
		System.err.println("To cite NGSEP, please include in your references the following manuscript:");
		System.err.println();
		System.err.println("Tello D, Gil J, Loaiza CD, Riascos JJ, Cardozo N, and Duitama J. (2019)");
		System.err.println("NGSEP3: accurate variant calling across species and sequencing protocols.");
		System.err.println("Bioinformatics 35(22): 4716-4723. http://doi.org/10.1093/bioinformatics/btz275");
		System.err.println();
		System.err.println("See the README.txt file for papers describing the algorithms implemented in NGSEP and supporting packages");
		System.err.println();
	}
	
	/**
	 * Loads the optional fields of a program
	 * @param programInstance Object of a program implementing one command
	 * @param args Arguments sent by the user
	 * @return int Next index to be processed in the arguments array
	 */
	public int loadOptions(Object programInstance, String [] args ) throws Exception {
		if (args.length == 0 || (args.length==1 && args[0].equals("-h") || args[0].equals("--help"))){
			printHelp(programInstance.getClass());
			System.exit(1);
		}
		Command c = commandsByClass.get(programInstance.getClass().getName());
		int i = 0;
		while(i<args.length && args[i].charAt(0)=='-') {
			if("-".equals(args[i])) break;
			CommandOption o = c.getOption(args[i].substring(1));
			if (o==null) {
				System.err.println("Unrecognized option "+args[i]);
				printHelp(programInstance.getClass());
				System.exit(1);
			} else if (o.isDeprecated()) {
				System.err.println("WARN: Deprecated option "+args[i]);
				System.err.println(o.getDescription());
				i++;
				if(!CommandOption.TYPE_BOOLEAN.equals(o.getType())) i++;
				continue;
			}
			Method setter = o.findSetMethod(programInstance);
			Object value=null;
			if(CommandOption.TYPE_BOOLEAN.equals(o.getType())) {
				value = true;
			} else {
				i++;
				if(setter.getParameterTypes()[0].equals(String.class)) {
					value = args[i];
				} else {
					try {
						value = o.decodeValue(args[i]);
					} catch (NumberFormatException e) {
						System.err.println("Error loading value \""+args[i]+"\" for option \""+o.getAttribute()+"\" of type "+o.getType()+": "+e.getMessage());
						printHelp(programInstance.getClass());
						System.exit(1);
					}
				}
			}
			try {
				setter.invoke(programInstance, value);
			} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
				System.err.println("Error setting value \""+value+"\" for option \""+o.getId()+"\" of type: "+o.getType()+". "+e.getMessage());
				e.printStackTrace();
				System.exit(1);
			}
			i++;
		}
		return i;
	}
	public String getCurrentCommandId(String formerId) {
		return formerCommandIds.get(formerId);
	}
	
}
