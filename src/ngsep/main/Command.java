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

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class Command {
	private String id;
	private Class<?> program;
	private String title;
	private String intro;
	private String description;
	private String groupId;
	private boolean printHelp = true;
	
	private Map<String,Boolean> arguments = new LinkedHashMap<String,Boolean>();
	private Map<String, CommandOption> options = new LinkedHashMap<String,CommandOption>();
	public Command(String id, Class<?> program) {
		super();
		this.id = id;
		this.program = program;
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public Class<?> getProgram() {
		return program;
	}
	public void setProgram(Class<?> program) {
		this.program = program;
	}
	public String getTitle() {
		return title;
	}
	public void setTitle(String title) {
		this.title = title;
	}
	public String getIntro() {
		return intro;
	}
	public void setIntro(String intro) {
		this.intro = intro;
	}
	public String getDescription() {
		return description;
	}
	public void setDescription(String description) {
		this.description = description;
	}
	public String getGroupId() {
		return groupId;
	}
	public void setGroupId(String groupId) {
		this.groupId = groupId;
	}
	public boolean isPrintHelp() {
		return printHelp;
	}
	public void setPrintHelp(boolean printHelp) {
		this.printHelp = printHelp;
	}
	public List<String> getArguments() {
		List<String> answer = new ArrayList<String>();
		answer.addAll(arguments.keySet());
		return answer;
	}
	public boolean isMultiple(String argument) {
		Boolean mult = arguments.get(argument);
		if(mult==null) mult = false;
		return mult.booleanValue();
	}
	public Map<String, CommandOption> getOptions() {
		return options;
	}
	public List<CommandOption> getOptionsList() {
		List<CommandOption> optionsList = new ArrayList<CommandOption>();
		for(String optId:options.keySet()) {
			optionsList.add(options.get(optId));
		}
		return optionsList;
	}
	public void addArgument(String argument, boolean isMultiple) {
		arguments.put(argument,isMultiple);
	}
	public void addOption(CommandOption option) {
		if(options.containsKey(option.getId())) throw new IllegalArgumentException("Duplicated option id: "+option.getId());
		options.put(option.getId(), option);
	}
	public CommandOption getOption(String optionId) {
		return options.get(optionId);
	}
}
