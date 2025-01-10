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
package ngsep.hmm;

/**
 * @author Leonidas Villamil
 * @author Jorge Duitama
 */
public class ProfileAlignmentDomain {
	private final String query;
	private int start;
	private int end;
	private String hmmID;
	private double evalue;
	private String alignmentQuery;
	private String alignmentProfile;
	public ProfileAlignmentDomain(String query, int start, int end, String hmmID) {
		super();
		this.query = query;
		this.start = start;
		this.end = end;
		this.hmmID = hmmID;
	}
	public double getEvalue() {
		return evalue;
	}
	public void setEvalue(double evalue) {
		this.evalue = evalue;
	}
	public String getQuery() {
		return query;
	}
	public int getStart() {
		return start;
	}
	public int getEnd() {
		return end;
	}
	public String getHmmID() {
		return hmmID;
	}
	
	public String getAlignmentQuery() {
		return alignmentQuery;
	}
	public String getAlignmentProfile() {
		return alignmentProfile;
	}
	public void setAlignment(String alignmentQuery, String alignmentProfile) {
		this.alignmentQuery = alignmentQuery;
		this.alignmentProfile = alignmentProfile;
	}
}
