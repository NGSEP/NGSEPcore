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
package ngsep.genome;

/**
 * 
 * @author Jorge Gomez
 *
 */
public class HomologyEdge {
	private HomologyUnit queryUnit;
	private HomologyUnit subjectUnit;
	private double score;
	private double ks;
	private double ka;
	
	public HomologyEdge(HomologyUnit queryUnit, HomologyUnit subjectUnit, double score) {
		super();
		this.queryUnit = queryUnit;
		this.subjectUnit = subjectUnit;
		this.score = score;
	}

	public HomologyUnit getQueryUnit() {
		return queryUnit;
	}

	public HomologyUnit getSubjectUnit() {
		return subjectUnit;
	}

	public double getScore() {
		return score;
	}
	
	public void setScore(double score) {
		this.score = score;
	}

	public double getKs() {
		return ks;
	}

	public void setKs(double ks) {
		this.ks = ks;
	}

	public double getKa() {
		return ka;
	}

	public void setKa(double ka) {
		this.ka = ka;
	}

}
