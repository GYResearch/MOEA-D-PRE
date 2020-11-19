package jmetal.metaheuristics.MOEADPRE;

import java.util.ArrayList;
import java.util.List;

import jmetal.core.Problem;

public class InitWeight {
	public int numberofobjectives;
	public int populationsize;
	
	public InitWeight(int numberofobjectives,int populationsize) {
		super();
		this.numberofobjectives = numberofobjectives;
		this.populationsize =populationsize;
	}

	/**
	 * (m,H1,H2,size) (2,0,0,100/101),
	 * (3,13,0,105);(6,4,1,132);(8,3,2,156);(10,3,2,275); generate weight
	 * vectors
	 * 
	 * @return
	 */
	public double[][] createWeight() {
		int H1 = 0, H2 = 0;
		if (numberofobjectives == 2) {// popsize=101; or 100
			double[][] list = new double[populationsize][numberofobjectives];
			for (int i = 0; i < populationsize; i++) {
				list[i][0] = 1.0 * i / (populationsize - 1);;
				list[i][1] = 1 - list[i][0];
			}
			return list;
		} else if (numberofobjectives == 3) {
			H1 = 13;
			H2 = 0;
			NormalBoundaryIntersectionGenerateWeight nn = new NormalBoundaryIntersectionGenerateWeight(
					numberofobjectives, H1, H2);
			return exchange(nn.generate());
		} else if (numberofobjectives == 6) {
			H1 = 4;
			H2 = 1;
			NormalBoundaryIntersectionGenerateWeight nn = new NormalBoundaryIntersectionGenerateWeight(
					numberofobjectives, H1, H2);
			return exchange(nn.generate());
		} else if (numberofobjectives == 8) {
			H1 = 3;
			H2 = 2;
			NormalBoundaryIntersectionGenerateWeight nn = new NormalBoundaryIntersectionGenerateWeight(
					numberofobjectives, H1, H2);
			return exchange(nn.generate());
		} else if (numberofobjectives == 10) {
			H1 = 3;
			H2 = 2;
			NormalBoundaryIntersectionGenerateWeight nn = new NormalBoundaryIntersectionGenerateWeight(
					numberofobjectives, H1, H2);
			return exchange(nn.generate());
		} else {
			System.out.println("You need to update the InitWeight!!!!");
			return null;
		}

	}
	
	public double[][] exchange( List<double[]> list) {
		double[][] temp = new double[populationsize][numberofobjectives];
		for (int i = 0; i < list.size(); i++) {
			for (int j = 0; j < numberofobjectives; j++) {
				temp[i][j] = list.get(i)[j];
			}
		}
		return temp;
	}
	

}
