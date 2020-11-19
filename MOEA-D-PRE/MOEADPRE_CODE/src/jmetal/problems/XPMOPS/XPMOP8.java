
package jmetal.problems.XPMOPS;

import jmetal.core.*;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;
import jmetal.util.Configuration.*;

/**
 * Class representing problem PMOP8: g_8, l_2, k_3, h_2
 * 
    * You need to set the parameters 
 * int A = 6, B = 1, S = -2; // in knee function
 * int P = 1; // in underlying shape function
 * LINKAGE= TRUE OR FALSE, need linkage function or not...
 */
public class XPMOP8 extends Problem {
	
//	private boolean LINKAGE = true;
	private boolean LINKAGE = false;

	/**
	 * Creates a default PMOP1 problem (7 variables and 3 objectives)
	 * 
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public XPMOP8(String solutionType) throws ClassNotFoundException {
//		 this(solutionType, 7,2);
//		 this(solutionType, 12,3);//K=5 //n=m+k-1 k=5
//		 this(solutionType, 14,5);//(7,3)
		 this(solutionType, 17,8);//(7,3)
		// this(solutionType, 14,10);//(7,3)
	}

	/**
	 * 
	 * @param solutionType
	 * @param numberOfVariables
	 * @param numberOfObjectives
	 * @throws ClassNotFoundException
	 */
	public XPMOP8(String solutionType, Integer numberOfVariables, Integer numberOfObjectives)
			throws ClassNotFoundException {
		numberOfVariables_ = numberOfVariables.intValue();
		numberOfObjectives_ = numberOfObjectives.intValue();
		numberOfConstraints_ = 0;

		problemName_ = "PMOP8";

		lowerLimit_ = new double[numberOfVariables_];
		upperLimit_ = new double[numberOfVariables_];
		// the decision variables m-1 [0,1], others [0,10]
		for (int var = 0; var < numberOfObjectives_ - 1; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 1.0;
		} // for
		for (int var = numberOfObjectives_-1; var < numberOfVariables_; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 10.0;
		} // for

		if (solutionType.compareTo("BinaryReal") == 0)
			solutionType_ = new BinaryRealSolutionType(this);
		else if (solutionType.compareTo("Real") == 0)
			solutionType_ = new RealSolutionType(this);
		else {
			System.out.println("Error: solution type " + solutionType + " invalid");
			System.exit(-1);
		}
	}

	/**
	 * Evaluates a solution
	 * 
	 * @param solution
	 *            The solution to evaluate
	 * @throws JMException
	 */
	public void evaluate(Solution solution) throws JMException {
		Variable[] gen = solution.getDecisionVariables();

		double[] x = new double[numberOfVariables_];
		double[] f = new double[numberOfObjectives_];
		int k = numberOfVariables_ - numberOfObjectives_ + 1;

		for (int i = 0; i < numberOfVariables_; i++)
			x[i] = gen[i].getValue();

		double g = 0.0, tmp = 0.0, tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0;
		for(int i=0;i<k;i++) {
			if (LINKAGE==true) {// linkage function: linear 
				tmp1 = (1 + Math.cos(0.5 * Math.PI * (i+1)/k))*(x[numberOfObjectives_-1+i]-lowerLimit_[numberOfObjectives_-1+i])
						- x[0]*(upperLimit_[numberOfObjectives_-1+i]-lowerLimit_[numberOfObjectives_-1+i]);
			}else {
				tmp1 = x[numberOfObjectives_-1+i];
			}
			tmp2 += tmp1*tmp1;
			tmp3 += Math.cos(2*Math.PI*tmp1);
		}
		g=-20*Math.exp(-0.2*Math.sqrt(tmp2/k))-Math.exp(tmp3/k)+20 + Math.E; // g_8 function
		
		
		

		double k_ = 1.0;
		int A = 4, B = 1, S = 2; // in knee function
		int P = 1; // in underlying shape function
		double tmp_ =0.0;
		for (int i = 0; i < numberOfObjectives_ - 1; i++) {// knee function k_3 with no degerneration 
			tmp_ = 1 + Math.exp(Math.sin(A*Math.pow(x[i], B)*Math.PI+Math.PI*0.1/(0.2)))/ (Math.pow(2, S) * A);
			k_ =k_*tmp_;
		}

		k_ = k_ /  (numberOfObjectives_ - 1); //  k_3

		for (int i = 0; i < numberOfObjectives_; i++) // objective function :
														// g_8, k_3,
			f[i] = (1.0 + g) * k_;

		for (int i = 0; i < numberOfObjectives_; i++) {// g_8, k_3, l2, h2
			for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
				f[i] *= Math.cos(Math.pow(x[j], P)*Math.PI/2);
			if (i != 0) {
				int aux = numberOfObjectives_ - (i + 1);
				f[i] *= Math.sin(Math.pow(x[aux], P)*Math.PI/2);
			} // if
		} // for

		for (int i = 0; i < numberOfObjectives_; i++)
			solution.setObjective(i, f[i]);
	} // evaluate

}
