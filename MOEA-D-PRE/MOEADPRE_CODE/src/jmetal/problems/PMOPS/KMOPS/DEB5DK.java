package jmetal.problems.PMOPS.KMOPS;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

public class DEB5DK extends Problem{
	private static final long serialVersionUID = 2L;
	private int K= 3;

	public DEB5DK(String solutionType) throws ClassNotFoundException {
		 this(solutionType, 14,5);
	}

	public DEB5DK(String solutionType, Integer numberOfVariables, Integer numberOfObjectives) throws ClassNotFoundException{
		// TODO Auto-generated constructor stub
		numberOfVariables_ = numberOfVariables.intValue();
		numberOfObjectives_ = numberOfObjectives.intValue();
		numberOfConstraints_ = 0;

		problemName_ = "DEB5DK";

		lowerLimit_ = new double[numberOfVariables_];
		upperLimit_ = new double[numberOfVariables_];
		// the decision variables m-1 [0,1], others [0,10]
		for (int var = 0; var < numberOfVariables_; var++) {
			lowerLimit_[var] = 0;
			upperLimit_[var] = 1;
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
	
	@Override
	public void evaluate(Solution solution) throws JMException {
		// TODO Auto-generated method stub
				Variable[] gen = solution.getDecisionVariables();

				double[] x = new double[numberOfVariables_];
				double[] f = new double[numberOfObjectives_];
				
				for (int i = 0; i < numberOfVariables_; i++)
					x[i] = gen[i].getValue();
				
				double g = 0.0,r1=0.0,r2=0.0,r3=0.0,r4=0.0,r=0.0;
				for(int i=1;i<numberOfVariables_;i++) {
					g+=x[i];
				}
				g=1+9*g/(numberOfVariables_-1);
				
				r1=5+10*(x[0]-0.5)*(x[0]-0.5)+4*Math.cos(2*K*Math.PI*x[0])/K;
				r2=5+10*(x[1]-0.5)*(x[1]-0.5)+4*Math.cos(2*K*Math.PI*x[1])/K;
				r3=5+10*(x[2]-0.5)*(x[2]-0.5)+4*Math.cos(2*K*Math.PI*x[2])/K;
				r4=5+10*(x[3]-0.5)*(x[3]-0.5)+4*Math.cos(2*K*Math.PI*x[3])/K;
				r = (r1+r2+r3+r4)/4;
				
				f[0]= g*r*Math.sin(Math.PI*x[0]/2)*Math.sin(Math.PI*x[1]/2)*Math.sin(Math.PI*x[2]/2)*Math.sin(Math.PI*x[3]/2);
				f[1]= g*r*Math.sin(Math.PI*x[0]/2)*Math.sin(Math.PI*x[1]/2)*Math.sin(Math.PI*x[2]/2)*Math.cos(Math.PI*x[3]/2);
				f[2]= g*r*Math.sin(Math.PI*x[0]/2)*Math.sin(Math.PI*x[1]/2)*Math.cos(Math.PI*x[2]/2);
				f[3]= g*r*Math.sin(Math.PI*x[0]/2)*Math.cos(Math.PI*x[1]/2);
				f[4]= g*r*Math.cos(Math.PI*x[0]/2);
				
				
				for (int i = 0; i < numberOfObjectives_; i++)
					solution.setObjective(i, f[i]);
	}
	
	
	
}
