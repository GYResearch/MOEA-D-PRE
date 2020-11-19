package jmetal.problems.FMOPS;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

public class FMOP3 extends Problem { 
	
	private boolean LINKAGE = false;
 /** 
  * Creates a default PMOP1 problem (7 variables and 3 objectives)
  * @param solutionType The solution type must "Real" or "BinaryReal". 
  */
  public FMOP3(String solutionType) throws ClassNotFoundException {
	  this(solutionType, 7,2);
//	  this(solutionType, 7,3);//K=5    
//      this(solutionType, 12,3);//K=5    //n=m+k-1   k=5
//    this(solutionType, 9,5);//(7,3)
//    this(solutionType, 14,10);//(7,3)
  } 
  
  
  public FMOP3(String solutionType, Integer numberOfVariables, Integer numberOfObjectives) throws ClassNotFoundException {
	  numberOfVariables_  = numberOfVariables.intValue();
	  numberOfObjectives_ = numberOfObjectives.intValue();
	  numberOfConstraints_= 0;

	  problemName_        = "FMOP3";

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


@Override
public void evaluate(Solution solution) throws JMException {
	// TODO Auto-generated method stub
			Variable[] gen = solution.getDecisionVariables();
			
			double[] x = new double[numberOfVariables_];
			double[] f = new double[numberOfObjectives_];
			int k = numberOfVariables_ - numberOfObjectives_ + 1;
			
			for (int i = 0; i < numberOfVariables_; i++)
				x[i] = gen[i].getValue();
			
			double g = 0.0;
			double tmp1 = 0.0, temp = 0.0;
			
			for(int i=0;i<k;i++) {
				if (LINKAGE==true) {// linkage function: linear l1
					tmp1 = (1+ (i+1)/k)*(x[numberOfObjectives_-1+i]-lowerLimit_[numberOfObjectives_-1+i])
							- x[0]*(upperLimit_[numberOfObjectives_-1+i]-lowerLimit_[numberOfObjectives_-1+i]);
				}else {
					tmp1 = x[numberOfObjectives_-1+i];
				}
				
		    	g+= tmp1*tmp1-10*Math.cos(4*Math.PI*tmp1); //g_3 function
			}
			g+=1+10*(k);  //g_3 function
			
//			g=0; // to get the solutions on the PoF.
			
			double k_=0.0;
		    int A=6,B=1,S=-1; // in knee function
		    int P=1; // in underlying shape function
		    
		    for(int i=0;i<numberOfObjectives_-1;i++){// knee function k_3 with no degerneration
		    	k_+=1+Math.exp(Math.sin(A*Math.PI*Math.pow(x[i], B)+Math.PI/2))/(Math.pow(2, S)*A);
		    }
		    
		    k_ =   k_ / (numberOfObjectives_ - 1);// knee function k_3
		    
		    for (int i = 0; i < numberOfObjectives_; i++) // objective function : g_3, k_3,
		        f[i] = (1.0 + g) * k_;
		    
		    for (int i = 0; i < numberOfObjectives_; i++) {// g_3, k_3, l1, h3
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					f[i] *= 1-Math.cos(Math.pow(x[j], P)*Math.PI/2);
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					f[i] *= 1-Math.sin(Math.pow(x[aux], P)*Math.PI/2);
				} // if
			} // for

		    
		    for (int i = 0; i < numberOfObjectives_; i++)
		      solution.setObjective(i,f[i]); 
			
}

}
