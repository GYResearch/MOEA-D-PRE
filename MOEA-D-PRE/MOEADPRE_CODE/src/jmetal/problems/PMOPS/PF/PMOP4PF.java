
package jmetal.problems.PMOPS.PF;

import jmetal.core.*;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;
import jmetal.util.Configuration.*;

/** 
 * Class representing problem PMOP4: g_4, l_1, k_4, h_2
 */
public class PMOP4PF extends Problem {   
 /** 
  * Creates a default PMOP1 problem (7 variables and 3 objectives)
  * @param solutionType The solution type must "Real" or "BinaryReal". 
  */
  public PMOP4PF(String solutionType) throws ClassNotFoundException {
//	  this(solutionType, 7,2);
	  this(solutionType, 7,3);//K=5    
//      this(solutionType, 12,3);//K=5    //n=m+k-1   k=5
//    this(solutionType, 9,5);//(7,3)
//    this(solutionType, 14,10);//(7,3)
  } 
    
  /**
   * 
   * @param solutionType
   * @param numberOfVariables
   * @param numberOfObjectives
   * @throws ClassNotFoundException
   */
  public PMOP4PF(String solutionType, 
               Integer numberOfVariables, 
  		         Integer numberOfObjectives) throws ClassNotFoundException {
    numberOfVariables_  = numberOfVariables.intValue();
    numberOfObjectives_ = numberOfObjectives.intValue();
    numberOfConstraints_= 0;
    
    problemName_        = "PMOP4";
        
    lowerLimit_ = new double[numberOfVariables_];
    upperLimit_ = new double[numberOfVariables_];       
    // the decision variables m-1 [0,1], others [0,10]
    for (int var = 0; var < numberOfObjectives_-1; var++){
      lowerLimit_[var] = 0.0;
      upperLimit_[var] = 1.0;
    } //for
    for (int var = numberOfObjectives_; var < numberOfVariables_; var++){
        lowerLimit_[var] = 0.0;
        upperLimit_[var] = 10.0;
      } //for
        
    
    
    if (solutionType.compareTo("BinaryReal") == 0)
    	solutionType_ = new BinaryRealSolutionType(this) ;
    else if (solutionType.compareTo("Real") == 0)
    	solutionType_ = new RealSolutionType(this) ;
    else {
    	System.out.println("Error: solution type " + solutionType + " invalid") ;
    	System.exit(-1) ;
    }            
  }            
 
  /** 
  * Evaluates a solution 
  * @param solution The solution to evaluate
   * @throws JMException 
  */    
  public void evaluate(Solution solution) throws JMException {
    Variable[] gen  = solution.getDecisionVariables();
                
    double [] x = new double[numberOfVariables_];
    double [] f = new double[numberOfObjectives_];
    int k = numberOfVariables_ - numberOfObjectives_ + 1;
        
    for (int i = 0; i < numberOfVariables_; i++)
      x[i] = gen[i].getValue();
        
    double g = 0.0,tmp1=0.0,temp=0.0;
    for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++){//The first variable is x[0]
    	tmp1=(1+0.1*(i-numberOfObjectives_+2)/(0.1*k))*(x[i]-lowerLimit_[i])-x[0]*(upperLimit_[i]-lowerLimit_[i]); // l_1 function: transform x[i] to tmp1
    	g+= (tmp1-0.5)*(tmp1-0.5)-Math.cos(20*Math.PI*(tmp1-0.5)); //g_4 function
    }
    g=100*(k+g);  //g_4 function
    
    double k_=0.0;
    int A=6,B=2,S=-2; // in knee function
    int P=1; // in underlying shape function
    
    for(int i=0;i<numberOfObjectives_-1;i++){// knee function k_4 with no degerneration
    	double tmp=2+Math.abs(Math.sin(A*Math.pow(x[i], P))-Math.cos(A*Math.pow(x[i], P)-Math.PI/4))/(Math.pow(2, S));
    	k_+=tmp;
    }
    
    k_ = 0.1 * k_ / (0.1 * (numberOfObjectives_ - 1));// knee function k_4
    
    
    for (int i = 0; i < numberOfObjectives_; i++) // objective function : g_4, k_4,
      f[i] = (1.0 + g) * k_;
        
        
    for (int i = 0; i < numberOfObjectives_; i++) {// g_4, k_4, l1, h2
		for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
			f[i] *= Math.cos(Math.pow(x[j], P)*Math.PI/2);
		if (i != 0) {
			int aux = numberOfObjectives_ - (i + 1);
			f[i] *= Math.sin(Math.pow(x[aux], P)*Math.PI/2);
		} // if
	} // for

    
    for (int i = 0; i < numberOfObjectives_; i++)
      solution.setObjective(i,f[i]);        
  } // evaluate   
  
}

