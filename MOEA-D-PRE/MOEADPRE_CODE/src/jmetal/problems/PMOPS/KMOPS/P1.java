
package jmetal.problems.PMOPS.KMOPS;

import jmetal.core.*;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;
import jmetal.util.Configuration.*;

/** 
 * Class representing problem P1 
 */
public class P1 extends Problem {   
 /** 
  * Creates a default P1 problem (7 variables and 3 objectives)
  * @param solutionType The solution type must "Real" or "BinaryReal". 
  */
  public P1(String solutionType) throws ClassNotFoundException {
	  this(solutionType, 4,2);
  }
    
  /** 
  * Creates a P1 problem instance
  * @param numberOfVariables Number of variables
  * @param numberOfObjectives Number of objective functions
  * @param solutionType The solution type must "Real" or "BinaryReal". 
  */
  public P1(String solutionType, 
               Integer numberOfVariables, 
  		         Integer numberOfObjectives) throws ClassNotFoundException {
    numberOfVariables_  = numberOfVariables.intValue();
    numberOfObjectives_ = numberOfObjectives.intValue();
    numberOfConstraints_= 0;
    
    problemName_        = "P1";
        
    lowerLimit_ = new double[numberOfVariables_];
    upperLimit_ = new double[numberOfVariables_];       
    // the decision variables 
    for (int var = 0; var < numberOfVariables_; var++){
      lowerLimit_[var] = 0.0;
      upperLimit_[var] = 1.0;
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
    
    for (int i = 0; i < numberOfVariables_; i++)
        x[i] = gen[i].getValue();
    
    
    double g=0.0;
    for(int i=1;i<numberOfVariables_;i++){
    	g+=(x[i]-0.5)*(x[i]-0.5);
    }
    g=1+g;
    
    f[0] = x[0];  
    
    double y = 0.0;
    if((x[0]>=0)&&(x[0]<=0.5/1.4)){
    	y= 24-Math.sqrt(20-(14*x[0]-5)*(14*x[0]-5));
    	f[1]=(y-3.5)*g/20.5;
    	solution.setObjective(0,f[0]);   
        solution.setObjective(1,f[1]); 
    }else if ((x[0]>0.5/1.4)&&(x[0]<=0.6/1.4)) {
    	y= 18+Math.sqrt(1-(14*x[0]-5)*(14*x[0]-5));
    	f[1]=(y-3.5)*g/20.5;
    	solution.setObjective(0,f[0]);   
        solution.setObjective(1,f[1]); 
	}else if ((x[0]>0.6/1.4)&&(x[0]<=0.7/1.4)) {
    	y= 18-Math.sqrt(1-(14*x[0]-7)*(14*x[0]-7));
    	f[1]=(y-3.5)*g/20.5;
    	solution.setObjective(0,f[0]);   
        solution.setObjective(1,f[1]); 
	}else if ((x[0]>0.7/1.4)&&(x[0]<=0.8/1.4)) {
    	y= 59-84*x[0];
    	f[1]=(y-3.5)*g/20.5;
    	solution.setObjective(0,f[0]);   
        solution.setObjective(1,f[1]); 
	}else if ((x[0]>0.8/1.4)&&(x[0]<=1.0/1.4)) {
    	y= -3.5*x[0]+13;
    	f[1]=(y-3.5)*g/20.5;
    	solution.setObjective(0,f[0]);   
        solution.setObjective(1,f[1]); 
	}else if ((x[0]>1.0/1.4)&&(x[0]<=1.1/1.4)) {
    	y= 50.5-56*x[0];
    	f[1]=(y-3.5)*g/20.5;
    	solution.setObjective(0,f[0]);   
        solution.setObjective(1,f[1]); 
	}else if((x[0]>1.1/1.4)&&(x[0]<=1)){
    	y= 17.5-14*x[0];
    	f[1]=(y-3.5)*g/20.5;
    	solution.setObjective(0,f[0]);   
        solution.setObjective(1,f[1]); 
	}
    
  } // evaluate   
  
}

