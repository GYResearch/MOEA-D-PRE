
package jmetal.problems.robust;

import jmetal.core.*;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;
import jmetal.util.Configuration.*;

/** 
 * Class representing problem P1 
 */
public class R1 extends Problem {   
 /** 
  * Creates a default P1 problem (7 variables and 3 objectives)
  * @param solutionType The solution type must "Real" or "BinaryReal". 
  */
  public R1(String solutionType) throws ClassNotFoundException {
	  this(solutionType, 4,2);
  }
    
  /** 
  * Creates a P1 problem instance
  * @param numberOfVariables Number of variables
  * @param numberOfObjectives Number of objective functions
  * @param solutionType The solution type must "Real" or "BinaryReal". 
  */
  public R1(String solutionType, 
               Integer numberOfVariables, 
  		         Integer numberOfObjectives) throws ClassNotFoundException {
    numberOfVariables_  = numberOfVariables.intValue();
    numberOfObjectives_ = numberOfObjectives.intValue();
    numberOfConstraints_= 0;
    
    problemName_        = "R1";
        
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
  public void evaluate(Solution solution,double disturbation) throws JMException {
    Variable[] gen  = solution.getDecisionVariables();
                
    double [] x = new double[numberOfVariables_];
    double [] f = new double[numberOfObjectives_];
    
    for (int i = 0; i < numberOfVariables_; i++)
        x[i] = gen[i].getValue()+disturbation;
    
    
    double g=0.0;
    for(int i=1;i<numberOfVariables_;i++){
    	g+=10+x[i]*x[i]-10*Math.cos(4*Math.PI*x[i]);
    }
      
    f[0] = x[0];  
    
    double a=1, b=1;
    double h=0.0,s=0.0;
    h = 1-x[0]*x[0];
    s = a/(0.2+x[0])+b*x[0]*x[0];
    
    double y = 0.0;
    	f[1]=y= h+g*s;
    	solution.setObjective(0,f[0]);   
        solution.setObjective(1,f[1]); 
  } // evaluate   

@Override
public void evaluate(Solution solution) throws JMException {
	// TODO Auto-generated method stub
	
}
  
}

