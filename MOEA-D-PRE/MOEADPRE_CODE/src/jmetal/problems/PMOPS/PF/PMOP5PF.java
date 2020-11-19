
package jmetal.problems.PMOPS.PF;

import jmetal.core.*;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;
import jmetal.util.Configuration.*;

/** 
 * Class representing problem PMOP5: g_5, l_1, k_5, h_1
 */
public class PMOP5PF extends Problem {   
 /** 
  * Creates a default PMOP1 problem (7 variables and 3 objectives)
  * @param solutionType The solution type must "Real" or "BinaryReal". 
  */
  public PMOP5PF(String solutionType) throws ClassNotFoundException {
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
  public PMOP5PF(String solutionType, 
               Integer numberOfVariables, 
  		         Integer numberOfObjectives) throws ClassNotFoundException {
    numberOfVariables_  = numberOfVariables.intValue();
    numberOfObjectives_ = numberOfObjectives.intValue();
    numberOfConstraints_= 0;
    
    problemName_        = "PMOP5";
        
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
        
    double g = 0.0,tmp1=0.0,tmp2=0.0,temp=0.0; // g_5: use numberOfVariables_-1 not numberOfVariables_  
    for (int i = numberOfVariables_ - k; i < numberOfVariables_-1; i++){//The first variable is x[0]
    	tmp1=(1+0.1*(i-numberOfObjectives_+2)/(0.1*k))*(x[i]-lowerLimit_[i])-x[0]*(upperLimit_[i]-lowerLimit_[i]); // l_1 function : transform x[i] to tmp1
    	tmp2=(1+0.1*((i-numberOfObjectives_+2)+1)/(0.1*k))*(x[i+1]-lowerLimit_[i+1])-x[0]*(upperLimit_[i+1]-lowerLimit_[i+1]); // l_1 function : transform x[i+1] to tmp2
    	temp= 100*(tmp1*tmp1-tmp2)*(tmp1*tmp1-tmp2)+(tmp1-1)*(tmp1-1); //g_5 function
    }
     g=temp;
     
    double k_=0.0;
    int A=6,B=2,S=-2,l=6; // in knee function  l>=3
    int P=1; // in underlying shape function
    for(int i=0;i<numberOfObjectives_-1;i++){// knee function k_5 with no degerneration
    	double tmp=2+Math.min(Math.sin(2*A*Math.pow(x[i], P)*Math.PI), Math.cos(2*A*Math.pow(x[i], P)*Math.PI-Math.PI/l))/(Math.pow(2, S)*A);
    	k_+=tmp;
    }
    k_=0.1*k_/(0.1*(numberOfObjectives_-1)); // k_5
    
    
    for (int i = 0; i < numberOfObjectives_; i++) // objective function : g_5, k_5,
      f[i] = (1.0 + g) * k_;
        
    
    for (int i = 0; i < numberOfObjectives_; i++){// g_5, k_5, l1, h1
      for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)            
        f[i] *= Math.pow(x[j],P);                
        if (i != 0){
          int aux = numberOfObjectives_ - (i + 1);
          f[i] *= 1 - Math.pow(x[aux],P);
        } //if
    }//for
        
    for (int i = 0; i < numberOfObjectives_; i++)
      solution.setObjective(i,f[i]);        
  } // evaluate   
  
}

