package jmetal.util.comparators;

import java.util.Comparator;

import jmetal.core.Solution;

public class AngleComparator implements Comparator{
	
	/** 
	   * stores a comparator for check the OverallConstraintComparator
	   */
	  private static final Comparator overallConstraintViolationComparator_ =
	                              new OverallConstraintViolationComparator();

	@Override
	public int compare(Object object1, Object object2){
		// TODO Auto-generated method stub
		if (object1==null)
		      return 1;
		    else if (object2 == null)
		      return -1;
		
		Solution solution1 = (Solution)object1;
	    Solution solution2 = (Solution)object2;

	    int dominate1 ; // dominate1 indicates if some objective of solution1 
	                    // dominates the same objective in solution2. dominate2
	    int dominate2 ; // is the complementary of dominate1.

	    dominate1 = 0 ; 
	    dominate2 = 0 ;
	    
	    int flag; //stores the result of the comparison

	    if (solution1.getOverallConstraintViolation()!= 
	        solution2.getOverallConstraintViolation() &&
	       (solution1.getOverallConstraintViolation() < 0) ||         
	       (solution2.getOverallConstraintViolation() < 0)){            
	      return (overallConstraintViolationComparator_.compare(solution1,solution2));
	    }
	                                                
	    // Equal number of violated constraints. Applying a dominance Test then
	    double value1, value2;
	    double theta = Math.PI/3; // pi/3.  
	    
	    // -1 : solution1 better; 1 solution2 better; 0 equal.
	    
	    // AB
	    int objs = solution1.numberOfObjectives();
	    double[] standard = new double[objs];
	    double[] AB = new double[objs];
	    // BA
	    double[] BA = new double[objs];  
	    for (int i=0;i<objs;i++) {
	    	AB[i] = solution2.getObjective(i) - solution1.getObjective(i);
	    	standard[i] = 1;
	    	BA[i] = -AB[i];
	    }
	    
	    if(calcAngle(AB, standard)<theta) { // solution 1 dominates
	    	return -1;
	    }else if (calcAngle(BA, standard)<theta) { // solution 2 dominates
			return 1;
		}else {
			return 0; 
		}
	    
} // compare
	
	/**
	 * Calc the angel between two vectors
	 * @param A
	 * @param standard
	 * @return
	 */
	public double calcAngle(double[] A, double[] standard){
		int size = A.length;
		double lenA = 0.0, lenSta = 0.0;
		double theta = 0.0;
		
		for(int i=0;i<size;i++) {
			lenA += A[i]*A[i];
			lenSta += standard[i]*standard[i];
		}
		lenA = Math.sqrt(lenA);
		lenSta = Math.sqrt(lenSta);
		
		double tmp = 0.0;
		for(int j=0;j<size;j++) {
			tmp += A[j]*standard[j];
		}
		
		theta = Math.acos(tmp/(lenA*lenSta));
		
		return theta;	
	}
	

}
