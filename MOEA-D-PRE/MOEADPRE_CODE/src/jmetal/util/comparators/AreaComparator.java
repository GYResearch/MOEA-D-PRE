/**
 * a new comparator created by lursonkj:
 */
package jmetal.util.comparators;
import jmetal.core.Solution;
import java.util.Comparator;

/**
 *
 * @author lursonkj
 */
public class AreaComparator implements Comparator{
    
//    public int popsize;//传入种群大小。
//
//    public AreaComparator(int popsize) {
//        this.popsize = popsize;
//    }
    
  /** 
   * stores a comparator for check the OverallConstraintComparator
   */
    private static final Comparator overallConstraintViolationComparator_ =
                              new OverallConstraintViolationComparator();
     
  /**
   * cal_sum:计算目标值之和。
   * @param a
   * @return sum
   */
   
    public static double cal_sum(Solution a){
        int i = 0,j = 0;
        double sum=0.0,temp=0.0;
        for(;i<a.numberOfObjectives();i++){
            temp = a.getObjective(i);
            sum = sum+temp;
        }
        return sum;
    }
    
    /**
     * 新发现：SUM在求面积支配的时候是没有用的，等比例时最终为定值
    * cal_ind_area:求对应的个体的面积。
    * @param a
    * @return 
    */    
    
    public static double cal_ind_area(Solution a){                                         
        double PI = 3.1415926;
        double tempa=0.0,area=0.0;
//        double sum = 0.0 ;
        double lastx,lasty;
        int nobj = a.numberOfObjectives();
        double x[] = new double[nobj];
        double y[] = new double[nobj];
        int i,j;
//        sum = cal_sum(a);
        for(i=0;i<nobj;i++){
            x[i] = a.getObjective(i)*Math.cos(2*PI*i/nobj);
            y[i] = a.getObjective(i)*Math.sin(2*PI*i/nobj);			
        }
         for(j = 0;j<nobj-1;j++){
                tempa  = tempa + x[j]*y[j+1]-x[j+1]*y[j];
         }
	area = (tempa+x[nobj-1]*y[0]-y[nobj-1]*x[0])/2.0;
	return area;
    }

 
   /**
    * Compares two solutions.
    * @param object1 Object representing the first <code>Solution</code>.
    * @param object2 Object representing the second <code>Solution</code>.
    * @return -1, or 0, or 1 if solution1 dominates solution2, both are 
    * non-dominated, or solution1  is dominated by solution22, respectively.
    */   
    @Override
    public int compare(Object object1, Object object2) {   
        double temp1,temp2;
        double epson ;
	double tempa,tempb ;
        temp1 = cal_ind_area((Solution)object1);
        temp2 = cal_ind_area((Solution)object2);
        epson = Math.abs(temp1-temp2)/100;  //100 = popsize
        if (object1==null)
            return 1;
        else if (object2 == null)
             return -1;
        
	tempa = cal_ind_area((Solution)object1);
	tempb = cal_ind_area((Solution)object2);
	if(tempa<(tempb-epson))
		{
			return 1;
		}
        else if((tempa-epson)>tempb){
			return -1;
		}
	else {
            return 0;
        }
        
    }     
}
    


   