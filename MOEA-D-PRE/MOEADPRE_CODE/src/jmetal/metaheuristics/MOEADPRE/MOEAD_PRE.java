package jmetal.metaheuristics.MOEADPRE;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.LinkedList;
import java.util.StringTokenizer;
import jmetal.util.*;

import java.util.Vector;

import jmetal.core.*;
import jmetal.util.PseudoRandom;

public class MOEAD_PRE extends Algorithm {

  private int populationSize_;
  /**
   * Stores the population
   */
  private SolutionSet population_;
  
  /*************************************************************************************************/
  private final double epson=0.2;
  /*************************************************************************************************/
  double[] refer;
  int ngeneration;
  /*************************************************************************************************/
  private final double sita = 5.0;
  /**
   * Z vector (ideal point)
   */
  double[] z_;
  /**
   * Lambda vectors
   */
  //Vector<Vector<Double>> lambda_ ; 
  double[][] lambda_;
  /**
   * T: neighbour size
   */
  int T_;
  /**
   * Neighborhood
   */
  int[][] neighborhood_;
  private int[] frequency_;
  /**
   * delta: probability that parent solutions are selected from neighbourhood
   */
  double delta_; 
  /**
   * nr: maximal number of solutions replaced by each child solution
   */
  int nr_; 
  Solution[] indArray_;
  String functionType_;
  int evaluations_;
  private final int nsize = 1000;// 200   500
  /**
   * Operators
   */
  Operator crossover_;
  Operator mutation_;

  String dataDirectory_;

  /** 
   * Constructor
   * @param problem Problem to solve
   */
  public MOEAD_PRE(Problem problem) {
    super (problem) ;
   
//    functiontype="WEIGHT_SUM"\"_TCHE1"\"AREA"\"PBI"
//            functionType_ = "_TCHE1";
//        functionType_ = "WEIGHT_SUM";
            functionType_ = "PBI";
     System.out.println("Algorithm:"+this.toString());
    System.out.println("Function"+functionType_);
  } // DMOEA

  public SolutionSet execute() throws JMException, ClassNotFoundException {
    int maxEvaluations;

    evaluations_ = 0;
    ngeneration = 0;
    maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations")).intValue();
    populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
//    dataDirectory_ = this.getInputParameter("dataDirectory").toString();
    System.out.println("POPSIZE: "+ populationSize_) ;

    population_ = new SolutionSet(populationSize_); //construct a population
    indArray_ = new Solution[problem_.getNumberOfObjectives()]; //      save obj

    T_ = 10;    // the neighbor
    delta_ = 0.9;
    nr_ = 2;    // 
/*
    T_ = (int) (0.1 * populationSize_);
    delta_ = 0.9;
    nr_ = (int) (0.01 * populationSize_);
*/
    neighborhood_ = new int[populationSize_][T_];   //

    z_   = new double[problem_.getNumberOfObjectives()];  // ieal point
    /*************************************************************************************************/
    refer =new double[problem_.getNumberOfObjectives()]; //
    /*************************************************************************************************/
    
    //lambda_ = new Vector(problem_.getNumberOfObjectives()) ;
    lambda_ = new double[populationSize_][problem_.getNumberOfObjectives()];// 
    crossover_ = operators_.get("crossover"); // default: DE crossover
    mutation_ = operators_.get("mutation");  // default: polynomial mutation
    /*************************************************************************************************/
//    frequency_ = new int[populationSize_];  //
//    for (int i = 0; i < population_.size(); i++) {
//            frequency_[i] = 0;
//        }
    /*************************************************************************************************/
    // STEP 1. Initialization
    // STEP 1.1. Compute euclidean distances between weight vectors and find T

    //    initUniformWeight();    
    initWeight();
    initNeighborhood();     //

    // STEP 1.2. Initialize population
    initPopulation();      

    // STEP 1.3. Initialize z_
//    initIdealPoint();       
    // STEP 2. Update
    do {
      int[] permutation = new int[populationSize_];
      Utils.randomPermutation(permutation, populationSize_);

      for (int i = 0; i < populationSize_; i++) {
//        int n = permutation[i]; // or int n = i;
        int n = i ;
        int type;
        double rnd = PseudoRandom.randDouble();
//        frequency_[n]++;
        // STEP 2.1. Mating selection based on probability
        if (rnd < delta_) // if (rnd < realb)     delta_=0.9
        {
          type = 1;   // neighborhood
        } else {
          type = 2;   // whole population
        }
        Vector<Integer> p = new Vector<Integer>();
        matingSelection(p, n, 2, type);

        // STEP 2.2. Reproduction
        Solution child;
        Solution[] parents = new Solution[3];
        
//        Solution[] parents = new Solution[2];
        parents[0] = population_.get(p.get(0));
        parents[1] = population_.get(p.get(1));
        parents[2] = population_.get(n);

        // Apply DE crossover 
        child =  (Solution) crossover_.execute(new Object[]{parents[2], parents});
//        child = (Solution[]) crossover_.execute(parents);
        // Apply mutation
        mutation_.execute(child);

        // Evaluation
        problem_.evaluate(child);      
        
        evaluations_++;

        // STEP 2.3. Repair. Not necessary

        // STEP 2.4. Update z_
//        updateReference(child);

            // STEP 2.5. Update of solutions
        updateProblem(child, n, type);
      } // for 
      ngeneration++;
      System.out.println("Generations: "+ngeneration);
/************************************************************************************/     
    } while (evaluations_ < maxEvaluations);
    return population_;
  }

public void initWeight() {
                int nobj = problem_.getNumberOfObjectives();
		lambda_ = new double[populationSize_][nobj];
		double []refleck = new double[nobj];
                        initRefer();
			refleckPoint(refer,nobj,refleck);
			refleckPoint(refer,nobj,lambda_[0]);	
		if(nobj<3){
			double steplenth = epson/(populationSize_-1);
			for(int i=1;i<populationSize_;i++){
				lambda_[i][0] = refleck[0]-epson/2+steplenth*i;
				lambda_[i][1] = 1-lambda_[i][0];
			}
		}	
		if(nobj>=3){	
			double[][] bj = new double[nobj][nobj];
				getBj(refleck, nobj, bj);
                        int weightcount = 1;
	        	for(int i=0;i<nobj;i++){
	        		for(int j=0;j<nobj;j++){
	        			lambda_[i+1][j] = bj[i][j];	
	        		}
	        	}
	        	weightcount+=nobj;
			
			int newCount = initBegin(lambda_,weightcount,nobj,populationSize_);
			weightcount+=newCount;
	
			while(weightcount<populationSize_){
				newCount = creat(lambda_,weightcount-newCount,weightcount,populationSize_,nobj);
				weightcount+=newCount;
			}
		}
        print(lambda_, nobj);
    }

	 /**
	  * @param array
	  * @param begin
	  * @param end
	  * @param limit
	  * @param nobj
	  */
	 public int creat(double[][] array,int begin,int end,int limit,int nobj) {
		int k=0;
		 for(int i=0;i<begin;i++){
			for(int j=begin;j<end;j++){
				if(end+k<limit){
				array[end+k] = getlamb(array[i], array[j], nobj);
				k++;
				}else {
					break;
				}
				
			}
		}
		 return k;
	}
	  /**
	   * @param a
	   * @param b
	   * @param nobj
	   * @return 
	   */
	 private double[] getlamb(double[] a,double[] b,int nobj){
	      double[] lamb = new double [nobj];
	      for(int i=0;i<nobj;i++){
	          lamb[i] = (a[i]+b[i])/2;
	      }
	      return lamb;
	  }
	 /**
	  * @param a
	  * @param end
	  * @param b
	  * @param nobj
	  */
	 public int initBegin(double[][] old,int end,int nobj,int limitSize) {
		int k=0;
		 for (int i = 0; i <end; i++) {
			for (int j = i+1; j < end; j++) {
				if (k+end<limitSize) {
					old[k+end] = getlamb(old[i],old[j],nobj);
					k++;
				}else {
					break;
				}
				
			}
		}
		 return k;
	}
	 
	/**
	 * set reference vectors
	 */
	 public void initRefer() {
               int nobj = problem_.getNumberOfObjectives();
               if(nobj==2){
             /**
              * zdt1 : 
              */             
                   //1
//                   refer[0] = 0.3;
//                   refer[1] = 0.45;
                   //2 
                   refer[0] = 0.25;
                   refer[1] = 0.75;     
             /**
              * zdt2 :
              */
                    //1 
//                   refer[0] = 0.6;
//                   refer[1] = 0.64;
                   //2 
//                   refer[0] = 0.7;
//                   refer[1] = 0.8;  
             /**
              * zdt3 :
              */
                    //1 
//                   refer[0] = 0.24;
//                   refer[1] = 0.28;
                   //2 
//                   refer[0] = 0.4;
//                   refer[1] = 0.4;
             /**
              * zdt4 :
              */
                    //1 
//                   refer[0] = 0.3;
//                   refer[1] = 0.45;
                   //2 
//                   refer[0] = 0.5;
//                   refer[1] = 0.5;
             /**
              * zdt6 :
              */
                    //1 
//                   refer[0] = 0.6;
//                   refer[1] = 0.64;
//                   2 
//                   refer[0] = 0.7;
//                   refer[1] = 0.8;
               }
        if(nobj==3){   
        	
             /**
              * dtlz1:
              */
                   //1 
                   refer[0] = 0.1;
                   refer[1] = 0.2;
                   refer[2] = 0.2;
//                   //2 
//                   refer[0] = 0.25;
//                   refer[1] = 0.25;
//                   refer[2] = 0.25;
            /**
              * dtlz2:
              */
                   //1 
//                   refer[0] = 0.4;
//                   refer[1] = 0.8;
//                   refer[2] = 0.45;
//                   //2 
//                   refer[0] = 0.8;
//                   refer[1] = 0.8;
//                   refer[2] = 0.8;
           /**
             * dtlz3:
             */
                   //1
//                   refer[0] = 0.4;
//                   refer[1] = 0.8;
//                   refer[2] = 0.45;
//                   //2 
//                   refer[0] = 0.8;
//                   refer[1] = 0.8;
//                   refer[2] = 0.8;
             /**
             * dtlz4:
             */
                   //1 
//                   refer[0] = 0.5;
//                   refer[1] = 0.5;
//                   refer[2] = 0.7;
//                   //2 
//                   refer[0] = 0.6;
//                   refer[1] = 0.6;
//                   refer[2] = 0.8; 
             /**
             * dtlz5:
             */
                   //1 
//                   refer[0] = 0.4;
//                   refer[1] = 0.4;
//                   refer[2] = 0.82;
//                   //2 
//                   refer[0] = 0.7;
//                   refer[1] = 0.7;
//                   refer[2] = 0.9;    
              /**
             * dtlz6:
             */
                   //1 
//                   refer[0] = 0.3;
//                   refer[1] = 0.3;
//                   refer[2] = 0.9;
//                   //2 
//                   refer[0] = 0.7;
//                   refer[1] = 0.6;
//                   refer[2] = 0.6;    
           
               }if (nobj == 5) {
            	   refer[0] = 0.001;
                   refer[1] = 0.001;
                   refer[2] = 0.001;
                   refer[3] = 0.001;
                   refer[4] = 0.796;
            	   
//            	   refer[0] = 0.7*(1.0/5.0);
//                   refer[1] = 0.7*(1.0/5.0);
//                   refer[2] = 0.7*(1.0/5.0);
//                   refer[3] = 0.7*(1.0/5.0);
//                   refer[4] = 0.7*(1.0/5.0);
               }if (nobj == 10) {
            	   refer[0] = 0.001;
                   refer[1] = 0.001;
                   refer[2] = 0.001;
                   refer[3] = 0.001;
                   refer[4] = 0.001;
                   refer[5] = 0.001;
                   refer[6] = 0.001;
                   refer[7] = 0.001;
                   refer[8] = 0.001;
                   refer[9] = 0.791;
//            	   refer[0] = 0.7*(1.0/10.0);
//                   refer[1] = 0.7*(1.0/10.0);
//                   refer[2] = 0.7*(1.0/10.0);
//                   refer[3] = 0.7*(1.0/10.0);
//                   refer[4] = 0.7*(1.0/10.0);
//                   refer[5] = 0.7*(1.0/10.0);
//                   refer[6] = 0.7*(1.0/10.0);
//                   refer[7] = 0.7*(1.0/10.0);
//                   refer[8] = 0.7*(1.0/10.0);
//                   refer[9] = 0.7*(1.0/10.0);
			}
//               if(nobj==10){
//                    refer[0] = 0.30;
//                    refer[1] = 0.30;
//                    refer[2] = 0.30;
//                    refer[3] = 0.10;
//                    refer[4] = 0.30;
//                    refer[5] = 0.55;
//                    refer[6] = 0.35;
//                    refer[7] = 0.35;
//                    refer[8] = 0.25;
//                    refer[9] = 0.45;
// 
//               }
      
	  }
	  /**
	   * @param refleck
	   * @param nobj
	   * @return 
	   */
	 private void getBj(double[] refleck,int nobj,double[][] bj2){
          double[] bj1 = new double[nobj];
	      for(int i=0;i<nobj;i++){
	          for(int j=0;j<nobj;j++){
	              if(j==i){
	                bj1[i] += Math.pow(1.0-refleck[j], 2);
	              }
	              else{
		            bj1[i] += Math.pow(refleck[j], 2);
	              }
	          }
	          bj1[i] = Math.sqrt(bj1[i]);
	      }
	      /**
	       */
	      for(int i=0;i<nobj;i++){
	          for(int j=0;j<nobj;j++){
	              if(j==i){
	            	bj2[i][j] = epson*1.0/bj1[i]+(1.0-epson/bj1[i])*refleck[j];  
	              }else{
                	bj2[i][j] = (1.0-epson/bj1[i])*refleck[j];
	              }
	          }
	      }
	  }
	/**
	 * @param refer
	 * @param nobj
	 * @param refleck
	 */
	 private void refleckPoint(double []refer,int nobj,double[] refleck){
	      double sum = 0.0;
	      for(int m=0;m<nobj;m++)
	      {
	              sum+=refer[m];
	      }
	      for(int i=0;i<nobj;i++){
	              refleck[i]=refer[i]/sum;            
	          }
	  }
	 /**
	  * @param refleck
	  * @param nobj
	  * @param begin
	  */
	 public void initbegin(double [] refleck,int nobj,double[] begin){
	      for(int i=0;i<nobj;i++){
	          begin[i] = refleck[i]-epson/2;
	      }
	  }
	 /**
	  * @param weiDian
	  * @param begin
	  * @param stepNum
	  * @param step
	  * @param nobj
	  */
	 public void getWeiDian(double[][] weiDian,double[] begin,int stepNum,double step,int nobj){
	      for(int i=0;i<nobj;i++){
	          for(int j=0;j<stepNum;j++){
	              weiDian[i][j] = begin[i]+step*j;
	          }
	      }
	  }
	 /**
	  * 打印权重�?
	  * @param a
	  * @param len
	  */
	 private void print(double a[][],int len){
	      for(int i=0;i<a.length;i++){
	          for(int j=0;j<len;j++){
	              System.out.print(a[i][j]+"\t");
	          }
	          System.out.println();
	      }
	  }

  
  
  
/**********************************************************************************************************************************************/
  /**
   * initUniformWeight
   */
  public void initUniformWeight() {
    if ((problem_.getNumberOfObjectives() == 2) && (populationSize_ <= 300)) {
      for (int n = 0; n < populationSize_; n++) {
        double a = 1.0 * n / (populationSize_ - 1);
        lambda_[n][0] = a;
        lambda_[n][1] = 1-a;
        
      } // for
    } // if
    else {
      String dataFileName;
      dataFileName = "W" + problem_.getNumberOfObjectives() + "D_" +
        populationSize_ + ".dat";
   
      try {
        // Open the file
        FileInputStream fis = new FileInputStream(dataDirectory_ + "/" + dataFileName);
        InputStreamReader isr = new InputStreamReader(fis);
        BufferedReader br = new BufferedReader(isr);

        int numberOfObjectives = 0;
        int i = 0;
        int j = 0;
        String aux = br.readLine();
        while (aux != null) {
          StringTokenizer st = new StringTokenizer(aux);
          j = 0;
          numberOfObjectives = st.countTokens();
          while (st.hasMoreTokens()) {
            double value = (new Double(st.nextToken())).doubleValue();
            lambda_[i][j] = value;
            //System.out.println("lambda["+i+","+j+"] = " + value) ;
            j++;
          }
          aux = br.readLine();
          i++;
        }
        br.close();
      } catch (Exception e) {
        System.out.println("initUniformWeight: failed when reading for file: " + dataDirectory_ + "/" + dataFileName);
        e.printStackTrace();
      }
    } // else

    //System.exit(0) ;
  } // initUniformWeight
  
/**********************************************************************************************************************************************/  
  /**
   * 
   */
  public void initNeighborhood() {
    double[] x = new double[populationSize_];
    int[] idx = new int[populationSize_];

    for (int i = 0; i < populationSize_; i++) {
      // calculate the distances based on weight vectors
      for (int j = 0; j < populationSize_; j++) {
        x[j] = Utils.distVector(lambda_[i], lambda_[j]);
        //x[j] = dist_vector(population[i].namda,population[j].namda);
        idx[j] = j;
      //System.out.println("x["+j+"]: "+x[j]+ ". idx["+j+"]: "+idx[j]) ;
      } // for

      // find 'niche' nearest neighboring subproblems
      Utils.minFastSort(x, idx, populationSize_, T_);
      //minfastsort(x,idx,population.size(),niche);

      for (int k = 0; k < T_; k++) {
        neighborhood_[i][k] = idx[k];
      //System.out.println("neg["+i+","+k+"]: "+ neighborhood_[i][k]) ;
      }
    } // for
  } // initNeighborhood

  /**
   * 
   */
  public void initPopulation() throws JMException, ClassNotFoundException {
    for (int i = 0; i < populationSize_; i++) {
      Solution newSolution = new Solution(problem_);

      problem_.evaluate(newSolution);
      evaluations_++;
      population_.add(newSolution) ;
    } // for
  } // initPopulation

  /**
   * 
   */
  void initIdealPoint() throws JMException, ClassNotFoundException {
    for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
      z_[i] = 1.0e+30;
      indArray_[i] = new Solution(problem_);
      problem_.evaluate(indArray_[i]);
      evaluations_++;
    } // for

    for (int i = 0; i < populationSize_; i++) {
      updateReference(population_.get(i));
    } // for
  } // initIdealPoint

  public void matingSelection(Vector<Integer> list, int cid, int size, int type) {
    // list : the set of the indexes of selected mating parents
    // cid  : the id of current subproblem
    // size : the number of selected mating parents
    // type : 1 - neighborhood; otherwise - whole population
    int ss;
    int r;
    int p;

    ss = neighborhood_[cid].length;
    while (list.size() < size) {
      if (type == 1) {
        r = PseudoRandom.randInt(0, ss - 1);
        p = neighborhood_[cid][r];
      //p = population[cid].table[r];
      } else {
        p = PseudoRandom.randInt(0, populationSize_ - 1);
      }
      boolean flag = true;
      for (int i = 0; i < list.size(); i++) {
        if (list.get(i) == p) // p is in the list
        {
          flag = false;
          break;
        }
      }

      //if (flag) list.push_back(p);
      if (flag) {
        list.addElement(p);
      }
    }
  } // matingSelection

  /**
   * 
   * @param individual
   */
  void updateReference(Solution individual) {
    for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
      if (individual.getObjective(n) < z_[n]) {
        z_[n] = individual.getObjective(n);

        indArray_[n] = individual;
      }
    }
  } // updateReference

  /**
   * @param individual
   * @param id
   * @param type
   */
  void updateProblem(Solution indiv, int id, int type) {
    // indiv: child solution
    // id:   the id of current subproblem
    // type: update solutions in - neighborhood (1) or whole population (otherwise)
    int size;
    int time;

    time = 0;

    if (type == 1) {
      size = neighborhood_[id].length;
    } else {
      size = population_.size();
    }
    int[] perm = new int[size];

    Utils.randomPermutation(perm, size);

    for (int i = 0; i < size; i++) {
      int k;
      if (type == 1) {
        k = neighborhood_[id][perm[i]];
      } else {
        k = perm[i];      // calculate the values of objective function regarding the current subproblem
      }
      double f1, f2;

      f1 = fitnessFunction(population_.get(k), lambda_[k]);
      f2 = fitnessFunction(indiv, lambda_[k]);

      if (f2 < f1) {
        population_.replace(k, new Solution(indiv));
        //population[k].indiv = indiv;
        time++;
      }
      // the maximal number of solutions updated is not allowed to exceed 'limit'
      if (time >= nr_) {
        return;
      }
    }
  } // updateProblem

  /**
   * @param individual
   * @param lambda
   * @return 
   */
    double fitnessFunction(Solution individual,double[] lambda){
      double fitness = 0.0;
        if (functionType_.equals("_TCHE1")) {
            fitness = _TCHE1(individual,lambda);
        }else if(functionType_.equals("WEIGHT_SUM")){
            fitness = WEIGHT_SUM(individual,lambda);
        }else if(functionType_.equals("PBI")){
            fitness = PBI(individual,lambda);
        }       
      return fitness;
  }
  
  /**
   * functionType = "_TCHE1"
   * @param individual
   * @param lambda
   * @return 
   */
  double _TCHE1(Solution individual, double[] lambda) {
    double fitness;
    fitness = 0.0;
    if (functionType_.equals("_TCHE1")) {
      double maxFun = -1.0e+30;

      for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
        double diff = Math.abs(individual.getObjective(n) - z_[n]);

        double feval;
        if (lambda[n] == 0) {
          feval = 0.0001 * diff;
        } else {
          feval = diff * lambda[n];
        }
        if (feval > maxFun) {
          maxFun = feval;
        }
      } // for

      fitness = maxFun;
    } // if
    else {
      System.exit(-1);
    }
    return fitness;
  } // fitnessEvaluation
  
 /**
  * functionType="WEIGHT_SUM"
  * @param individual
  * @param lambda
  * @return 
  */
  double WEIGHT_SUM(Solution individual, double[] lambda){
      double fitness = 0;
      if (functionType_.equals("WEIGHT_SUM")) {
          for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
              fitness+=individual.getObjective(i)*lambda[i];
          }
      }else{
          System.exit(-1);
      }
      return fitness;
  }
  
 /**
  * functionType = "PBI"
  * @param individual
  * @param lambda
  * @return 
  */
  double PBI(Solution individual,double[] lambda){
      double fitness = 0.0;
      double d1 = 0.0;
      double d2 = 0.0;
      double sum=0.0;
      int nobj = problem_.getNumberOfObjectives();
      //double[] a = null;
      for(int i=0;i<nobj;i++){
          sum +=Math.pow(lambda[i], 2); 
      }
      sum = Math.sqrt(sum);
      for(int j = 0;j<nobj;j++){
          d1+=Math.abs((individual.getObjective(j) -z_[j])*lambda[j]/sum);
      }
      for(int k = 0;k<nobj;k++){
          d2+=Math.pow((individual.getObjective(k)-(z_[k]+d1*lambda[k])),2);
     }
      d2 = Math.sqrt(d2);
      fitness = d1+sita*d2;
      return fitness;
  }
  
} // MOEAD

