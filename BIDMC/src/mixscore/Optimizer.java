package mixscore;

import java.util.Arrays;
import java.util.stream.IntStream;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.FastMath;

public class Optimizer {
	static final int maxEvalualtion = 100000;
	static final double maxOddsRatio = 20; //maximum odds ratio.
	static double maxOmiga = 50;
	static double sigmaF = 0.001; //for allele frequency. 0.001
	static double sigmaR = 0.01; //for odds ratio. 0.01
	static double freLeft = 0.00001;
	static double freRight = 0.99999;
	
//	static final CMAESOptimizer optimizer = new CMAESOptimizer(maxEvalualtion, 0, true, 10, 0, new MersenneTwister(), false, null);

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new Optimizer().testAPP();
	}
	
	public void testAPP() {
		MultivariateFunction f = new MultivariateFunction() {
			@Override
			public double value(double[] x) {
				double xs = x[0] * x[0];
				return (100 * FastMath.pow((x[1] -xs), 2)) + FastMath.pow(1-x[0], 2);
			}
		};
		
		BOBYQAOptimizer optim = new BOBYQAOptimizer(4);
		PointValuePair results = optim.optimize(new MaxEval(100000),
								new ObjectiveFunction(f),
								new InitialGuess(new double[]{0,0}),
								GoalType.MINIMIZE,
								new SimpleBounds(new double[]{-10, -10}, new double[]{10,10})
								);
		
		System.out.println(optim.computeObjectiveValue(results.getPoint()));
		System.out.println(optim.getIterations());
		System.out.println(optim.getEvaluations());
		System.out.println(results.getValue());
		System.out.println(Arrays.toString(results.getPoint()));

	}
	
	
	
	/**
	 * optimize totally four parameters. pa0, pe0, ra, re
	 * @param CaseGeno      Case genotype count.
	 * @param ControlGeno   Control genotype count.
	 * @param initial, set inital if has, otherwise optimize by program itself.
	 * @return              Estimated parameters: logLikelyhood, pa0, pe0, ra, re
	 * 
	 * Genotype count format.
	 *   	2	1	0 Copy of ref. allele 
	 *    RR	RV	VV
	 * AA
	 * AE
	 * EE
	 */
	public static double[] fullOptimizer(int[][] CaseGeno, int[][] ControlGeno, double[]... initial) {
		
		double[] case_coeff = coefficient(CaseGeno);
		double[] ctrl_coeff = coefficient(ControlGeno);
		
		MultivariateFunction f = new MultivariateFunction() {
			@Override
			/**
			 * All four parameters.
			 * pa0: ref. allele frequency in Africa population.
			 * pe0: ref. allele frequency in European population.
			 * ra:  Odds ratio in Africa population.
			 * re:  Odds ratio in European population. 
			 * parameters: pa0, pe0, ra, re
			 */
			public double value(double[] parameters) {  // 4 parameters: pa0, pe0, ra, re
				double pa0 = parameters[0];       //control 0.
				double pe0 = parameters[1];
				double ra  = parameters[2];
				double re  = parameters[3];
				
				double pa1 = ra * pa0 /(1 - pa0 + ra * pa0);  //case 1.
				double pe1 = re * pe0 /(1 - pe0 + re * pe0);
				
				return groupLikehood(ctrl_coeff, pa0, pe0) + groupLikehood(case_coeff, pa1, pe1);
			}	
		};
		
//		BOBYQAOptimizer optimizer = new BOBYQAOptimizer(6);
//		PointValuePair results = optimizer.optimize(new MaxEval(maxEvalualtion),
//				new ObjectiveFunction(f),
//				new InitialGuess(getPAPE(ControlGeno)), //best guess allele as it was estimated from control data.
//				GoalType.MINIMIZE,
//				new SimpleBounds(new double[]{0,0,0,0}, new double[]{1,1,maxOddsRatio, maxOddsRatio})
//				);
//		
//		System.err.println(optimizer.getEvaluations());
		double[] bestGuess = new double[4];
		if (initial.length == 0) {
			System.arraycopy(getPAPE(ControlGeno), 0, bestGuess, 0, 2);
			bestGuess[2] = 1.0;
			bestGuess[3] = 1.0;
		}else {
//			System.err.println("initial:" + Arrays.toString(initial[0]));
			bestGuess = initial[0];
		}
		
		CMAESOptimizer optimizer = new CMAESOptimizer(maxEvalualtion, 0, true, 10, 0, new MersenneTwister(), false, null);
		PointValuePair results = optimizer.optimize(new MaxEval(maxEvalualtion),
				new ObjectiveFunction(f),
				new InitialGuess(bestGuess), //best guess allele as it was estimated from control data.
				GoalType.MINIMIZE,
				new SimpleBounds(new double[]{0,0,0,0}, new double[]{1,1,maxOddsRatio,maxOddsRatio})
				,new CMAESOptimizer.PopulationSize(5),
                new CMAESOptimizer.Sigma(new double[]{sigmaF,sigmaF,sigmaR,sigmaR})
				);
		
		double[] re = new double[4 + 1];
		re[0] = results.getValue();
		for (int i = 1; i < re.length; i++) {
			re[i] = results.getPointRef()[i-1];
		}
		return re;
	}
	
	/**
	 * Optimize pa, pe, and same r (odds ratio) for African and European.
	 * @param CaseGeno
	 * @param ControlGeno
	 * @param initial, set inital if has, otherwise optimize by program itself.
	 * @return logLikelyhood, pa0, pe0, r
	 */
	public static double[] sameOddsOptimizer(int[][] CaseGeno, int[][] ControlGeno, double[]... initial) {
		
		double[] case_coeff = coefficient(CaseGeno);
		double[] ctrl_coeff = coefficient(ControlGeno);
		
		MultivariateFunction f = new MultivariateFunction() {
			@Override
			/**
			 * All four parameters.
			 * pa0: ref. allele frequency in Africa population.
			 * pe0: ref. allele frequency in European population.
			 * ra:  Odds ratio in Africa population.
			 * re:  Odds ratio in European population. 
			 * parameters: pa0, pe0, r
			 */
			public double value(double[] parameters) {  // 3 parameters: pa0, pe0, r, same odds ratio for African and European.
				double pa0 = parameters[0];       //control 0.
				double pe0 = parameters[1];
				double ra  = parameters[2];
				double re  = parameters[2];
				
				double pa1 = ra * pa0 /(1 - pa0 + ra * pa0);  //case 1.
				double pe1 = re * pe0 /(1 - pe0 + re * pe0);
				
				return groupLikehood(ctrl_coeff, pa0, pe0) + groupLikehood(case_coeff, pa1, pe1);
			}	
		};
		
//		System.err.println("AF:" + Arrays.toString(getPAPE(ControlGeno)));
		double[] bestGuess = new double[3];
		if (initial.length == 0) {
			System.arraycopy(getPAPE(ControlGeno), 0, bestGuess, 0, 2);
			bestGuess[2] = 1.0;
		}else {
//			System.err.println("initial:" + Arrays.toString(initial[0]));
			bestGuess = initial[0];
		}
		
//		BOBYQAOptimizer optimizer = new BOBYQAOptimizer(6);
//		PointValuePair results = optimizer.optimize(new MaxEval(maxEvalualtion),
//				new ObjectiveFunction(f),
//				new InitialGuess(bestGuess), //best guess allele as it was estimated from control data.
//				GoalType.MINIMIZE,
//				new SimpleBounds(new double[]{0,0,0}, new double[]{1,1,maxOddsRatio})
//				);
		CMAESOptimizer optimizer = new CMAESOptimizer(maxEvalualtion, 0, true, 10, 0, new MersenneTwister(), false, null);
		PointValuePair results = optimizer.optimize(new MaxEval(maxEvalualtion),
				new ObjectiveFunction(f),
				new InitialGuess(bestGuess), //best guess allele as it was estimated from control data.
				GoalType.MINIMIZE,
				new SimpleBounds(new double[]{0,0,0}, new double[]{1,1,maxOddsRatio})
				,new CMAESOptimizer.PopulationSize(5),
                new CMAESOptimizer.Sigma(new double[]{sigmaF,sigmaF,sigmaR})
				);
//		System.err.println("sigmaF: " + sigmaF);
//		System.err.println("sigmaR:" + sigmaR);
		//System.err.println("sameOddsOptimizer:" + optimizer.getEvaluations());
//		System.err.println("OptimizedBest: " + results.getValue());
//		System.err.println("Initialed: " + optimizer.computeObjectiveValue(bestGuess));
		
		double[] re = new double[3 + 1];
		re[0] = results.getValue();
		for (int i = 1; i < re.length; i++) {
			re[i] = results.getPointRef()[i-1];
		}
		return re;
	}
	
	/**
	 * Very basic optimizer, only optimize allele frequency, set odds ratio constant as 1.
	 * @param CaseGeno
	 * @param ControlGeno
	 * @return logLikelyhood, pa0, pe0
	 */
	public static double[] onlyFrencyOptimizer(int[][] CaseGeno, int[][] ControlGeno) {
		
		double[] case_coeff = coefficient(CaseGeno);
		double[] ctrl_coeff = coefficient(ControlGeno);
		
		MultivariateFunction f = new MultivariateFunction() {
			@Override
			/**
			 * All four parameters.
			 * pa0: ref. allele frequency in Africa population.
			 * pe0: ref. allele frequency in European population.
			 * ra:  Odds ratio in Africa population.
			 * re:  Odds ratio in European population. 
			 * parameters: pa0, pe0, 1
			 */
			public double value(double[] parameters) {  // 2 parameters: pa0, pe0, r, set odds ratio as 1.
				double pa0 = parameters[0];       //control 0.
				double pe0 = parameters[1];
				double ra  = 1;
				double re  = 1;
				
				double pa1 = ra * pa0 /(1 - pa0 + ra * pa0);  //case 1.
				double pe1 = re * pe0 /(1 - pe0 + re * pe0);
				
				return groupLikehood(ctrl_coeff, pa0, pe0) + groupLikehood(case_coeff, pa1, pe1);
			}	
		};
		
//		System.err.println("AF:" + Arrays.toString(getPAPE(ControlGeno)));
//		BOBYQAOptimizer optimizer = new BOBYQAOptimizer(4);
//		PointValuePair results = optimizer.optimize(new MaxEval(maxEvalualtion),
//				new ObjectiveFunction(f),
//				new InitialGuess(getPAPE(ControlGeno)), //best guess allele as it was estimated from control data.
//				GoalType.MINIMIZE,
//				new SimpleBounds(new double[]{0,0}, new double[]{1,1})
//				);
		
		CMAESOptimizer optimizer = new CMAESOptimizer(30000, 0, true, 10, 0, new MersenneTwister(), false, null);
		PointValuePair results = optimizer.optimize(new MaxEval(maxEvalualtion),
				new ObjectiveFunction(f),
				new InitialGuess(getPAPE(ControlGeno)), //best guess allele as it was estimated from control data.
				GoalType.MINIMIZE,
				new SimpleBounds(new double[]{0,0}, new double[]{1,1})
				,new CMAESOptimizer.PopulationSize(5),
                new CMAESOptimizer.Sigma(new double[]{sigmaF,sigmaF})
				);
		
		//System.err.println("onlyFrencyOptimizer:" + optimizer.getEvaluations());
//		System.err.println("Frequency:" + Arrays.toString(getPAPE(ControlGeno)));
		
		double[] re = new double[2 + 1];
		re[0] = results.getValue();
		for (int i = 1; i < re.length; i++) {
			re[i] = results.getPointRef()[i-1];
		}
		return re;
	}
	
	/**
	 * Estimated allele frequency for pa and pe.
	 * @param genoCount
	 * @return
	 */
	private static double[] getPAPE(int[][] genoCount) {
//		for (int[] is : genoCount) {
//			System.err.println("genoCount: "+Arrays.toString(is));
//		}
		
		double AECount = genoCount[1][0] + genoCount[1][1] + genoCount[1][2]; // A, E admixture part.
		
		double pACount = 2*genoCount[0][0] + genoCount[0][1] + genoCount[1][0] + 0.5 * genoCount[1][1]; //ref. allele count.
		double pATotal = 2*(genoCount[0][0] + genoCount[0][1] + genoCount[0][2]) +  AECount;
		
		double pECount = 2*genoCount[2][0] + genoCount[2][1] + genoCount[1][0] + 0.5 * genoCount[1][1]; //ref. allele count.
		double pETotal = 2*(genoCount[2][0] + genoCount[2][1] + genoCount[2][2]) + AECount;
//		System.err.println("pECount: " + pECount);
//		System.err.println("pETotal: " + pETotal);
		return new double[]{pACount/pATotal, pECount/pETotal};
	}
	
	/**
	 * Get ref allele frequency of Pop1 and Pop2 from genotype configuration.
	 * @param genoCount
	 * @return
	 */
	private static double[] getHaplotypeP0P1(int[][] genoCount) {
		return new double[]{genoCount[0][0]*1.0/(genoCount[0][0] + genoCount[0][1]), genoCount[1][0] * 1.0/(genoCount[1][0] + genoCount[1][1])};
	}
	
	/**
	 * Compute likelihood for a group based on genotype count and allele frequency.
	 * @param coefficient, please ref. method: coefficient
	 * @param pa
	 * @param pe
	 * @return
	 * 
	 */
	private static double groupLikehood(double[] coefficient, double pa, double pe) {
//		return (2*gCount[0][0]+gCount[0][1])*FastMath.log(pa) + (2*gCount[0][2] + gCount[0][1])*FastMath.log(1-pa)
//				+ (gCount[1][0]+0.5*gCount[1][1])*FastMath.log(pa) + (gCount[1][2] + 0.5*gCount[1][1])*FastMath.log(1-pa)
//				+ (gCount[1][0]+0.5*gCount[1][1])*FastMath.log(pe) + (gCount[1][2] + 0.5*gCount[1][1])*FastMath.log(1-pe)
//				+ (2*gCount[2][0]+gCount[2][1])*FastMath.log(pe) + (2*gCount[2][2] + gCount[2][1])*FastMath.log(1-pe);
		
		double log_pa = FastMath.log(pa) ; double log_1_pa = FastMath.log(1-pa);
		double log_pe = FastMath.log(pe) ; double log_1_pe = FastMath.log(1-pe);
		return (	 coefficient[0] * log_pa + coefficient[1] * log_1_pa
			   + coefficient[2] * log_pa + coefficient[3] * log_1_pa
			   + coefficient[2] * log_pe + coefficient[3] * log_1_pe
			   + coefficient[4] * log_pe + coefficient[5] * log_1_pe) * -1.0;
	}
	
	/**
	 * Compute log likelihood for haplotype configurations.
	 * @param alleleCounts 
	 * @param pa
	 * @param pe
	 * @return
	 *  AlleleCounts
	 *  	ref		nonRef
	 *  A	
	 *  E
	 */
	private static double haplotypeGroupLikehood(int[][] alleleCounts, double pa, double pe) {
		double log_pa = FastMath.log(pa) ; double log_1_pa = FastMath.log(1-pa);
		double log_pe = FastMath.log(pe) ; double log_1_pe = FastMath.log(1-pe);
		return (	 
					alleleCounts[0][0] * log_pa + alleleCounts[0][1] * log_1_pa
					+ alleleCounts[1][0] * log_pe + alleleCounts[1][1] * log_1_pe 
			   ) * -1.0;
	}
	
	/**
	 * Convert genotype count to coefficient.
	 * @param gCount
	 * @return
	 * * * Genotype count format.
	 *   	2	1	0 Copy of ref. allele 
	 *    RR	RV	VV
	 * AA
	 * AE
	 * EE
	 */
	private static double[] coefficient(int[][] gCount){
		double[] results = new double[6];
		results[0] = 2*gCount[0][0] + gCount[0][1];
		results[1] = 2*gCount[0][2] + gCount[0][1];
		results[2] = gCount[1][0]   + 0.5*gCount[1][1];
		results[3] = gCount[1][2]   + 0.5*gCount[1][1];
		results[4] = 2*gCount[2][0] + gCount[2][1];
		results[5] = 2*gCount[2][2] + gCount[2][1];
		
//		System.err.println("Results: "+ Arrays.toString(results));
		return results;
	}
	
	
	/**
	 * Compute log likelihood and optimized parameters from haplotype counts.
	 * @param CaseGeno 2*2 table for allele counts with African and European ancestry.
	 * 			Ref  nonRef		
	 * 	 A(pop1)
	 *   E(pop2)
	 * @param ControlGeno
	 * @param initial	Initial gess values for pa, pe, odds ratio of Africa, odds ratio for European.
	 * @return logLikelihood, pRefPop1, pRefPop2, odds ratio for Pop1, odds ratio for pop2.
	 */
	public static double[] haplotypeTwoRatiosOptimizer(int[][] CaseGeno, int[][] ControlGeno, double[]... initial) {
		
		MultivariateFunction f = new MultivariateFunction() {
			@Override
			/**
			 * All four parameters.
			 * pa0: ref. allele frequency in Africa population.
			 * pe0: ref. allele frequency in European population.
			 * ra:  Odds ratio in Africa population.
			 * re:  Odds ratio in European population. 
			 * parameters: pa0, pe0, ra, re
			 */
			public double value(double[] parameters) {  // 4 parameters: pa0, pe0, ra, re
				double pa0 = parameters[0];       //control 0. a for pop1, e for pop2.
				double pe0 = parameters[1];
//				double ra  = parameters[2];
//				double re  = parameters[3];
				
//				double pa1 = ra * pa0 /(1 - pa0 + ra * pa0);  //case 1.
//				double pe1 = re * pe0 /(1 - pe0 + re * pe0);
				
				//direct optimize allele frequency in case, convert to odds ratio later.
				double pa1 = parameters[2];  //case 1, a for pop1, e for pop2.
				double pe1 = parameters[3];  
				
				return haplotypeGroupLikehood(ControlGeno, pa0, pe0) + haplotypeGroupLikehood(CaseGeno, pa1, pe1);
			}	
		};
		
		double[] bestGuess = new double[4];
		if (initial.length == 0) {
			System.arraycopy(getHaplotypeP0P1(ControlGeno), 0, bestGuess, 0, 2);
			System.arraycopy(getHaplotypeP0P1(CaseGeno), 0, bestGuess, 2, 2);
//			bestGuess[2] = 1.0;
//			bestGuess[3] = 1.0;
		}else {
			bestGuess = initial[0];
		}
		
//		CMAESOptimizer optimizer = new CMAESOptimizer(maxEvalualtion, 0, true, 10, 0, new MersenneTwister(), false, null);
//		PointValuePair results = optimizer.optimize(new MaxEval(maxEvalualtion),
//				new ObjectiveFunction(f),
//				new InitialGuess(bestGuess), //best guess allele as it was estimated from control data.
//				GoalType.MINIMIZE,
//				new SimpleBounds(new double[]{0,0,0,0}, new double[]{1,1,maxOddsRatio,maxOddsRatio})
//				,new CMAESOptimizer.PopulationSize(5),
//                new CMAESOptimizer.Sigma(new double[]{sigmaF,sigmaF,sigmaR,sigmaR})
//				);
		PointValuePair results = Optimizer.optimize(f, bestGuess);
		
		double[] re = new double[4 + 1]; //log likelihood + 4 estimated likelihood.
		re[0] = results.getValue();
		for (int i = 1; i < re.length; i++) {
			re[i] = results.getPointRef()[i-1];
		}
//		System.err.println(Arrays.toString(results.getPointRef()));
		//convert allele frequency to odds ratio.
		re[3] = ( re[3]/(1-re[3]) ) / (re[1]/(1-re[1])); //odds = case/control
		re[4] = ( re[4]/(1-re[4]) ) / (re[2]/(1-re[2]));
		
		return re;
	}
	
	/**
	 * 
	 * @param f
	 * @param bestGuess  allele frequencies.
	 * @return
	 */
	private static PointValuePair optimize(MultivariateFunction f, double[] bestGuess) {
//		System.err.println(Arrays.toString(bestGuess));
		double[] left = new double[bestGuess.length];
		double[] right = new double[bestGuess.length];
		double[] sigma = new double[bestGuess.length];
//		Arrays.fill(left, 0);
//		Arrays.fill(right, 1);
		Arrays.fill(left, freLeft);
		Arrays.fill(right, freRight);
		Arrays.fill(sigma, sigmaF);
		for (int i = 0; i < bestGuess.length; i++) {
			if (bestGuess[i] < freLeft) {
				bestGuess[i] = freLeft;
			}
			if (bestGuess[i] > freRight) {
				bestGuess[i] = freRight;
			}
		}
		
		CMAESOptimizer optimizer = new CMAESOptimizer(maxEvalualtion, 0, true, 10, 0, new MersenneTwister(1), false, null);
		PointValuePair results = optimizer.optimize(new MaxEval(maxEvalualtion),
				new ObjectiveFunction(f),
				new InitialGuess(bestGuess), //best guess allele as it was estimated from control data.
				GoalType.MINIMIZE,
				new SimpleBounds(left, right)
				,new CMAESOptimizer.PopulationSize(5),
                new CMAESOptimizer.Sigma(sigma)
				);
		
		return results;
	}
	
	/**
	 * Compute log likelihood and optimized parameters from haplotype counts.
	 * @param CaseGeno 2*2 table for allele counts with African and European ancestry.
	 * 			Ref  nonRef		
	 * 	 A(pop1)
	 *   E(pop2)
	 * @param ControlGeno
	 * @param initial	Initial gess values for pa, pe, odds ratio of Africa, odds ratio for European.
	 * @return logLikelihood, pRefPop1, pRefPop2, odds ratio (same odds ratio for Pop1 and pop2).
	 */
	public static double[] haplotypeOneRatiosOptimizer(int[][] CaseGeno, int[][] ControlGeno, double[]... initial) {
		
		MultivariateFunction f = new MultivariateFunction() {
			@Override
			/**
			 * All four parameters.
			 * pa0: ref. allele frequency in Africa population.
			 * pe0: ref. allele frequency in European population.
			 * ra:  Odds ratio in Africa population.
			 * re:  Odds ratio in European population. 
			 * parameters: pa0, pe0, ra, re
			 */
			public double value(double[] parameters) {  // 4 parameters: pa0, pe0, ra, re
				double pa0 = parameters[0];       //control 0.
				double pe0 = parameters[1];
				double ra  = parameters[2]; //same odds ratio for African and European.
				double re  = parameters[2];
				
				double pa1 = ra * pa0 /(1 - pa0 + ra * pa0);  //case 1.
				double pe1 = re * pe0 /(1 - pe0 + re * pe0);
				
				return haplotypeGroupLikehood(ControlGeno, pa0, pe0) + haplotypeGroupLikehood(CaseGeno, pa1, pe1);
			}	
		};
		
		double[] bestGuess = new double[3];
		if (initial.length == 0) {
			System.arraycopy(getHaplotypeP0P1(ControlGeno), 0, bestGuess, 0, 2);
			bestGuess[2] = 1.0;
		}else {
			bestGuess = initial[0];
		}
		
		CMAESOptimizer optimizer = new CMAESOptimizer(maxEvalualtion, 0, true, 10, 0, new MersenneTwister(), false, null);
		PointValuePair results = optimizer.optimize(new MaxEval(maxEvalualtion),
				new ObjectiveFunction(f),
				new InitialGuess(bestGuess), //best guess allele as it was estimated from control data.
				GoalType.MINIMIZE,
				new SimpleBounds(new double[]{0,0,0}, new double[]{1,1,maxOddsRatio})
				,new CMAESOptimizer.PopulationSize(5),
                new CMAESOptimizer.Sigma(new double[]{sigmaF,sigmaF,sigmaR})
				);
		
		double[] re = new double[3 + 1]; //log likelihood + 4 estimated likelihood.
		re[0] = results.getValue();
		for (int i = 1; i < re.length; i++) {
			re[i] = results.getPointRef()[i-1];
		}
		return re;
	}
	
	/**
	 * Compute log likelihood and optimized parameters from haplotype counts.
	 * @param CaseGeno 2*2 table for allele counts with African and European ancestry.
	 * 			Ref  nonRef		
	 * 	 A(pop1)
	 *   E(pop2)
	 * @param ControlGeno
	 * @param initial	Initial gess values for pa, pe, odds ratio of Africa, odds ratio for European.
	 * @return logLikelihood, pa, pe.
	 */
	public static double[] haplotypeOnlyFrequencyOptimizer(int[][] CaseGeno, int[][] ControlGeno) {
//		for (int[] is : ControlGeno) {
//			System.err.println("ControlGeno" + Arrays.toString(is));
//		}
//		for (int[] is : CaseGeno) {
//			System.err.println("CaseGeno" + Arrays.toString(is));
//		}
		
		MultivariateFunction f = new MultivariateFunction() {
			@Override
			/**
			 * All four parameters.
			 * pa0: ref. allele frequency in Africa population.
			 * pe0: ref. allele frequency in European population.
			 * ra:  Odds ratio in Africa population.
			 * re:  Odds ratio in European population. 
			 * parameters: pa0, pe0, ra, re
			 */
			public double value(double[] parameters) {  // 4 parameters: pa0, pe0, ra, re
				double pa0 = parameters[0];       //control 0.
				double pe0 = parameters[1];
//				double ra  = 1.0; //set odds ratio constant as 1. only optimize allele frequency.
//				double re  = 1.0;
				
//				double pa1 = ra * pa0 /(1 - pa0 + ra * pa0);  //case 1.
//				double pe1 = re * pe0 /(1 - pe0 + re * pe0);
				
				//odds ratio 1. same allele frequency in case and control
				double pa1 = pa0;
				double pe1 = pe0;
				
				return haplotypeGroupLikehood(ControlGeno, pa0, pe0) + haplotypeGroupLikehood(CaseGeno, pa1, pe1);
			}	
		};
		
		double[] bestGuess = new double[2];
		bestGuess = getHaplotypeP0P1(ControlGeno);

		
//		CMAESOptimizer optimizer = new CMAESOptimizer(maxEvalualtion, 0, true, 10, 0, new MersenneTwister(), false, null);
//		PointValuePair results = optimizer.optimize(new MaxEval(maxEvalualtion),
//				new ObjectiveFunction(f),
//				new InitialGuess(bestGuess), //best guess allele as it was estimated from control data.
//				GoalType.MINIMIZE,
//				new SimpleBounds(new double[]{0,0}, new double[]{1,1})
//				,new CMAESOptimizer.PopulationSize(5),
//                new CMAESOptimizer.Sigma(new double[]{sigmaF,sigmaF})
//				);
		
		PointValuePair results = Optimizer.optimize(f, bestGuess);
		
		double[] re = new double[bestGuess.length + 1]; //log likelihood + 2 estimated likelihood.
		re[0] = results.getValue();
		for (int i = 1; i < re.length; i++) {
			re[i] = results.getPointRef()[i-1];
		}
		return re;
	}
	
	
	
	/**
	 * Optimize ancestral omiga.
	 * @param ancestry [array, 0 or 1, 0 for pop1, 1 for pop2. each person two successive elements.]
	 * @param idAveProportion [array, each element: individual level of pop1 proportion].
	 * @return 2*IncreasedLogLikelihood, optimized omiga. *** Result from this optimizer is different with other optimizers. 
	 */
	public static double[] ancestryOptimizer(int[] ancestry, double[] idAveProportion) {
		
		//******CMAESOptimizer *****
		
//		MultivariateFunction f = new MultivariateFunction() {
//			@Override
//			/**
//			 * Only one parameter. omiga, in ADM model.
//			 */
//			public double value(double[] parameters) {  
//				double omiga = parameters[0];
//				return IntStream.range(0, idAveProportion.length)
//					.mapToDouble(i -> {
//						// number of pop1 ancestry. 0 for pop1, 1 for pop2. one individual two haps.
//						int copy = 2 - (ancestry[2*i] + ancestry[2*i +1]);
//						return LogAMDQik(idAveProportion[i], copy, omiga);
//					}).sum() * -1.0;
//			}
//		};
				
//		CMAESOptimizer optimizer = new CMAESOptimizer(maxEvalualtion, 0, true, 10, 0, new MersenneTwister(), false, null);
//		PointValuePair results = optimizer.optimize(new MaxEval(maxEvalualtion),
//				new ObjectiveFunction(f),
//				new InitialGuess(new double[]{1.0}), //best odds is 1.
//				GoalType.MINIMIZE,
//				new SimpleBounds(new double[]{0}, new double[]{maxOddsRatio})
//				,new CMAESOptimizer.PopulationSize(5),
//                new CMAESOptimizer.Sigma(new double[]{sigmaR})
//				);
		
		//****** BrentOptimizer *******
		UnivariateFunction f = new UnivariateFunction() {
			@Override
			public double value(double x) {
				double omiga = x;
				return IntStream.range(0, idAveProportion.length)
					.mapToDouble(i -> {
						// number of pop1 ancestry. 0 for pop1, 1 for pop2. one individual two haps.
						int copy = 2 - (ancestry[2*i] + ancestry[2*i +1]);
						return LogAMDQik(idAveProportion[i], copy, omiga);
					}).sum() * -1.0;
			}
		};
		
		BrentOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);
		UnivariatePointValuePair results = optimizer.optimize(
						 new MaxEval(maxEvalualtion)
						,new UnivariateObjectiveFunction(f)
						,GoalType.MINIMIZE
						,new SearchInterval(0, maxOmiga, 1.0) //low, high, initial starts.
						);
		
		double[] re = new double[1 + 1]; // delta_log_likelihood + estimated omiga. *** Result is different with other optimizer.
		double optimizedLog = results.getValue();
		double rawlog = f.value(1.0); //no optimize.
		re[0] = 2 * (rawlog - optimizedLog); //2*logdiff
		re[1] = results.getPoint();
		return re;
	}
	
	/**
	 * Compute log qi,k for AMD model. i_th individual, k copy of POP1 allele. 
	 * @param proption. Average proportion of POP1 allele for a individual.
	 * @param copy
	 * @return
	 */
	private static double LogAMDQik(double proption, int copy, double omiga) {
		double p1 = 1.0 - proption;
		double denominator = FastMath.pow(omiga*proption, 2) + 2*omiga*proption*p1 + FastMath.pow(p1, 2);
		switch (copy) {
		case 0:
			return FastMath.log(FastMath.pow(p1, 2) / denominator);
		case 1:
			return FastMath.log(2 * omiga * proption * p1 / denominator);
		default:
			return FastMath.log(FastMath.pow(omiga * proption, 2) / denominator);
		}
	}
}
