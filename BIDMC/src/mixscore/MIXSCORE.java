package mixscore;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.util.MathArrays;
import org.docopt.Docopt;

/**
 * 
 * @author wavefancy@gmail.com
 * 
 * @version 3.0
 * 1. Change input as haps/sample format.
 * 
 * @version 3.0
 * 1. Use Brent's method to optimize univariate function.
 * 
 * @version 3.7.
 * 2. Add frequency bounds [0.00001, 0.99999], avoid optimize fail for log(fre).
 *
 */
public class MIXSCORE {
	private static String anctryFile;
	private static String phenoFile;
	private static String hapsFile;
	private static String idPropfile;
	private static String algorithm;
	private static int permuteNum = -1; //number of permutation.
	
	private static int[][] ancestryArr;
	private static int[][] genoArr; //**** Store the number of non-ref. allele. Different with input. *****
	private static int[][] hapsArr; // Store haplotypes from haps file.
	private static String[] idArr; //store snp id.
	private static double[] idProp; //the proportion of pop1 ancestry for each individual.
	private static int[] caseIndex;
	private static int[] ctrlIndex;
	private static final DecimalFormat formater = new DecimalFormat("#.####");
	
	
	private static final String doc =
		    "MIXSCORE.\n"
		    + "\n"
		    + "Usage:\n"
		    + "  mixscore -a algorithm -s sample -n ancestray [-p idProp] [-g haps] [--permute N] [-t NUMcpus] [--sf sigmaF] [--sr sigmaR]\n"
		    + "  mixscore (-h | --help)\n"
		    + "  mixscore --version\n"
		    + "\n"
		    + "Options:\n"
		    + "  -a algorithm  Algorithm name: SNP1|HET1|HET2|HAP1|HAP2.\n"
		    + "                SNP1: Association test for snp effects only, \n"
		    + "                      compare model with single odds ratio in two ancestries \n"
		    + "                      with odds ratio as 1, DF 1.\n"
		    + "                HET1: Association test for snp effects only, \n"
		    + "                      compare model with two different odds ratio to \n"
		    + "                      model with single odds ratio in two ancestries, DF 1.\n"
		    + "                HET2: Association test for snp effects only, \n"
		    + "                      compare model with two different odds ratio to \n"
		    + "                      model with odds ratio as 1 in two ancestries, DF 2.\n"
		    + "                HAP1: Association test based on [haplotype data].\n"
		    + "                      Compare model with single odds ratio to \n"
		    + "                      model with odds ratio as 1 in two ancestries, DF 1.\n"
		    + "                HAP2: Association test based on [haplotype data], \n"
		    + "                      compare model with two different odds ratio to \n"
		    + "                      model with odds ratio as 1 in two ancestries, DF 2.\n"
		    + "                ADM: Admixture association using cases only.\n"
		    + "				   ADMHAP2: Compute the chisquare sum of ADM and HAP2 test.\n"
		    + "  -s sample     Sample file, read binary phenotype from 'PHENO' column, 0 control, 1 case.\n"
		    + "  -g haps       Haps file.\n"
		    + "  -n ancestry   Ancestry file, 1/2: 1 for POP1, 2 for POP2.\n"
		    + "  -p idProp     Individual level POP1 proporition. One line one individual, order as sample file.\n"
		    + "  --permute N   Permute case and control N times, and output permutation results.\n"
		    + "  --sf sigmaF   Sigma parameter(step size) for allele frequency in ACM-ES, default 0.001, works well.\n"
		    + "  --sr sigmaR   Sigma parameter(step size) for odds ratio in ACM-ES, default 0.01, works well.\n"
		    + "  -t cpus       Number of cpus for computing.\n"
		    + "  -h --help     Show this screen.\n"
		    + "  --version     Show version.\n"
		    + "\n";
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Map<String, Object> opts =
		     new Docopt(doc).withVersion("MIXSCORE 3.7").parse(args);
//		     System.err.println(opts);
		if(opts.get("-t") != null){
			System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", (String) opts.get("-t"));
		}
		if(opts.get("--sf") != null){
			Optimizer.sigmaF = Double.parseDouble((String) opts.get("--sf"));
		}
		if(opts.get("--sr") != null){
			Optimizer.sigmaR = Double.parseDouble((String) opts.get("--sr"));
		}
		if (opts.get("-g") != null) {
			hapsFile = (String) opts.get("-g");
		}
		if (opts.get("-p") != null) {
			idPropfile = (String) opts.get("-p");
		}
		if (opts.get("--permute") != null) {
			permuteNum = Integer.parseInt((String) opts.get("--permute"));
		}
		anctryFile = (String) opts.get("-n");
		phenoFile = (String) opts.get("-s");
		Set<String> algorithmSet = new HashSet<String>();
		algorithmSet.add("SNP1");
		algorithmSet.add("HET2");
		algorithmSet.add("HAP1");
		algorithmSet.add("HAP2");
		algorithmSet.add("ADM");
		algorithmSet.add("ADMHAP2");
		String a = ((String) opts.get("-a")).toUpperCase(Locale.ENGLISH); //get algorithm names.
		if (algorithmSet.contains(a)) {
			algorithm = a;
		}else {
			System.err.println("ERROR: please set proper parameter for -a.");
			System.exit(-1);
		}
		
		Set<String> haplotypeAlgorithmSet = new HashSet<>();
		haplotypeAlgorithmSet.add("HAP1");
		haplotypeAlgorithmSet.add("HAP2");
		haplotypeAlgorithmSet.add("ADM");
		haplotypeAlgorithmSet.add("ADMHAP2");
		if(haplotypeAlgorithmSet.contains(a)){
			new MIXSCORE().runHaplotypeAlgorithms();
		}else{
//			new MIXSCORE().runAPP();
		}
		
		
	}
	
//	public void runAPP() {
//		try {
//			
//			BufferedReader reader = new BufferedReader(new FileReader(phenoFile));
//			String pheno = reader.readLine();
//			reader.close(); reader = null;
//			
//			caseIndex = IntStream.range(0, pheno.length())
//					.filter(s -> pheno.charAt(s) == '1')
//					.toArray();
//			ctrlIndex = IntStream.range(0, pheno.length())
//					.filter(s -> pheno.charAt(s) == '0')
//					.toArray();
//			
//			//read ancestry file.
//			final LinkedList<String> tList = new LinkedList<>();
//			Files.lines(Paths.get(anctryFile))
//				.forEach(s -> tList.add(s));
//			ancestryArr = new int[tList.size()][];
//			int index = 0;
//			for (String string : tList) {
//				ancestryArr[index++] = IntStream.range(0, string.length())
//									   .map(s -> Integer.valueOf(string.substring(s, s+1))) //0 for pop1, 1 pop2.
//									   .toArray();
//			}
//			tList.clear();
//			
//			//read haps file.
//			Files.lines(Paths.get(genoFile))
//				.forEach(s -> tList.add(s));
//			genoArr = new int[tList.size()][];
//			index = 0;
//			for (String string : tList) {
//				//Convert the number of ref. allele to the number of non-ref allele.
//				//for the convenience of count genotype frequency.
//				genoArr[index++] = IntStream.range(0, string.length())
//									.map(s ->{
//										if (string.charAt(s) == '2') {
//											return 0;
//										}else if (string.charAt(s) == '0') {
//											return 2;
//										}else {
//											return 1;
//										}
//									}).toArray();
//			}
//			tList.clear();
//			
////			for (int[] is : ancestryArr) {
////				System.err.println("ancestryArr:" + Arrays.toString(is));
////			}
////			for (int[] is : genoArr) {
////				System.err.println("genoArr:" + Arrays.toString(is));
////			}
//			
//			switch (algorithm) {
//			case "SNP1":
//				SNP1Test();
//				break;
//			case "HET2":
//				HET2Test();
//				break;
//				
//			default:
//				break;
//			}
//			
//			
//			
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		
//		
//	}
	
	/**
	 * Group of haplotype based algorithms.
	 */
	public void runHaplotypeAlgorithms() {
		try {
			//read from sample file.
			final LinkedList<String> tList = new LinkedList<>();
			Files.lines(Paths.get(phenoFile))
				.forEach(s -> tList.add(s));
			final int phenoIndex = Arrays.asList(tList.getFirst().toUpperCase(Locale.ENGLISH).split("\\s+")).indexOf("PHENO");
			if (phenoIndex < 0) {
				System.err.println("ERROR: Can't find 'PHENO' column from sample file(-s), please check!");
				System.exit(-1);
			}
			List<Character> phenos = new ArrayList<Character>(tList.size() -2);
			tList.stream()
				.skip(2) //skip two title lines.
				.forEach(s -> {
					char c = s.split("\\s+")[phenoIndex].charAt(0);
					if(c == '0' || c == '1'){
						phenos.add(c);
					}else {
						System.err.println("ERROR: Phenotype should be coded as 0/1 !");
						System.err.println("       Phenotype code ERROR at: " + s);
						System.exit(-1);
					}
				});
			caseIndex = IntStream.range(0, phenos.size())
					.filter(s -> phenos.get(s) == '1')
					.toArray();
			ctrlIndex = IntStream.range(0, phenos.size())
					.filter(s -> phenos.get(s) == '0')
					.toArray();
			
			//read from ancestry file.
			tList.clear();
			Files.lines(Paths.get(anctryFile))
				.forEach(s -> tList.add(s));
			ancestryArr = new int[tList.size()][];
			int index = 0;
//			for (String string : tList) {  //convert 1/2 to 0/1 presentation. 1 for POP1, 2 for POP2.
//				String[] sArr = string.split("\\s+");
//				ancestryArr[index++] = IntStream.range(0, sArr.length)  // 0/1: 0 for pop1, 1 for pop2.
//									   .map(s -> Integer.valueOf(sArr[s]) -1)
//									   .toArray();
//			}
			while(tList.size() >0){ //convert 1/2 to 0/1 presentation. 1 for POP1, 2 for POP2.
				String[] sArr = tList.pollFirst().split("\\s+");
				ancestryArr[index++] = IntStream.range(0, sArr.length)  // 0/1: 0 for pop1, 1 for pop2.
						   .map(s -> Integer.valueOf(sArr[s]) -1)
						   .toArray();
			}
			
			tList.clear();
			if (ancestryArr[0].length/2 != phenos.size()) {
				System.err.println("ERROR: The number of individuals in sample file is different with that in ancestry file.");
				System.exit(-1);
			}
			
			//read from haps file.
			if (hapsFile != null) {
				Files.lines(Paths.get(hapsFile))
					.forEach(s -> tList.add(s));
				hapsArr = new int[tList.size()][];
				idArr = new String[tList.size()];
				index = 0;
//				for (String string : tList) {
//					//ref 0, non-ref 1.
//					//for the convenience of count genotype frequency.
//					String[] sArr = string.split("\\s+");
//					idArr[index] = sArr[1] + "\t" + sArr[2];
//					hapsArr[index++] = IntStream.range(5, sArr.length)
//										.map(s -> Integer.valueOf(sArr[s]))
//										.toArray();
//				}
				while(tList.size() >0) {
					//ref 0, non-ref 1.
					//for the convenience of count genotype frequency.
					String[] sArr = tList.pollFirst().split("\\s+");
					idArr[index] = sArr[1] + "\t" + sArr[2];
					hapsArr[index++] = IntStream.range(5, sArr.length)
										.map(s -> Integer.valueOf(sArr[s]))
										.toArray();
				}
				tList.clear();
			}
			
			if (idPropfile != null) {
				tList.clear();
				Files.lines(Paths.get(idPropfile))
					.forEach(s -> tList.add(s));
				idProp = tList.stream()
						.mapToDouble(s -> Double.parseDouble(s))
						.toArray();
				if (idProp.length != phenos.size()) {
					System.err.println("ERROR: The number of individuals in sample file is different with that in idPropFile.");
					System.exit(-1);
				}
				tList.clear();
			}
			
			
//			for (int[] is : ancestryArr) {
//				System.err.println("ancestryArr:" + Arrays.toString(is));
//			}
//			for (int[] is : genoArr) {
//				System.err.println("genoArr:" + Arrays.toString(is));
//			}
			
			switch (algorithm) {
			case "HAP2":
				HAP2Test();
				break;
			case "ADM":
				ADMTest();
				break;
			case "ADMHAP2":
				if (permuteNum > 0) {
					ADMHAP2Permute();
				}else {
					ADMHAP2Test();
				}
				break;
				
			default:
				break;
			}
			
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	/**
	 * Compute genotype count according to ancesty label.
	 * @param row
	 * @param colArr
	 * @param ancestry
	 * @param genotype
	 * @return
	 *    RefRef.   Ref-nonRef  nonRef-nonRef.
	 * AA
	 * AE
	 * EE
	 */
	private int[][] genotypeCount(int row, int[] colArr, int[][] ancestry, int[][] genotype){
		int[][] results = new int[3][3];  //ancestry: 0 1 2, 0 1 2 copy of eur. ancestry
		for (int col : colArr) {          //genotype: 0 1 2, Copy of non-ref allele.
			results[ ancestry[row][col] ][genotype[row][col]] ++;
		}
		
//		System.err.println("colArr:" + Arrays.toString(colArr));
//		for (int[] is : results) {
//			System.err.println("results:"+ Arrays.toString(is));
//		}
		
		return results;
	}
	
	private void SNP1Test() {
		double[][] ratioR = new double[genoArr.length][];
		double[][] freR   = new double[genoArr.length][];
		
		IntStream.range(0, ratioR.length)
			.parallel()
			.forEach(row -> {
				int[][] caseGeno = genotypeCount(row, caseIndex, ancestryArr, genoArr);
				int[][] ctrlGeno = genotypeCount(row, ctrlIndex, ancestryArr, genoArr);
				
				freR[row] = Optimizer.onlyFrencyOptimizer(caseGeno, ctrlGeno);
				double[] initial = new double[3];
				System.arraycopy(freR[row], 1, initial, 0, 2);
				initial[2] = 1.0;
				ratioR[row] = Optimizer.sameOddsOptimizer(caseGeno, ctrlGeno, initial);
			});
		
		StringBuilder sBuilder = new StringBuilder();
		sBuilder.append("LogSameRatio\tpa\tpe\tr\tLogOnlyFre\tpa\tpe\tDF\tchisq");
		System.out.println(sBuilder.toString());
		for (int i = 0; i < ratioR.length; i++) {
			sBuilder.setLength(0);
			for (int c = 0; c < ratioR[0].length; c++) {
				sBuilder.append(formater.format(ratioR[i][c])).append("\t");
			}
			for (int c = 0; c < freR[0].length; c++) {
				sBuilder.append(formater.format(freR[i][c])).append("\t");
			}
			sBuilder.append("1").append("\t");
			sBuilder.append(formater.format(2* (freR[i][0] - ratioR[i][0])));
			System.out.println(sBuilder.toString());
		}
	}
	
	private void HET2Test() {
		double[][] ratioR = new double[genoArr.length][]; //two different odds ratio.
		double[][] freR   = new double[genoArr.length][]; //same odds ratio, contant as 1.
		
		IntStream.range(0, ratioR.length)
			.parallel()
			.forEach(row -> {
				int[][] caseGeno = genotypeCount(row, caseIndex, ancestryArr, genoArr);
				int[][] ctrlGeno = genotypeCount(row, ctrlIndex, ancestryArr, genoArr);
				
				freR[row] = Optimizer.onlyFrencyOptimizer(caseGeno, ctrlGeno);
				double[] initial = new double[4];  //four parameters for full optimizer.
				System.arraycopy(freR[row], 1, initial, 0, 2);
				initial[2] = 1.0; initial[3] = 1.0;
				ratioR[row] = Optimizer.fullOptimizer(caseGeno, ctrlGeno, initial);
			});
		
		StringBuilder sBuilder = new StringBuilder();
		sBuilder.append("LogSameRatio\tpa\tpe\tra\tre\tLogOnlyFre\tpa\tpe\tDF\tchisq");
		System.out.println(sBuilder.toString());
		for (int i = 0; i < ratioR.length; i++) {
			sBuilder.setLength(0);
			for (int c = 0; c < ratioR[0].length; c++) {
				sBuilder.append(formater.format(ratioR[i][c])).append("\t");
			}
			for (int c = 0; c < freR[0].length; c++) {
				sBuilder.append(formater.format(freR[i][c])).append("\t");
			}
			sBuilder.append("2").append("\t");
			sBuilder.append(formater.format(2* (freR[i][0] - ratioR[i][0])));
			System.out.println(sBuilder.toString());
		}
	}
	
	/**
	 * Haplotype based test, optimize for 2 odds ratios.
	 */
	private void HAP2Test() {
		double[][][] results = new double[hapsArr.length][2][];
		
	    //optimize and output results.
		IntStream.range(0, results.length)
			.parallel()
			.forEach(row -> {
				results[row] = hap2Optimizer(row, caseIndex, ctrlIndex);
			});
		
		StringBuilder sBuilder = new StringBuilder();
		sBuilder.append("SNPID\tPOS\tLogHap2Ratios\tpRefPop1Ctrl\tpRefPop2Ctrl\tr1\tr2\tLogOnlyFre\tpRefPop1\tpRefPop2\tDF\tchisq");
		System.out.println(sBuilder.toString());
		for (int i = 0; i < results.length; i++) {
			sBuilder.setLength(0);
			sBuilder.append(idArr[i]).append("\t");
			
			for (int c = 0; c < results[i][1].length; c++) {
				if (results[i][1][0] <= -1000 || Double.isNaN(results[i][1][0])) {
					sBuilder.append("NA").append("\t");
				}else {
					sBuilder.append(formater.format(results[i][1][c])).append("\t");
				}
//					sBuilder.append(formater.format(results[i][1][c])).append("\t");
//				sBuilder.append(results[i][1][c]).append("\t");
			}
			for (int c = 0; c < results[i][0].length; c++) {
				if (results[i][0][0] <= -1000 || Double.isNaN(results[i][0][0])) {
					sBuilder.append("NA").append("\t");
				}else {
					sBuilder.append(formater.format(results[i][0][c])).append("\t");
				}
//				sBuilder.append(formater.format(results[i][0][c])).append("\t");
			}
			if (results[i][0][0] <= -1000) {
				sBuilder.append("NA\tNA\tWARNING:Convergence_failed,_Probably_no_admixture_at_this_site_in_case/control_or_both!");
			}else if (Double.isNaN(results[i][0][0]) || Double.isNaN(results[i][1][0]) ) {
				sBuilder.append("NA\tNA\tWARNING:MAF_is_0_at_this_site_in_case/control_or_both!");
			}else {
				sBuilder.append("2").append("\t");
				sBuilder.append(formater.format(2* (results[i][0][0] - results[i][1][0])));
			}
			
			
			System.out.println(sBuilder.toString());
		}
	}
	
	/**
	 * Compute optimized results for hap2.
	 * @param row          SNP index.
	 * @param caseIDindex  Case individual index.
	 * @param ctrlIDindex  Control individual index.
	 * @return double[2][] 
	 *      logLikelihood, pRefPop1, pRefPop2. (constant odds ratio 1.) |
	 *      logLikelihood, pRefPop1, pRefPop2, oddsRatioPop1, oddsRatioPop2.
	 */
	private double[][] hap2Optimizer(int row, int[] caseIDindex, int[] ctrlIDindex) {
		int[] caseHapIndex = IDindex2HapIndex(caseIDindex);
		int[] ctrlHapIndex = IDindex2HapIndex(ctrlIDindex);
		int[][] caseAlleles = alleleCounts2(row, caseHapIndex, ancestryArr, hapsArr);
		int[][] ctrlAlleles = alleleCounts2(row, ctrlHapIndex, ancestryArr, hapsArr);
		
		double[][] results = new double[2][];
		
		try {
			results[0] = Optimizer.haplotypeOnlyFrequencyOptimizer(caseAlleles, ctrlAlleles);
			double[] initial = new double[4];
			System.arraycopy(results[0], 1, initial, 0, 2);
			System.arraycopy(results[0], 1, initial, 2, 2);
			
			results[1] = Optimizer.haplotypeTwoRatiosOptimizer(caseAlleles, ctrlAlleles, initial);
			return results;
			
		} catch (MaxCountExceededException e) {
			results[0] = new double[]{-1000,0,0}; //convergence failed, probably missing admixture at this sites.
			results[1] = new double[]{-1000,0,0,0,0};
		}
		return results;
	}
	
	/**
	 * Convert IDindex to haplotype index. Each individual two successive haplotypes.
	 * @param IDindex
	 * @return
	 */
	private int[] IDindex2HapIndex(int[] IDindex) {
		int[] results = new int[IDindex.length * 2];
		int index = 0;
		for (int i : IDindex) {
			results[index++] = 2*i;
			results[index++] = 2*i + 1;
		}
		return results;
	}
	
	/**
	 * Compute allele counts according to ancestry label.
	 * @param row
	 * @param colArr
	 * @param ancestry
	 * @param genotype
	 * @return
	 *    Ref.   nonRef.
	 * A
	 * E
	 */
	private int[][] alleleCounts(int row, int[] colArr, int[][] ancestry, int[][] genotype){
		int[][] results = new int[2][2];  //ancestry: 0 1, 0 for pop1, 1 for pop2
		for (int col : colArr) {          //genotype: 0 1, Copy of non-ref allele.
			results[ ancestry[row][col] ][genotype[row][col]] ++;
		}
		
//		System.err.println("colArr:" + Arrays.toString(colArr));
//		for (int[] is : results) {
//			System.err.println("results:"+ Arrays.toString(is));
//		}
		
		return results;
	}
	
	/**
	 * Compute allele counts according to ancestry label.
	 * @param row		snp index
	 * @param colArr    haplotype indexes to compute.
	 * @param ancestry  [0 for pop1, 1 for pop2]
	 * @param hapsArray [0 for ref, 1 for non-ref.]
	 * @return  
	 * 		ref		nonRef
	 * pop1
	 * pop2
	 */
	private int[][] alleleCounts2(int row, int[] colArr, int[][] ancestry, int[][] hapsArray) {
		int[][] results = new int[2][2];
		for (int c : colArr) {
			results[ancestry[row][c]][hapsArray[row][c]] ++;
		}
		return results;
	}
	
	/**
	 * Compute the chisquare sum of hap2 association and amd association.
	 */
	private void ADMHAP2Test() {
		//row -> results for (haplotype constant odds ratio 1. haplotype two odds ratio, ancestry test)
		double[][][] results = new double[hapsArr.length][3][];
		
	    //optimize and output results.
		IntStream.range(0, results.length)
			.parallel()
			.forEach(row -> {
				double[][] hap = hap2Optimizer(row, caseIndex, ctrlIndex);
				results[row][0] = hap[0];
				results[row][1] = hap[1];
				results[row][2] = amdOptimizer(row, caseIndex);
			});
		
		StringBuilder sBuilder = new StringBuilder();
		sBuilder.append("SNPID\tPOS\tLogHap2Ratios\tpRefPop1Ctrl\tpRefPop2Ctrl\tr1\tr2\tLogOnlyFre\tpRefPop1\tpRefPop2\tchisqHap2\tchisqAMD\tomiga\tDF\ttotalChisq");
		System.out.println(sBuilder.toString());
		for (int i = 0; i < results.length; i++) {
			sBuilder.setLength(0);
			sBuilder.append(idArr[i]).append("\t");
			
//			for (double[] ds : results[i]) {
//				System.out.println(Arrays.toString(ds));
//			}
			
			
			for (int c = 0; c < results[i][1].length; c++) {
				if (results[i][1][0] <= -1000 || Double.isNaN(results[i][1][0])) {
					sBuilder.append("NA").append("\t");
				}else {
					sBuilder.append(formater.format(results[i][1][c])).append("\t");
				}
//					sBuilder.append(formater.format(results[i][1][c])).append("\t");
//				sBuilder.append(results[i][1][c]).append("\t");
			}
			for (int c = 0; c < results[i][0].length; c++) {
				if (results[i][0][0] <= -1000 || Double.isNaN(results[i][0][0])) {
					sBuilder.append("NA").append("\t");
				}else {
					sBuilder.append(formater.format(results[i][0][c])).append("\t");
				}
//				sBuilder.append(formater.format(results[i][0][c])).append("\t");
			}
			
			if (results[i][0][0] <= -1000 || Double.isNaN(results[i][0][0]) || Double.isNaN(results[i][1][0]) ) {
				sBuilder.append("NA\t");
			}else {
//				sBuilder.append("3").append("\t");
				sBuilder.append(formater.format(2* (results[i][0][0] - results[i][1][0])));
				sBuilder.append("\t");
			}
			
			for (int c = 0; c < results[i][2].length; c++) {
				if (results[i][2][0] <= -1000 || Double.isNaN(results[i][2][0])) {
					sBuilder.append("NA").append("\t");
				}else {
					sBuilder.append(formater.format(results[i][2][c])).append("\t");
				}
			}
			
			if (results[i][0][0] <= -1000) {
				sBuilder.append("NA\tNA\tWARNING:Convergence_failed,_Probably_no_admixture_at_this_site_in_case/control_or_both!");
			}else if (Double.isNaN(results[i][0][0]) || Double.isNaN(results[i][1][0]) ) {
				sBuilder.append("NA\tNA\tWARNING:MAF_is_0_at_this_site_in_case/control_or_both!");
			}else {
				sBuilder.append("3").append("\t");
				sBuilder.append(formater.format(2* (results[i][0][0] - results[i][1][0]) + results[i][2][0]));
			}
			
			System.out.println(sBuilder.toString());
		}
	}
	
	private void ADMHAP2Permute() {
		String[] results = new String[hapsArr.length];
		//all individual ids.
		int[] ids = new int[caseIndex.length + ctrlIndex.length];
		System.arraycopy(caseIndex, 0, ids, 0, caseIndex.length);
		System.arraycopy(ctrlIndex, 0, ids, caseIndex.length, ctrlIndex.length);
		
	    //optimize and output results.
		IntStream.range(0, results.length)
//			.parallel()
			.forEach(row -> {
				PermuteResults pr = new PermuteResults();
				pr.setOriginal(getADMHAP2OnePermute(row, caseIndex, ctrlIndex));
				
				//start permutation.
				IntStream.range(0, permuteNum)
					.parallel()
					.anyMatch(i ->{
						//copy before shuffle, cause later in-place shuffle.
						int[] temp_ids = new int[ids.length];
						System.arraycopy(ids, 0, temp_ids, 0, ids.length);
						
						int[][] case_ctrl = shuffledCaseControl(temp_ids);
						double val = getADMHAP2OnePermute(row, case_ctrl[0], case_ctrl[1]);
						return ! pr.addAndCheckLarger(val);
					});
				results[row] = pr.toString();
			});
		
		System.out.println("Original2*LogDiff\tsmallerNum\tEqualorLargerNum");
		for (String string : results) {
			System.out.println(string);
		}
	}
	
	/**
	 * Shuffle the input id arrays.(*** in-place), and return shuffled case and control id index.
	 * @param allIdIndex
	 * @return int[2][] caseIDIndex, controlIDIndex array.
	 */
	private int[][] shuffledCaseControl(int[] allIdIndex) {
		MathArrays.shuffle(allIdIndex);
		int[][] re = new int[2][];
		re[0] = new int[caseIndex.length];
		re[1] = new int[ctrlIndex.length];
		
		System.arraycopy(allIdIndex, 0, re[0], 0, caseIndex.length);
		System.arraycopy(allIdIndex, caseIndex.length, re[1], 0, ctrlIndex.length);
		return re;
	}
	
	/**
	 * Get one permute results for ADMHAP2 test.
	 * @param row
	 * @param caseIDindex
	 * @param ctrlIDindex
	 * @return total 2*total_log_diff.
	 */
	private double getADMHAP2OnePermute(int row, int[] caseIDindex, int[] ctrlIDindex) {
		double[][] hap = hap2Optimizer(row, caseIDindex, ctrlIDindex);
		return (hap[0][0] - hap[1][0]) * 2 + amdOptimizer(row, caseIDindex)[0];
	}
	
	
	/**
	 * Optimize ancestry association omiga.
	 * @param row          row index.
	 * @param caseIDindex  case individual index.
	 * @return 2*IncreasedLogLikelihood, optimized omiga.
	 */
	private double[] amdOptimizer(int row, int[] caseIDindex) {
		int[] caseAncestryArr = new int[caseIDindex.length * 2];
		double[] caseAveProportion = new double[caseIDindex.length]; //average pop1 proportion.
		int index = 0; int index2 = 0;
		for (int i : caseIDindex) {
			caseAncestryArr[index++] = ancestryArr[row][2*i];
			caseAncestryArr[index++] = ancestryArr[row][2*i +1];
			caseAveProportion[index2++] = idProp[i];
		}
		
		return Optimizer.ancestryOptimizer(caseAncestryArr, caseAveProportion);
	}
	
	private void ADMTest() {
		double[][] results = new double[ancestryArr.length][];
		//case ancestry array.
		int[][] caseAncestry = new int[ancestryArr.length][];
		double[] caseIDAve = IntStream.range(0, caseIndex.length)
									.mapToDouble(i -> idProp[i])
									.toArray();
		
//		for (int[] is : ancestryArr) {
//			System.err.println("ancestryArr: " + Arrays.toString(is));
//		}
//		System.out.println("caseIndex: " + Arrays.toString(caseIndex));
		
		//case ancestry pop1 copys.
		IntStream.range(0, caseAncestry.length)
				.parallel()
				.forEach(row -> {
					int[] cArr = new int[2*caseIndex.length];
					int index = 0;
					for (int i : caseIndex) { //each individual two columns.
						cArr[index++] = ancestryArr[row][2*i];
						cArr[index++] = ancestryArr[row][2*i +1];
					}
					caseAncestry[row] = cArr;
				});
//		for (int[] is : caseAncestry) {
//			System.err.println("caseAncestry: " + Arrays.toString(is));
//		}
		
		
		IntStream.range(0, results.length)
			.parallel()
			.forEach(row ->{
				results[row] = Optimizer.ancestryOptimizer(caseAncestry[row], caseIDAve);
			});
		
		//output results.
		StringBuilder sBuilder = new StringBuilder();
		System.out.println("2LogDiff\tOmiga4Pop1\tDF");
		for (int i = 0; i < results.length; i++) {
			sBuilder.setLength(0);
//			sBuilder.append(idArr[i]);
			for (double dd : results[i]) {
				sBuilder.append(formater.format(dd)).append("\t");
			}
			sBuilder.append("1");
			System.out.println(sBuilder.toString());
		}
	}
}
