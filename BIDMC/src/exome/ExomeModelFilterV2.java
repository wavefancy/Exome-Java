package exome;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

/**
 * Filter variants according to inherent model.
 * 
 * @dominant
 * 1. All case members of a family share the alt. allele.
 *    (If any family member has an alt allele, all case members of this family should share this alt. allele, otherwise remove this as candidate.)
 * 2. Control individual can't have this alt. allele. Assume 100 penetrance.
 * 
 * @Version 2.0
 * 1. Remove the condition 2. in Version1.0. For the penetrance may not be 100%.
 * 
 * @version 2.1
 * 1. Output the genotype of the whole candidate family, not only the cases.
 * 
 * @version 2.2
 * 1. Add recessive model filter.
 *      Filter by recessive model, 100% penetrate. 
 *      1. all unaffected individual carrying non-refHomo genotype.
 *      2. In a candidate family, all affected individual carry refHomo genotype.
 * 
 * @version 2.3.
 * 1. add function for dominant model.
 * * Filter by dominant model, since 2.3.
         * Dominant (assume 100% penetrate).
         * 1. unaffected should be 00.
         * 2. In a candidate family affected should be 01 or 11.
 * @version 2.3.1
 * 1. separate individual by ';', other than ','
 * 
 * @version 2.4
 * 1. Add function to remove candidate sites for dominant model if ref-homo was 
 *    observed on chr1-chr22 or chrX for female.
 * 
 * @version 2.4.1
 * 1. Recognize 0|. or .|0 for dominant model, treat this as ref-homo.
 * 
 * @author wallace
 *
 */
public class ExomeModelFilterV2 {
	
//	private static int nCPUs = 0;
	private static boolean isData = false;
	private static Map<String, ArrayList<String>> caseFamilies = new HashMap<String, ArrayList<String>>(); //familyName->(idname1,idname2,...)
	private static Map<String, ArrayList<String>> controlFamilies = new HashMap<String, ArrayList<String>>();
	private static Map<String, ArrayList<String>> allFamilies = new HashMap<String, ArrayList<String>>();
	private static int ColStart = 0;
	private static String[] nameArr;
	private static Map<String, Integer> nameIndexMap = new HashMap<String, Integer>(); //idname -> Array_index
	private static Set<String> pedIdSet = new HashSet<>(); 
    private static List<String> ctrlNames = new LinkedList<String>(); //name list for all control.
    private static final DecimalFormat formater = new DecimalFormat(".####");
    
    //version 2.4
    static boolean checkDomAltHomo = true;
    //male 1, female 2.
    private static final Map<String, String> genderMap = new HashMap<>();
	
	private static String imodel = "dom"; //inherent model.
	
	public static void main(String[] args) {
		List<String> argList = new ArrayList<>(5);
		//check arguments.
		for (int i = 0; i < args.length; i++) {
			switch (args[i]) {
			case "-t":
				i++; System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", args[i]);
				break;			
			case "-mod":
				i++;
					switch (args[i]) {
					case "dom":
						imodel = args[i];
						break;
					case "rec":
						imodel = args[i];
						break;
					case "com":
						imodel = args[i]; break;
					default:
						System.err.println("The -mod parameter should be dom|rec|com for dominant, recessive or compound heterozygosity model.");
						System.exit(-1);
					}
				break;
			default:
				argList.add(args[i]);
			}
		}
		
		if(argList.size() != 1){
			help();
		}

//		ColStart = Integer.parseInt(argList.get(1)) -1; //shift to 0 based.
                ColStart = 9; //start column or vcf file, 0 based.
		new ExomeModelFilterV2().runApp(argList.get(0));
	}
	
	private static void help() {
		System.out.println("--------------------------------");
		System.out.println("    ExomeModelFilter    version: 2.4.1     Author:wavefancy@gmail.com");
		System.out.println("--------------------------------");
		System.out.println("Usages: \nparameter1: ped file."
//				+ "\nparameter2(int): Column index for individual seq. starts(Inclusive)."
				+ "\nparameter(-t  int, optional): number of cpus, default all available."
				+ "\nparameter(-mod String, optional): dom|rec|com for dominant, recessive or compound heterozygous model, default: dom."
				);
		System.out.println("Notes:"
				+ "\n1. Read vcf file from stdin and output to stdout."
				+ "\n3. Column index starts from 1." );
//				+ "\n4. Output PED file to stdout, MAP file to stderr.");
		System.out.println("--------------------------------");
		System.exit(-1);
	}
	
	/**
	 * Put one individual in corresponding family map.
	 * @param familyMap
	 * @param fName
	 * @param idName
	 */
	private void addToFamily(Map<String, ArrayList<String>> familyMap, String fName, String idName ) {
		if(! familyMap.containsKey(fName)){
			familyMap.put(fName, new ArrayList<>(5));
		}
		familyMap.get(fName).add(idName);
	}
	
	private void runApp(String pedFile) {
		try {
			Files.lines(Paths.get(pedFile))
				.map(String::trim)
				.filter(s-> !s.isEmpty())
				.forEach(s -> {
					String[] ss = s.split("\\s+");
					
					switch (ss[5]){
					case "2":
						pedIdSet.add(ss[1]);
						addToFamily(caseFamilies, ss[0], ss[1]);
						addToFamily(allFamilies, ss[0], ss[1]);
						break;
					case "1":
						pedIdSet.add(ss[1]);
						addToFamily(controlFamilies, ss[0], ss[1]);
						addToFamily(allFamilies, ss[0], ss[1]);
                                                ctrlNames.add(ss[1]);
                                                        
						break;
					default:
						System.err.println("Warnning: Skip no phenotype individual:" + ss[1]);
					}
                    
                    //recode gender info.
                    if (checkDomAltHomo) {
                        genderMap.put(ss[1], ss[4]); //male 1, female 2.
                    }
				});
			
			
			// Read vcf from stdin and output stdout.
			BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
			//BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(System.out));
			in.lines()
				.map(String::trim)
				.filter(s-> !s.isEmpty())
				.forEach(s-> processOne(s));
			in.close();
			
                        System.out.flush();
                        System.err.flush();
//			writer.flush();
//			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	/**
	 * Check if any body having alt. allele.
	 * @param names
	 * @param genotypes
	 * @return
	 */
//	private static boolean anyAltAllele(List<String> names, String[] genotypes) {
//		for (String string : names) {
//			String geno = genotypes[nameIndexMap.get(string)];
//			if(geno.charAt(0) != '.' && (geno.charAt(0) == '1' || geno.charAt(2) == '1'))
//				return true;
//		}
//		
//		return false; //all missing also return false.
//	}
	
	/**
	 * Check if all member having alt. allele,(genotype 11 or 01).
         * return false if all member are missing.
	 * @param names
	 * @param genotypes
	 * @return
	 */
	private static boolean allAltAlleleAllMissingFalse(List<String> names, String[] genotypes){
		String[] nonMissingN = names.stream()
				.filter(s -> {					
					return genotypes[nameIndexMap.get(s)].charAt(0) != '.';})
				.toArray(String[]::new);
		
		if (nonMissingN.length <= 0) { //all missing.
			return false;
		}
		
		//check non-missing individual, if somebody without alt. allele.
		for (String string : nonMissingN) {
			String geno = genotypes[nameIndexMap.get(string)]; //check 00 or 0|. or .|0.
			if((geno.charAt(0) == '0' || geno.charAt(0) == '.' ) && ( geno.charAt(2) == '0' || geno.charAt(2) == '.')) //no alt. allele.
				return false;
		}
		
		return true;
	}
	
	/**
	 * Count the number of nonmissing individuals.
	 * @param names
	 * @param genotypes
	 * @return
	 */
	private static long nonMissingCount(List<String> names, String[] genotypes) {
		return names.stream()
			.filter(s -> genotypes[nameIndexMap.get(s)].charAt(0) != '.')
			.count();
	}
        
        /**
         * Get nonmissing id names.
         * @param names
         * @param genotypes
         * @return 
         */
    private static String[] nonMissingNames(List<String> names, String[] genotypes) {
		return names.stream()
			.filter(s -> genotypes[nameIndexMap.get(s)].charAt(0) != '.')
			.toArray(String[]::new);
	}
	
	/**
	 * Check if all individuals in a family are missing.
	 * @param names
	 * @param genotypes
	 * @return
	 */
	private static boolean allMissing(List<String> names, String[] genotypes){
		return names.stream()
			.allMatch(s -> genotypes[nameIndexMap.get(s)].charAt(0) == '.');
	}
        
        /**
         * Check whether any individual from 'names' has refHomo[11] genotype.
         * @param names
         * @param genotypes
         * @return true if 
         */
    private static boolean hasRefHomo(List<String> names, String[] genotypes){
//            names.stream()
//              .forEach(s -> System.out.println(s + ":" + genotypes[nameIndexMap.get(s)]));
            
            boolean hasRefHomo = names.stream()        // anyMath return false if stream is empty.
               .filter(s -> genotypes[nameIndexMap.get(s)].charAt(0) != '.') //non-missing individuals.
               .anyMatch(s -> genotypes[nameIndexMap.get(s)].charAt(0) != '0'
                           && genotypes[nameIndexMap.get(s)].charAt(2) != '0'
               );
            return hasRefHomo;           
        }
        
        /**
         * Check whether all individuals carrying altHomo[11] genotype.
         * ** if all individuals are missing, return false.
         * @param names
         * @param genotypes
         * @return 
         */
        private static boolean allAltHomoAllMissingFalse(List<String> names, String[] genotypes){
            String[] nonMissing = nonMissingNames(names, genotypes);
            if(nonMissing.length == 0){
                return false;
            }
            
            long cc =  Arrays.stream(nonMissing)
              .filter(s-> genotypes[nameIndexMap.get(s)].charAt(0) != '0'
                         && genotypes[nameIndexMap.get(s)].charAt(2) != '0'
              ).count();
            
            if(cc == nonMissing.length){
                return true;
            }else{
                return false;
            }
        }
        
         /**
         * Check whether all individuals carrying refHomo genotype.
         * ** 1. if all individuals are missing, return true.
         * ** 2. if all individuals are refHomo (00) return true.
         * @param names
         * @param genotypes
         * @return 
         */
        private static boolean allRefHomoAllMissingTrue(List<String> names, String[] genotypes){
//            for(String s : names){
//                System.out.println(genotypes[nameIndexMap.get(s)]);
//            }
            
            String[] nonMissing = nonMissingNames(names, genotypes);
            if(nonMissing.length == 0){
                return true;
            }
            
            long cc =  Arrays.stream(nonMissing)
              .filter(s-> genotypes[nameIndexMap.get(s)].charAt(0) == '0'
                         && genotypes[nameIndexMap.get(s)].charAt(2) == '0'
              ).count();
            
            if(cc == nonMissing.length){
                return true;
            }else{
                return false;
            }
        }
        
        /**
         * Get an individual's genotype and coverage.
         * @param name
         * @param genotypes
         * @return 
         */
        private static String getGenoAndCov(String name, String[] genotypes){
            String[] ss = genotypes[nameIndexMap.get(name)].split(":");
            if(ss.length >= 2){
                return ss[0] +":" + ss[1];
            }else{
                return ss[0];
            }
           
        }
	
	/**
	 * Determine and output one line.
	 * @param line
	 * @param writer
	 */
	private static void processOne(String line) {
		
            if(isData == false){
                    if (line.startsWith("#CHROM")) {
                            isData = true;
                            nameArr = line.split("\\s+");
                            //map idname to index
                            boolean ee = false;
                            for (int j = ColStart; j < nameArr.length; j++) {
                                    if (nameIndexMap.keySet().contains(nameArr[j])) { //check duplicate.
                                            ee = true;
                                            System.err.println("Error: Duplicate individuals were detected in VCF head. Name: " + nameArr[j]);
                                    }else{
                                            nameIndexMap.put(nameArr[j], j);
                                    }

                                    if(! pedIdSet.contains(nameArr[j])){
                                            ee = true;
                                            System.err.println("Error: This individual wasn't decleared in ped file or missing phenotype: " + nameArr[j]);
                                    }
                            }
                            //Check whether all individuals in pedIDSet are indexed.
                            for (String nameInPed : pedIdSet) {
                                    if (! nameIndexMap.containsKey(nameInPed)) {
                                            ee = true;
                                            System.err.println("Error: Individual decleared in ped file but not in VCF file: " + nameInPed);
                                    }
                            }

                            if (ee) {
                                    System.exit(-1);
                            }

                            //output title
                            StringJoiner sJoiner = new StringJoiner("\t");
//					for (int i = 0; i < 8; i++) {
//						sJoiner.add(nameArr[i]);
//					}

                            sJoiner.add(nameArr[0]);
                            sJoiner.add(nameArr[1]);
                            sJoiner.add(nameArr[3]);
                            sJoiner.add(nameArr[4]);

                            sJoiner.add("#CandidateFamilies");
                            sJoiner.add("#AverageFamilySize");
                            sJoiner.add("#CandidatesGenotype");
                            System.out.println(sJoiner.toString());
                            //sJoiner.add("#TotalNonMissingFamily");
//					sJoiner.add("#TotalNonMissingControlFamily");
//					writer.write(sJoiner.toString());
//					writer.newLine();
                    }

            }else {
                    String[] ss = line.split("\\s+");

                    switch (imodel) {
                    case "dom":
                            dominantModel(ss);
                            break;

                    case "rec":
                            recessiveModel(ss);
                            //System.err.println("Recessive model has not been implemented!");
                            //System.exit(-1);
                            break;
                    case "com":
                            System.err.println("Compound heterozygous model has not been implemented!");
                            System.exit(-1);
                            break;
                    default:
                            break;
                    }

            }
		
	}
        
        /**
         * Filter by recessive model, 100% penetrate. 
         * 1. all unaffected individual carrying 00/01 [refHom or Hetero.] genotype.
         * 2. In a candidate family, all affected individual carry [11]refHomo genotype.
         * @param oneLineArr 
         */
        private static void recessiveModel(String[] oneLineArr){
//            System.err.println(ctrlNames);
//            System.err.println(hasRefHomo(ctrlNames, oneLineArr));
            if (hasRefHomo(ctrlNames, oneLineArr)) {
                //has refHomo[11] in control, assume 100% penetrate, skip this site.
                //return;
            }else{
                //check candiates families. all unaffectd individual at this site carring non-refHomo genotype [00/01].
                //candidate family genotype, and total of individal in these families, including unaffected.
                int total = 0;
                int cfCount = 0; //candiate family count.
                //List<String> candidateFamilies = new LinkedList<String>();
                
                StringJoiner cfGeno = new StringJoiner(","); //genotype and coverage info. for candidate family.
                
                for (String cfname :  caseFamilies.keySet()) {
                    if(allAltHomoAllMissingFalse(caseFamilies.get(cfname), oneLineArr)){ //candidate family.
                        
                        //candidateFamilies.add(cfname);
                        cfCount += 1;
                        total += caseFamilies.get(cfname).size();
                        //affected geno in a family.
                        StringBuilder sb = new StringBuilder();
                        sb.append("("+cfname+"[");
                        for(String n: caseFamilies.get(cfname)){
                            sb.append(n).append(":");
                            sb.append(getGenoAndCov(n, oneLineArr));
                            sb.append(";"); //separate individual by ';'
                        }
                        sb.setLength(sb.length()-1); //remove last ",".
                        sb.append("]");
                        //unaffected geno in a family.
                        if(controlFamilies.containsKey(cfname)){
                            total += controlFamilies.get(cfname).size();
                            for(String n: controlFamilies.get(cfname)){
                                sb.append(";"); ////separate individual by ';'
                                sb.append(n).append(":");
                                sb.append(getGenoAndCov(n, oneLineArr));
                            }
                        }
                        sb.append(")");
                        cfGeno.add(sb.toString());
                    }
                }
                
                //skip if no candidate family.
                if(cfCount <= 0){
                    return;
                }
               
                //output final results.
                //chr pos ref alt #candidateFamilies #averageSize #Genotypes
                StringJoiner sj = new StringJoiner("\t");
                sj.add(oneLineArr[0]);
                sj.add(oneLineArr[1]);
                sj.add(oneLineArr[3]);
                sj.add(oneLineArr[4]);
                
                sj.add(Integer.toString(cfCount));
                sj.add(formater.format(total * 1.0 / cfCount));
                
                sj.add(cfGeno.toString());
            
                //output final results.
                System.out.println(sj.toString());
            }
        }
        
        /**
         * Filter by dominant model, since 2.3.
         * Dominant (assume 100% penetrate).
         * 1. unaffected should be 00.
         * 2. In a candidate family affected should be 01 or 11.
         * @param oneLineArr 
         */
        private static void dominantModel(String[] oneLineArr){
//            System.out.println("exome.ExomeModelFilterV2.dominantModel()");
            String[] ss = oneLineArr;
//            System.out.println(ctrlNames);
//            System.out.println("allRefHomoAllMissingTrue(ctrlNames, ss)" + allRefHomoAllMissingTrue(ctrlNames, ss));
            if(allRefHomoAllMissingTrue(ctrlNames, ss) == false){
                // pass this variants, unmet conditon 1.0.
            }else{// all ref homo[00] or all missing.
                //check candiates families. all unaffectd individual at this site carring 00.
                //candidate family genotype, and total of individal in these families, including unaffected.
                int total = 0;   //total individuals in candidate family, including unaffected.
                int cfCount = 0; //candiate family count.
                //List<String> candidateFamilies = new LinkedList<String>();
                
                StringJoiner cfGeno = new StringJoiner(","); //genotype and coverage info. for candidate family.
                
                for (String cfname :  caseFamilies.keySet()) {
                    //Check if all member having alt. allele,(genotype 11 or 01).
                    if(allAltAlleleAllMissingFalse(caseFamilies.get(cfname), oneLineArr)
                            
                            && checkAltHomo4Dom(caseFamilies.get(cfname), oneLineArr)
                            ){ //candidate family.
                        
                        //candidateFamilies.add(cfname);
                        cfCount += 1;
                        total += caseFamilies.get(cfname).size();
                        //affected geno in a family.
                        StringBuilder sb = new StringBuilder();
                        sb.append("("+cfname+"[");
                        for(String n: caseFamilies.get(cfname)){
                            sb.append(n).append(":");
                            sb.append(getGenoAndCov(n, oneLineArr));
                            sb.append(";"); //separate individual by ';'
                        }
                        sb.setLength(sb.length()-1); //remove last ",".
                        sb.append("]");
                        //unaffected geno in a family.
                        if(controlFamilies.containsKey(cfname)){
                            total += controlFamilies.get(cfname).size();
                            for(String n: controlFamilies.get(cfname)){
                                sb.append(";"); //separate individual by ';'
                                sb.append(n).append(":");
                                sb.append(getGenoAndCov(n, oneLineArr));
                            }
                        }
                        sb.append(")");
                        cfGeno.add(sb.toString());
                    }
                }
                
                //skip if no candidate family.
                if(cfCount <= 0){
                    return;
                }
               
                //output final results.
                //chr pos ref alt #candidateFamilies #averageSize #Genotypes
                StringJoiner sj = new StringJoiner("\t");
                sj.add(oneLineArr[0]);
                sj.add(oneLineArr[1]);
                sj.add(oneLineArr[3]);
                sj.add(oneLineArr[4]);
                
                sj.add(Integer.toString(cfCount));
                sj.add(formater.format(total * 1.0 / cfCount));
                
                sj.add(cfGeno.toString());
            
                //output final results.
                System.out.println(sj.toString());
            }
            
            
        }
        
        /**
         * True if all non-altHomo. or altHomo on male X. Otherwise false.
         * @param names
         * @param oneLineArr
         * @return 
         */
        private static boolean checkAltHomo4Dom(List<String> caseNames, String[] oneLineArr){
            for (String name : caseNames) {
                int index = nameIndexMap.get(name);
                //skip missing and gender info missing individuals.
                if(oneLineArr[index].charAt(0) == '.' || genderMap.get(name).equalsIgnoreCase("0")){
                    continue;
                }
                
                if (oneLineArr[index].charAt(0) != '0' && oneLineArr[index].charAt(2) != '0') { //alt homo.
                    //only permit male, X for alt-homo.
                    if (! (oneLineArr[0].substring(oneLineArr[0].length()-1, oneLineArr[0].length()).equalsIgnoreCase("x")
                            && genderMap.get(name).equalsIgnoreCase("1"))
                            ) {
                        return  false;
                    }
                }
            }
            
            return true;
        }
                
}
