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
import java.util.Iterator;
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
 * 1. Add function to remove candidate sites for dominant model if alt-homo was 
 *    observed on chr1-chr22 or chrX for female.
 *  
 *    1. gender coding, 1 for male, 2 for female.
 *    2. only permit male or gender missing, X for alt-homo.
 * 
 * @version 2.4.1
 * 1. Recognize 0|. or .|0 for dominant model, treat this as ref-homo.
 * 
 * @version 2.4.2
 * 1. Add function to filter by maximum of candidate families.
 *      A true rare variants should not be shared by many families.
 *      *** Only implemented for Dominant model at this stage.
 * 
 * @version 2.5.
 * 1. Permit alt homo on male Y for dominant model.
 * 2. Add function for checking compound heterozyous model.
 * 
 * @version 2.5.1.
 * 1. Add function to remove MNP sites when perform COM HET checking.
 *  MNP: http://www.sequenceontology.org/miso/release_2.2/term/SO:0001013
 * 
 * @version 2.5.2
 * 1. Add function to estimate rare variant family sharing rate.
 *      Only use family with het. genotype. family with all ref-homo or at least 
 *      one alt-homo will be excluded.
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
    private static final DecimalFormat formater = new DecimalFormat("#.####");
    
    private static String geneAnnoFile = ""; // gene annotation file, for compound hetero. model.
    private static Map<String,List<String>> geneAnnoMap = new HashMap<>(); //geneName -> variant list.
    private static String[][] comDataMatrix = null;
    private static LinkedList<String[]> tempList = new LinkedList<>(); //temp data for checking compound hetero.
    
    //version 2.4
    static boolean checkDomAltHomo = true;
    //male 1, female 2.
    private static final Map<String, String> genderMap = new HashMap<>();
	
	private static String imodel = "dom"; //inherent model.
    private static int MaxCandiateFamilies = 5000;
    
    private static boolean removeMNP = true;
    private static int MNPlen = 3; //the distance of two sites to treated as MNP.
    //example MNPS
    // chr19-6693067-G-C and chr19-6693069-A-G, 
    // http://exac.broadinstitute.org/variant/19-6693069-A-G
    // chr15-62942342-G-A and chr15-62942343-C-A
    // http://exac.broadinstitute.org/variant/15-62942342-G-A

	
	public static void main(String[] args) {
		List<String> argList = new ArrayList<>(5);
		//check arguments.
		for (int i = 0; i < args.length; i++) {
			switch (args[i]) {
			case "-t":
				i++; System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", args[i]);
				break;	
            case "-m":
                i++;    MaxCandiateFamilies = Integer.parseInt(args[i]);
                break;
            case "-a":
                i++; geneAnnoFile = args[i];
                break;
            case "-cmnp":
                removeMNP = false;
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
                    case "eshare":
                        imodel = args[i]; break;
					default:
						System.err.println("The -mod parameter should be dom|rec|com|eshare for dominant, recessive or compound heterozygosity model.");
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
		System.out.println("    ExomeModelFilter    version: 2.5     Author:wavefancy@gmail.com");
		System.out.println("--------------------------------");
		System.out.println("Usages: \nparameter1: ped file."
//				+ "\nparameter2(int): Column index for individual seq. starts(Inclusive)."
				+ "\nparameter(-t  int, optional): number of cpus, default all available."
				+ "\nparameter(-mod String, optional): "
                + "\n          dom|rec|com for dominant, recessive or compound heterozygous model, default: dom."
                + "\n          eshare: estimate the rare variant sharing rate (currently version, only support auto-chromosome.)."
                + "\nparameter(-m  int, optional): maxmium number of candidate families, (Default 5000, only for dom model, <=)."
                + "\nparameter(-a file, optional, required for compound hetero.), gene annotation file."
                + "\nparameter(-cmnp, optional, only applicable for compound hetero.), close the function to remove MNP sites."
				);
		System.out.println("Notes:"
				+ "\n1. Read vcf file from stdin and output to stdout."
                + "\n2. Code for gender: 1 for male, 2 for female."
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
			
            //for compound heterozygous model.
            //Read gene annotation file. 
            if (imodel.equalsIgnoreCase("com")) { //for compound heterozyous model.
                try {
                    Files.lines(Paths.get(geneAnnoFile))
                        .map(String::trim)
                        .filter(s -> !s.isEmpty())
                        .forEach(s -> {
                            String[] ss = s.split("\\s+"); //chr pos ref alt geneName.
                            if (!geneAnnoMap.containsKey(ss[4])) {
                                geneAnnoMap.put(ss[4], new LinkedList<String>());
                            }
                            StringJoiner sj = new StringJoiner("-");
                            for (int i = 0; i <= 3; i++) {
                                 sj.add(ss[i]);
                            }
                            geneAnnoMap.get(ss[4]).add(sj.toString());
                        });
                } catch (Exception ioException) {
                    System.err.println("ERROR: Can not find gene annotation file[-a]: " + geneAnnoFile);
                    System.exit(-1);
                }
                
            }
            
			// Read vcf from stdin and output stdout.
			BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
			//BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(System.out));
			in.lines()
				.map(String::trim)
				.filter(s-> !s.isEmpty())
				.forEach(s-> processOne(s));
			in.close();
			
            
            if(imodel.equalsIgnoreCase("com")){
                compoundHeterozygousModel();
            }
            
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
     * Estimate the allele sharing rate.
     * Only use families with more than one non-missing member.
     * Remove family if any individual with alt-homo genotype. hard to estimate the rate.
     * 
     * @param names
     * @param genotypes
     * @return [number of individual with het genotype, total number of individuals in a family.]
     */
    private static long [] estimateShareRate(List<String> names, String[] genotypes){
        long[] re = {-1,-1};
        String[] nonMissingN = names.stream()
				.filter(s -> {					
					return genotypes[nameIndexMap.get(s)].charAt(0) != '.';})
				.toArray(String[]::new);
		
		if (nonMissingN.length <= 1) { //all missing or only 1 individual.
			return re;
		}
        
        //check any alt-homo
        boolean anyAltHomo = Arrays.stream(nonMissingN)
                  .anyMatch(s->{
                      return genotypes[nameIndexMap.get(s)].charAt(0) != '0' && genotypes[nameIndexMap.get(s)].charAt(2) != '0';
                  });
                    
        if(anyAltHomo){
            return re;
        }
        
        // number of individual with het genotype.
        re[0] = Arrays.stream(nonMissingN)
                  .filter(s -> {
                      return genotypes[nameIndexMap.get(s)].charAt(0) != '0' || genotypes[nameIndexMap.get(s)].charAt(2) != '0';
                  })
                  .count();
        
        re[1] = nonMissingN.length;
        return re;
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
                            if(imodel.equalsIgnoreCase("dom") || imodel.equalsIgnoreCase("rec")){
                                StringJoiner sJoiner = new StringJoiner("\t");
                                sJoiner.add(nameArr[0]);
                                sJoiner.add(nameArr[1]);
                                sJoiner.add(nameArr[3]);
                                sJoiner.add(nameArr[4]);

                                sJoiner.add("#CandidateFamilies");
                                sJoiner.add("#AverageFamilySize");
                                sJoiner.add("#CandidatesGenotype");
                                System.out.println(sJoiner.toString());
                            }else if (imodel.equalsIgnoreCase("com")) {
                           
                                StringJoiner sJoiner = new StringJoiner("\t");
                                sJoiner.add("GeneName");
                                sJoiner.add("FamilyName");
                                sJoiner.add("VariantID");
                                sJoiner.add("FamilyGenotype");
                                System.out.println(sJoiner.toString());
                            }else if (imodel.equalsIgnoreCase("eshare")) {
                                StringJoiner sJoiner = new StringJoiner("\t");
                                sJoiner.add(nameArr[0]);
                                sJoiner.add(nameArr[1]);
                                sJoiner.add(nameArr[3]);
                                sJoiner.add(nameArr[4]);
                                sJoiner.add("#Carrier");
                                sJoiner.add("#Total");
                                sJoiner.add("SharingRate");
                                System.out.println(sJoiner.toString());
                            }
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
                            tempList.add(ss);
                        
                            //System.err.println("Compound heterozygous model has not been implemented!");
                            //System.exit(-1);
                            break;
                    case "eshare":
                            eshare(ss);
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
                
                //skip if no candidate family, OR too many candidate families.
                if(cfCount <= 0 || cfCount > MaxCandiateFamilies){
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
         * Filter by compound heterozyous model, since version 2.5.
         */
        private static void compoundHeterozygousModel(){
            comDataMatrix = new String[tempList.size()][];
            Map<String,Integer> variantIndexMap = new HashMap<>();
            int t_index = 0;
            for (Iterator<String[]> iterator = tempList.iterator(); iterator.hasNext();) {
                String[] ss = iterator.next();
//                StringJoiner sj = new StringJoiner("-");
//                sj.add(ss[0]);
//                sj.add(ss[1]);
//                sj.add(ss[3]);
//                sj.add(ss[4]);
                variantIndexMap.put(getVariantkey(ss), t_index);
                comDataMatrix[t_index++] =ss;
            }
            tempList.clear();
            
            //System.out.println(geneAnnoMap.toString());
            
            //checking compound hetero model.
            for (Map.Entry<String, List<String>> entry : geneAnnoMap.entrySet()) { // Iterate by gene.
                String gene = entry.getKey();
                List<String> value = entry.getValue();
                int[] vIndex =  value.stream()
                           .filter(s -> variantIndexMap.get(s) != null) //varinat in annotation but not in vcf file.
                           .mapToInt(s->{
//                               if (variantIndexMap.get(s) == null) {
//                                   System.err.println("ERROR: Can not find this entry in vcf file, please check vcf and annotation file: " + s);
//                                   System.exit(-1);
//                                }
                               return variantIndexMap.get(s);})
                           .toArray();
                
                //System.out.println(gene);
                //System.out.println(Arrays.toString(vIndex));
                
                //Iterate by families.
                for (String caseFamily : caseFamilies.keySet()) {
//                    int[] idIndex = caseFamilies.get(caseFamily).stream()
//                                        .mapToInt(s -> nameIndexMap.get(s))
//                                        .toArray();
                    
                    //checking for candidate sites.
                    int[] passed = Arrays.stream(vIndex)
                                    //share alt allele, 11/01
                                    .filter(s -> allAltAlleleAllMissingFalse(caseFamilies.get(caseFamily), comDataMatrix[s]))
                                    //remove 11 genotype, except male x|y as 11.
                                    .filter(s -> checkAltHomo4Dom(caseFamilies.get(caseFamily), comDataMatrix[s]))
                                    .toArray();
                    //System.err.println("Fam: " +caseFamily);
                    //System.err.println(Arrays.toString(passed));
                    
                    //checking for MNP sites.
                    if (removeMNP && passed.length >= 2) {
                        Set<Integer> mnps = new HashSet<>();
                        long[] pos = Arrays.stream(passed)
                                .mapToLong(s -> Integer.parseInt(comDataMatrix[s][1]))
                                .toArray();
                        //checking for mnp
                        for (int i = 0; i < pos.length; i++) {
                            for (int j = i+1; j < pos.length; j++) {
                                if (Math.abs(pos[i] - pos[j]) <= MNPlen ) {
                                    mnps.add(i);
                                    mnps.add(j);
                                }
                            }
                        }
                        //update passed
                        int[] newpassed = new int[passed.length - mnps.size()];
                        int temp_index = 0;
                        for (int i = 0; i < passed.length; i++) {
                            if (!mnps.contains(i)) {
                                newpassed[temp_index++] = passed[i];
                            }
                        }
                        passed = newpassed;
                    }
                    
                    
                    StringJoiner sj = new StringJoiner("\t");
                    //candidate gene.
                    if (passed.length >= 2) {
                        sj.add(gene);
                        sj.add(caseFamily);
                        
                        //Iterate variant list.
                        for (int i : passed) {
                            StringJoiner out = new StringJoiner("\t");
                            out.add(sj.toString());
                            out.add(getVariantkey(comDataMatrix[i]));
                            
                            //genotype infor for this family at this sites.
                            StringBuilder sb = new StringBuilder();
                            sb.append("("+caseFamily+"[");
                            for(String n: caseFamilies.get(caseFamily)){
                                sb.append(n).append(":");
                                sb.append(getGenoAndCov(n, comDataMatrix[i]));
                                sb.append(";"); //separate individual by ';'
                            }
                            sb.setLength(sb.length()-1); //remove last ",".
                            sb.append("]");
                            //unaffected geno in a family.
                            if(controlFamilies.containsKey(caseFamily)){
                                for(String n: controlFamilies.get(caseFamily)){
                                    sb.append(";"); //separate individual by ';'
                                    sb.append(n).append(":");
                                    sb.append(getGenoAndCov(n, comDataMatrix[i]));
                                }
                            }
                            sb.append(")");
                            out.add(sb.toString());
                            System.out.println(out.toString());
                        }
                    }
                }
            }
        }
        
        /**
         * VCF variant key, chr-pos-ref-alt.
         * @param ss
         * @return 
         */
        private static String getVariantkey(String[] ss){
            StringJoiner sj = new StringJoiner("-");
                sj.add(ss[0]);
                sj.add(ss[1]);
                sj.add(ss[3]);
                sj.add(ss[4]);
            return sj.toString();
        }
        
        /**
         * True if all non-altHomo. or altHomo on male X or Y. Otherwise false.
         * @param names
         * @param oneLineArr
         * @return 
         */
        private static boolean checkAltHomo4Dom(List<String> caseNames, String[] oneLineArr){
            for (String name : caseNames) {
                int index = nameIndexMap.get(name);
                //skip missing
                //if(oneLineArr[index].charAt(0) == '.' || genderMap.get(name).equalsIgnoreCase("0")){
                if(oneLineArr[index].charAt(0) == '.'){
                    continue;
                }
                
                if (oneLineArr[index].charAt(0) != '0' && oneLineArr[index].charAt(2) != '0') { //alt homo.
                    //alt homo not on X|Y chromosome.  return false.
                    if(! (oneLineArr[0].substring(oneLineArr[0].length()-1, oneLineArr[0].length()).equalsIgnoreCase("x")
                            || oneLineArr[0].substring(oneLineArr[0].length()-1, oneLineArr[0].length()).equalsIgnoreCase("y")
                            )){
                        return false;
                    }
                    
                    //only permit male or gender missing, X|Y for alt-homo.
                    if (! ((oneLineArr[0].substring(oneLineArr[0].length()-1, oneLineArr[0].length()).equalsIgnoreCase("x")
                            || oneLineArr[0].substring(oneLineArr[0].length()-1, oneLineArr[0].length()).equalsIgnoreCase("y")
                            )
                            && (genderMap.get(name).equalsIgnoreCase("1") || genderMap.get(name).equalsIgnoreCase("0")))
                            ) {
                        return  false;
                    }
                }
            }
            
            return true;
        }
        
         private static void eshare(String[] oneLineArr){
//            System.out.println("exome.ExomeModelFilterV2.dominantModel()");
            String[] ss = oneLineArr;

//            if(allRefHomoAllMissingTrue(ctrlNames, ss) == false){ 
//                // pass this variants, unmet conditon 1.0.
//                // remove family with all missing, or no alt allele.
//            }else{
                //check candiates families. all unaffectd individual at this site carring 00.
                //candidate family genotype, and total of individal in these families, including unaffected.
                int total = 0;   //total individuals in candidate family, only affected.
                int hetCount = 0; //Number of indiviudal with het. genotype.
                //List<String> candidateFamilies = new LinkedList<String>();
                
                StringJoiner cfGeno = new StringJoiner(","); //genotype and coverage info. for candidate family.
                
                for (String cfname :  caseFamilies.keySet()) {
                    //all missing or no alt allele, do not need checking.
                    if(allRefHomoAllMissingTrue(caseFamilies.get(cfname), oneLineArr)){
                        continue;
                    }
                           
                    long [] re = estimateShareRate(caseFamilies.get(cfname), oneLineArr);
                    if (re[0] >= 0) {
                        hetCount += re[0];
                        total += re[1];
                    }
                }
             
                //chr pos ref alt #candidateFamilies #averageSize #Genotypes
                StringJoiner sj = new StringJoiner("\t");
                sj.add(oneLineArr[0]);
                sj.add(oneLineArr[1]);
                sj.add(oneLineArr[3]);
                sj.add(oneLineArr[4]);
                
                sj.add(Integer.toString(hetCount));
                sj.add(Integer.toString(total));
                
                if (total > 0) {
                    sj.add(formater.format(hetCount * 1.0 / total));
                }else{
                    sj.add("0");
                }
                
                sj.add(cfGeno.toString());
            
                //output final results.
                System.out.println(sj.toString());
//            }
        }
}
