/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bootstrapByFam;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import org.docopt.Docopt;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;

/**
 *
 * @author wallace
 */
public class BootstrapingMain {
    private static String model = "";
    static Map<String, Family> aFamMap = new HashMap<>(); // for affected famlies.
    static Family[] aFamArray = null;
    static Map<String, Family> uFamMap = new HashMap<>(); // for unaffected families.
    static Family[] uFamArray = null;
    
    static int ptime = 0; // number of permutations.
    
    static final Map<String, int[]> geneCountMap = new ConcurrentHashMap<>(); // Record the number of times of #Case>=Ctrl, #Case<Ctrl for each gene. 
    
	private static final String doc =
		    "BootstrappingRare.\n"
		    + "\n"
		    + "Usage:\n"
		    + "  BootstrappingRare -a afile -u ufile -m model -n ptimes \n"
		    + "  BootstrappingRare (-h | --help)\n"
		    + "  BootstrappingRare --version\n"
		    + "\n"
		    + "Options:\n"
		    + "  -a afile      Input file for affecteds. line by line: famName, geneName, MetaInfo \n"
		    + "  -u ufile      Input file for unaffecteds.\n"
		    + "  -m model      Inheritance model for bootstrapping, dominant|recessive|compound.\n"
            + "  -n ptimes     Number of times for permutations.\n"
		    + "  -t cpus       Number of cpus for computing.\n"
		    + "  -h --help     Show this screen.\n"
		    + "  --version     Show version.\n"
		    + "\n";
	
    /**
     * Read parameters and load data.
     * @param args 
     */
	public static void main(String[] args) {
        
		// TODO Auto-generated method stub
		Map<String, Object> opts =
		     new Docopt(doc).withVersion("BootstrappingRare 0.1").parse(args);
//		     System.err.println(opts);
		if(opts.get("-t") != null){
			System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", (String) opts.get("-t"));
		}
        
        ptime = Integer.parseInt((String) opts.get("-n"));
                 
        //Load data from file.
        try {
            Files.lines(Paths.get((String) opts.get("-a")))
                    .sequential()
                    .map(s -> s.trim())
                    .filter(s -> !s.isEmpty())
                    .forEach((String s) -> {
                        String[] ss = s.split("\\s+", 2);
                        if (! aFamMap.containsKey(ss[0])) {
                            aFamMap.put(ss[0], new Family(ss[0]));
                        }
                        //aFamMap.get(ss[0]).addGene(ss[1], ss[2]);
                        aFamMap.get(ss[0]).addGene(ss[1]);
                        if (!geneCountMap.containsKey(ss[1])) {
                            geneCountMap.put(ss[1], new int[2]);
                        }
                    });
            
            Files.lines(Paths.get((String) opts.get("-u")))
                    .sequential()
                    .map(s -> s.trim())
                    .filter(s -> !s.isEmpty())
                    .forEach((String s) -> {
                        String[] ss = s.split("\\s+", 2);
                        if (! uFamMap.containsKey(ss[0])) {
                            uFamMap.put(ss[0], new Family(ss[0]));
                        }
                        //uFamMap.get(ss[0]).addGene(ss[1], ss[2]);
                        uFamMap.get(ss[0]).addGene(ss[1]);
//                        if (!geneCountMap.containsKey(ss[1])) {
//                            geneCountMap.put(ss[1], new int[2]);
//                        }
                    });
        } catch (IOException ex) {
            Logger.getLogger(BootstrapingMain.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        HashSet<String> modelSet = new  HashSet<>();
        
		modelSet.add("dominant");
		modelSet.add("recessive");
		modelSet.add("compound");

		String a = ((String) opts.get("-m")).toLowerCase(Locale.ENGLISH); //get algorithm names.
		if (modelSet.contains(a)) {
			model = a;
		}else {
			System.err.println("ERROR: please set proper parameter for -m.");
			System.exit(-1);
		}
        
        //prepare data for computing.
        aFamArray = new Family[aFamMap.size()];
        int t_index  = 0;
        for (Map.Entry<String, Family> entry : aFamMap.entrySet()) {
            entry.getValue().listToArray();
            aFamArray[t_index++] = entry.getValue();
        }
        t_index = 0;
        uFamArray = new Family[uFamMap.size()];
        for (Map.Entry<String, Family> entry : uFamMap.entrySet()) {
            entry.getValue().listToArray();
            uFamArray[t_index++] = entry.getValue();
        }
        aFamMap = null;
        uFamMap = null;
        
        // permutation starts from here.
        IntStream.range(0, ptime)
                .parallel()
                .forEach(s -> {
                    permute();
                });
        
        // output results.
        DecimalFormat format = new DecimalFormat();
        System.out.println("Gene\tCaseBigger\tCtrolBiger\tPvalue");
        for (Map.Entry<String, int[]> entry : geneCountMap.entrySet()) {
            String key = entry.getKey();
            int[] value = entry.getValue();
            
            if ( value[0] == 0 && value[0] == value[1]) {
                System.out.println(key + "\t0\t0\tNA");
            }else{
                System.out.println(key + "\t" + Integer.toString(value[0]) +  
                        "\t" + Integer.toString(value[1]) + "\t" + format.format(value[1] * 1.0 / (value[0] + value[1])));
            }
        }
	}
    
    private static void geneCounter(Map<String,int[]> gMap, List<String> caseGeneArr, List<String> ctrlGeneArr){
        caseGeneArr.stream()
                .forEach((String s) -> {
                    gMap.get(s)[0] = gMap.get(s)[0] +1;
                });
        //Only count gene in candidate gene list.
        ctrlGeneArr.stream()
                .filter((String s) -> gMap.containsKey(s))
                .forEach((String s) -> {
                    gMap.get(s)[1] = gMap.get(s)[1] +1;
                });
//        System.out.println("------bootstrapByFam.BootstrapingMain.geneCounter()----");
//        System.out.println(caseGeneArr.toString());
//        System.out.println(ctrlGeneArr.toString());
    }
    
    private static void permute(){
        Random rndRandom = ThreadLocalRandom.current();
        
        Map<String,int[]> gMap = new HashMap<>(geneCountMap.size());  //geneName -> #AffectedFamily, #unaffectedFamily.
        for (Map.Entry<String, int[]> entry : geneCountMap.entrySet()) {
            String key = entry.getKey();
            int[] value = {0,0};
            gMap.put(key, value);
        }
        
        for (Family afam : aFamArray) {
            Family ufam = uFamArray[rndRandom.nextInt(uFamArray.length)];
            String[] aArr = afam.getGeneArrayCopy();
            String[] uArr = ufam.getGeneArrayCopy();
            List<String> aList = Arrays.asList(aArr);
            List<String> uList = Arrays.asList(uArr);
            if (aArr.length < uArr.length) {
                uList = FisherYatesArrayShuffle.Shuffle(uArr, aArr.length);
            }else if (aArr.length > uArr.length) {
                aList = FisherYatesArrayShuffle.Shuffle(aArr, uArr.length);
            }
            geneCounter(gMap, aList, uList);
        }
        
        
        //get summary data from permutation.
        for (Map.Entry<String, int[]> entry : gMap.entrySet()) {
            String key = entry.getKey();
            int[] value = entry.getValue();
            
            if (value[0] == 0 && value[0] == value[1]) {
                continue;
            }
            
            if (value[0] >= value[1]) {
                geneCountMap.get(key)[0] += 1;
            }else{
                geneCountMap.get(key)[1] += 1;
            }
        }
        
//        for (Map.Entry<String, int[]> entry : gMap.entrySet()) {
//            String key = entry.getKey();
//            int[] value = entry.getValue();
//            System.out.println(key + "\t" + Arrays.toString(value));
//        }
    }
}
