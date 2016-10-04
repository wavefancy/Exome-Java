/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bootstrapByFam;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import org.docopt.Docopt;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

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
    
    static final Map<String, int[]> geneCountMap = new HashMap<>(); // Record the number of times of #Case>=Ctrl, #Case<Ctrl for each gene. 
    
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

		String a = ((String) opts.get("-a")).toLowerCase(Locale.ENGLISH); //get algorithm names.
		if (modelSet.contains(a)) {
			model = a;
		}else {
			System.err.println("ERROR: please set proper parameter for -a.");
			System.exit(-1);
		}
        
        //prepare data for computing.
        aFamArray = new Family[aFamMap.size()];
        for (Map.Entry<String, Object> entry : aFamMap) {
            String key = entry.getKey();
            Object value = entry.getValue();
            
        }
        
	}
    
    private void permute(){
        
    }
}
