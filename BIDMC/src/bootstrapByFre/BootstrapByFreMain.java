/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bootstrapByFre;

import bootstrapByFam.BootstrapingMain;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.StringJoiner;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;
import org.docopt.Docopt;

/**
 * Bootstrapping from allele frequency.
 * @author wallace
 */
public class BootstrapByFreMain {
    static int nadaptIncrease = 1000; // each increasement of adpative permutation.
    static String imodel = ""; // inheritance model.
    static int ntotalPerm = (int) 1e5; // total number of permutation.
    static double pthreshold  = 0.05; // pvalue threshold for dropping adaptive permutation.
    static Map<String, List<Double>> geneFreMap = new HashMap<>();
    static final Map<String, double[]> geneFreMapArray = new HashMap<>();
    static final Map<String, Integer> geneCountMap = new  HashMap<>();
    static final DecimalFormat dformater = new DecimalFormat("0.00E00");
    static int ntotalIndividuals = 0;
    static double mutationRate = 1e-6; //mutation rate for denovo mutation.
    
    private static final String doc =
		    "BootstrapByFre.\n"
		    + "\n"
		    + "Usage:\n"
		    + "  BootstrapByFre -a afile -c cfile -m model -i nIDs [-n ptimes] [-t cpus] [-e etotal] [-r mrate] \n"
		    + "  BootstrapByFre (-h | --help)\n"
		    + "  BootstrapByFre --version\n"
		    + "\n"
		    + "Options:\n"
		    + "  -a afile      Input file for frequency annotation. line by line: geneName, siteFrequency, pos, MetaInfo \n"
		    + "  -c cfile      Input file for the number of candiate families for each gene.\n"
            + "                Line by line: geneName, candidateFamCount, MetaInfo.\n"
		    + "  -m model      Inheritance model for bootstrapping, dominant|recessive|compound.\n"
            + "  -i nIDs       Number of individuals.\n"
            + "  -n ptimes     Number of times each increase for addaptive permutation, default 1000.\n"
            + "  -e etotal     Total number of permutations, default 1e5.\n"
            + "  -r mrate      The mutation rate for de novo mutation, default 1e-6.\n"
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
		     new Docopt(doc).withVersion("BootstrapByFre 0.1").parse(args);
		     System.err.println(opts);
		if(opts.get("-t") != null){
			System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", (String) opts.get("-t"));
		}
        
        //ptime = Integer.parseInt((String) opts.get("-n"));
        ntotalIndividuals = Integer.parseInt((String) opts.get("-i"));
        if(opts.get("-e") != null){
			ntotalPerm = (int) Double.parseDouble((String) opts.get("-e"));
		}
        if(opts.get("-r") != null){
			mutationRate = Double.parseDouble((String) opts.get("-r"));
		}
              
        try {
            //Load data from annotation file.
            Files.lines(Paths.get((String) opts.get("-a")))
                    .sequential()
                    .map(s -> s.trim())
                    .filter(s -> !s.isEmpty())
                    .forEach((String s) -> {
                        String[] ss = s.split("\\s+", 3);
                        if (! geneFreMap.containsKey(ss[0])) {
                            geneFreMap.put(ss[0], new LinkedList<>());
                        }
                        try {
                            double fre = Double.parseDouble(ss[1]);
                            if (fre <= 0) {
                                fre = mutationRate;
                            }
                            geneFreMap.get(ss[0]).add(fre);
                        } catch (Exception e) {
                            System.err.println("SKIPPED: parse double(fre) error at: " + s);
                        }
                    });
            
            for (Map.Entry<String, List<Double>> e : geneFreMap.entrySet()) {
                double[] df = new double[e.getValue().size()];
                int index = 0;
                for (Iterator<Double> iterator = e.getValue().iterator(); iterator.hasNext();) {
                   df[index++] = iterator.next();
                }
                        
                geneFreMapArray.put(e.getKey(), df);
            }
            geneFreMap.clear();
            geneFreMap = null;
            
            //Load data from gene count file.
            Files.lines(Paths.get((String) opts.get("-c")))
                    .sequential()
                    .map(s -> s.trim())
                    .filter(s -> !s.isEmpty())
                    .forEach((String s) -> {
                        String[] ss = s.split("\\s+", 3);
                        int count = -1;
                        try {
                            count = Integer.parseInt(ss[1]);
                        } catch (Exception e) {
                            System.err.println("SKIPPED: parse int(count) error at: " + s);
                        }
                        if (count >= 0) {
                            if(! geneFreMapArray.containsKey(ss[0])){
                                System.err.println("ERROR: Can not find frequency annotation for gene: " + ss[0]);
                                System.exit(-1);
                            }
                            if (! geneCountMap.containsKey(ss[0])) {
                                geneCountMap.put(ss[0], count);
                            }else{
                                System.err.println("SKIPPED: Duplicated gene count entry: " + s);
                            }
                        }  
                    });
        } catch (IOException ex) {
            Logger.getLogger(BootstrapingMain.class.getName()).log(Level.SEVERE, null, ex);
        }
        
//        System.out.println(geneFreMap);
//        System.out.println(geneCountMap);
        
        //
        imodel = (String) opts.get("-m");
        imodel = imodel.toLowerCase(Locale.ENGLISH);
        System.out.println("GeneName\tCount\tPvalue\tTotalPerm");
        switch(imodel){
            case "compound":
                bootstrappingCompound();
                break;
            default:
                System.err.println("Can not recognize model parameter(-m): " + imodel);
        }
        
        System.out.flush();
        System.err.flush();
    }
    
    /**
     * bootstrapping according to compound heterozygosity model.
     */
    private static void bootstrappingCompound(){
        geneCountMap.entrySet().stream()
//            .parallel()
            .forEach(map -> {
                String gName = map.getKey();
                int count = map.getValue();
                StringJoiner sj = new  StringJoiner("\t");
                if (count == 0) {
                    sj.add(gName)
                      .add(Integer.toString(count))
                      .add(dformater.format(1.0))
                      .add("0");
                }else{ //adaptive bootstrapping starts from here.
                
                    int totalPerm = 0;
                    int totalSuccess = 0;
                    double pvalue = 0;
                    for (int i = 0; i < ntotalPerm; i += nadaptIncrease) {
                        totalPerm += nadaptIncrease;
                        totalSuccess += batchPvalueCompound(gName, count);
                        pvalue = totalSuccess*1.0/totalPerm;
                        if ( pvalue >= pthreshold) {
                            break;
                        }
                    }
                    
                    sj.add(gName)
                      .add(Integer.toString(count))
                      .add(dformater.format(pvalue))
                      .add(Integer.toString(totalPerm));
                }
                
                //output results:
                System.out.println(sj.toString());
            
            });
                
    }
    
    /**
     * Compute pvalue for one batch of bootstrapping.
     * @return 
     */
    private static long batchPvalueCompound(String geneName, int count){
        return  IntStream.range(0, nadaptIncrease)
                .parallel()
                .mapToLong(a -> { //bootstapping once.
                    
                    long cc = IntStream.range(0, ntotalIndividuals)
                        .mapToLong(i -> {
                        //generate psudo individual and check compound heterozygosity.
                        //generate genotype for individual i.
                        //System.out.println("ID:"+i);
                        Random random = ThreadLocalRandom.current();
//                        long hetSize = geneFreMap.get(geneName).stream()
                          long hetSize = Arrays.stream(geneFreMapArray.get(geneName))
                                .sequential()
                                .mapToInt(fre -> {
                                    int a1 = (random.nextDouble() <= fre) ? 1: 0;
                                    int a2 = (random.nextDouble() <= fre) ? 1: 0;

                                    return a1 + a2;
                                })
                                .filter(g -> g == 1)
                                .count();
                        
                        //System.err.println(geneName + "  " + hetSize);
                        return  hetSize;
                        })
                        .filter(c -> c >=2) //individual meet condition.
                        .count();
                    
                    //System.err.println(geneName + "  CP ID COUNT:" + cc);
                    return  cc; //number of sampled individuals met compound conditions.
                })
                .filter(c -> c >= count) 
                .count(); //count the number of extreme conditions.
    }
    
}
