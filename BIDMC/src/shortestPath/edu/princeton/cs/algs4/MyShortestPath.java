package shortestPath.edu.princeton.cs.algs4;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.stream.IntStream;
import org.docopt.Docopt;

/**
 * 
 * @author wavefancy@gmail.com
 * 
 * @version 1.0
 * Compute shortest path based on princeton's code base.
 * Compute the shortest path between nodes.
 *
 */
public class MyShortestPath {
    
    private static String[] topGenes = null;   //top gene list.
    private static String[] knownGenes = null; //known gene list.
    private static int bootTimes = -1; //bootstrap times.
    private static int RandomPickGenes = -1; //random pick this number of genes.
    private final static Map<String,Integer> geneNameMap = new HashMap(); //map genename -> int number.
    private static int  tempIndex = 0;
	
    private static final String doc =
                "MyUndirectedWeightedShortestPath.\n"
                + "\n"
                + "Usage:\n"
                + "  MyUndirectedWeightedShortestPath -k knownGene (-i topGenes | -r int -b int) [-t cpus] \n"
                + "  MyUndirectedWeightedShortestPath (-h | --help)\n"
                + "  MyUndirectedWeightedShortestPath --version\n"
                + "\n"
                + "---------------------------\n"
                + "Read network from stdin, each line three columns: GeneA GeneB weight.\n"
                + "If only two columns, all edges were assigned the equal weight of 1.\n"
                + "Compute the avarage shorest path from topGenes with knownGenes. \n"
                + "mean(min_i(KnownGenes), i iterate all topGenes.\n"
                + "min_i: the minimal distance between topGene_i with all knownGenes."
                + "---------------------------\n"
                + "Options:\n"
                + "  -k knownGene  Input known gene list, one line one gene.\n"
                + "  -i topGenes   Input top gene list, one line one gene.\n"
                + "  -r int        Number of random picked genes for bootstrapping.\n"
                + "  -b int        Number of bootstrappings.\n"
                + "  -t cpus       Number of cpus for computing.\n"
                + "  -h --help     Show this screen.\n"
                + "  --version     Show version.\n"
                + "\n";
    
    /**
     * Get the average min distance between testGeneIDs with targetGeneIDs.
     * Average(min(test_i_ShortestDistanceToAllTargetGenes), iterate i))
     * @param testGeneIDs
     * @param targetGeneIDs
     * @param G
     * @return 
     */
    private static OptionalDouble averageMindistance(int[] testGeneIDs, int[] targetGeneIDs, EdgeWeightedGraph G){

        return Arrays.stream(testGeneIDs)
                        .parallel()
                        .mapToDouble(s->{
                            //min(distance to all knwon genes.)
                            DijkstraUndirectedSP sp = new DijkstraUndirectedSP(G, s);
                            double[] minimalDis = Arrays.stream(targetGeneIDs)
                                    .filter(k->{return sp.hasPathTo(k);})
                                    .mapToDouble(k->{
                                        return sp.distTo(k);
                                    })
                                    .toArray();
                            if (minimalDis.length >0) {
                                return Arrays.stream(minimalDis).min().getAsDouble();
                            }else{
                                return -1;
                            }
                        })
                        .filter(s-> s>=0 )
                        .average();
    }

    public static void main(String[] args) {
            Map<String, Object> opts =
                 new Docopt(doc).withVersion("MIXSCORE 3.7").parse(args);
//		     System.err.println(opts);
            if(opts.get("-t") != null){
                    System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", (String) opts.get("-t"));
            }
            //read knwon genes.
            if(opts.get("-k") != null){
                try {
                   knownGenes = Files.lines(Paths.get((String) opts.get("-k")))
                        .filter(s -> s.length() > 0)
                        .toArray(String[]::new);
                   
                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
            }
            //read top genes
            if(opts.get("-i") != null){
                try {
                   topGenes = Files.lines(Paths.get((String) opts.get("-i")))
                        .filter(s -> s.length() > 0)
                        .toArray(String[]::new);
                   
                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
            }
            if(opts.get("-r") != null){
                RandomPickGenes = Integer.parseInt((String) opts.get("-r"));
            }
            if(opts.get("-b") != null){
                bootTimes = Integer.parseInt((String) opts.get("-b"));
            }
            
//            System.out.println(Arrays.toString(topGenes));
//            System.out.println(Arrays.toString(knownGenes));
            //read network from stdin.
            BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
            final LinkedList<String[]> tempEdgesList = new LinkedList();
            in.lines()
              .filter(s -> s.length() > 0)
              .forEach(s ->{
                  String[] ss = s.split("\\s+");
                  if (!geneNameMap.containsKey(ss[0])) {
                      geneNameMap.put(ss[0], tempIndex++);
                  }
                  if (!geneNameMap.containsKey(ss[1])) {
                      geneNameMap.put(ss[1], tempIndex++);
                  }
                  tempEdgesList.add(ss);
              });
            
//            System.out.println(tempEdgesList);
            //build graph.
            EdgeWeightedGraph G = new EdgeWeightedGraph(geneNameMap.size());
            boolean withWeight = true;
            if (tempEdgesList.get(0).length == 2) {
                withWeight = false;
            }
            final boolean wWeight = withWeight;
            tempEdgesList.stream().forEach(s->{
                int v = (int) geneNameMap.get(s[0]);
                int w = (int) geneNameMap.get(s[1]);
                double weight = 1.0;
                if (wWeight) {
                    weight = Double.parseDouble(s[2]);
                }
                Edge e = new Edge(v, w, weight);
                G.addEdge(e);
            });
            tempEdgesList.clear();
            
            //compute the distance with knwon genes
            final int[] knwonGeneIDs = Arrays.stream(knownGenes)
                     .mapToInt(s->{
                            Object t = geneNameMap.get(s);
                            if(t == null){
                                System.err.println("Error can not find this node in Graph, please check: " + s);
                                System.exit(-1);
                            }
                                    
                            return (int)t;
                     })
                    .toArray();
            
            //System.out.println(G.toString());
            DecimalFormat decimalFormat = new DecimalFormat("#.####");
            if(topGenes != null){
                int[] topGeneIDs = Arrays.stream(topGenes)
                        .parallel()
                        .mapToInt(s->{
                            Object t = geneNameMap.get(s);
                            if(t == null){
                                System.err.println("Error can not find this node in Graph, please check: " + s);
                                System.exit(-1);
                            }
                                    
                            return (int)t;
                        })
                        .toArray();
                //System.out.println(Arrays.toString(knwonGeneIDs));
//                System.out.println(Arrays.toString(topGeneIDs));
                OptionalDouble results = averageMindistance(topGeneIDs, knwonGeneIDs, G);
                //System.out.println(decimalFormat.format(results));
                if (results.isPresent()) {
                    System.out.println(decimalFormat.format(results.getAsDouble()));
                }else{
                    System.out.println("NA");
                }
                
                System.exit(0);
            }
            
            //compute bootstrapping values.
            //candidate gene list, do not including known genes.
            HashSet<String> knownSet = new HashSet(Arrays.asList(knownGenes));
            String[] candidateGenes = geneNameMap.keySet().stream()
                    .map(e->e.toString())
                    .filter(s->!knownSet.contains(s))
                    .toArray(String[]::new);
//            System.out.println(Arrays.toString(candidateGenes));
            
            if (bootTimes>0 && RandomPickGenes>0) {
                //number of bootstrap.
                IntStream.range(0, bootTimes)
                        .sequential()
                        .forEach(s->{
                            int[] randomGenes = FisherYatesArrayShuffle.Shuffle(candidateGenes, RandomPickGenes)
                                    .stream()
                                    .mapToInt(e->geneNameMap.get(e))
                                    .toArray();
                            
                            //System.out.println(decimalFormat.format(averageMindistance(randomGenes, knwonGeneIDs, G)));
                            OptionalDouble results = averageMindistance(randomGenes, knwonGeneIDs, G);
                            //System.out.println(decimalFormat.format(results));
                            if (results.isPresent()) {
                                System.out.println(decimalFormat.format(results.getAsDouble()));
                            }else{
                                System.out.println("NA");
                            }
                        });
            }
    }
}