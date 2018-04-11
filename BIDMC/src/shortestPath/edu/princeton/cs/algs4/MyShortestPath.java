package shortestPath.edu.princeton.cs.algs4;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.stream.IntStream;
import org.docopt.Docopt;


/**
 * 
 * @author wavefancy@gmail.com
 * 
 * @version 1.8
 * Compute shortest path based on princeton's code base.
 * http://algs4.cs.princeton.edu/44sp/
 * all code on github: https://github.com/kevin-wayne/algs4
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
    private static int outputClosestNKnownGenes = -1; //output closest N knownGene for each topGenes.
    private static String TASK = "";
    private final static HashMap<Integer,String> reverseIDNameMap = new HashMap<>(geneNameMap.size());
    private static final DecimalFormat decimalFormat = new DecimalFormat("#.####");
    private static final HashMap<String, Integer> GENELENGTH_MAP = new  HashMap<>();
    private static double PENALTY4NOCONNECTION = 1000.0;
    private static boolean OUTPUTMATCHMIN = false; //output the minal value of gene with all_known_genes, to stdout.
    private static boolean withSteps = false; // if true, the distance between two nodes is the total_weight + number_of_steps.
                                              // if two nodes are not connectable, set the distance as the --penalty value.
	
    private static final String DOC =
                "MyUndirectedWeightedShortestPath.\n"
                + "\n"
                + "Usage:\n"
                + "  MyUndirectedWeightedShortestPath [-k knownGene] (-i topGenes [-c int] | [-r int -b int] [--task string]) [-g gene] [-p int] [-t cpus] [-l file] [-n] [--omatchmin] [--penalty f] [--step] \n"
                + "  MyUndirectedWeightedShortestPath (-h | --help)\n"
                + "  MyUndirectedWeightedShortestPath --version\n"
                + "\n"
                + "---------------------------\n"
                + "Read network from stdin, each line three columns: GeneA GeneB weight.\n"
                + "If only two columns, all edges were assigned the equal weight of 1.\n"
                + "Compute the avarage shorest path from topGenes with knownGenes. \n"
                + "mean(min_i(KnownGenes)), i iterate all topGenes.\n"
                + "min_i: the minimal distance between topGene_i with all knownGenes.\n"
                + "---------------------------\n"
                + "\n"
                + "Options:\n"
                + "  --task string Assign the task want this program to do!\n"
                + "                 path: use with -g\n"
                + "                 proportion: use with -k -b -r and -p, pick -r number of genes and bootstrap -b time.\n"
                + "                     Compute the proportion of times for each knwon gene, \n"
                + "                     which can be ranked to the top -p closest genes amoung all the query genes.\n"
                + "  -p int        The number of top genes to pick for task 'proportion'\n"
                + "  -g gene       Two gene names, output the path from gene A to gene B. input example: A,B \n"
                + "  -k knownGene  Input known gene list, one line one gene.\n"
                + "  -l file       Gene length annotation file. Output total gene length of top/random gene set.Two columns: geneName, length.\n"
                + "  -i topGenes   Input top gene list, one line one gene.\n"
                + "  -c int        Output the closest [int] knownGenes for each topGenes.\n"
                + "  -r int        Number of random picked genes for bootstrapping.\n"
                + "  -b int        Number of bootstrappings.\n"
                + "  -n            Output the gene name list for bootstrapping, in the first -r columns.\n"
                + "  --penalty f   Set the penality for no connection between 2 genes, default 1000.0.\n"
                + "  --omatchmin   Output the min value of a gene with known genes to stderr, default no output.\n"
                + "  --step        Set the distance is the total_edge_weight + num_of_steps, if no path between nodes, the distance is the --penalty value.\n"
                + "  -t cpus       Number of cpus for computing.\n"
                + "  -h --help     Show this screen.\n"
                + "  --version     Show version.\n"
                + "\n";
    
    /**
     * Get the average min distance between testGeneIDs with targetGeneIDs.
     * Average(min(test_i_ShortestDistanceToAllTargetGenes), iterate i))
     * This function will be run in parallel model.
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
//                                    .filter(k->{return sp.hasPathTo(k);})
                                    .mapToDouble(k->{
                                        if(sp.hasPathTo(k)){
                                            if (withSteps==true){
                                                int step = 0; //make the edge is the weight + number_of_steps.
                                                Iterable<Edge> edges = sp.pathTo(k);
                                                for (Edge edge : edges) {
                                                    step +=1;
                                                }
//                                                System.out.println(reverseIDNameMap.get(k)+":step: " + step);
                                                return sp.distTo(k) + step;
                                            }else{
                                                return sp.distTo(k);
                                            }
                                        }else{
                                            return PENALTY4NOCONNECTION;
                                        }
                                    })
                                    .toArray();
//                            if (minimalDis.length >0) {
//                                return Arrays.stream(minimalDis).min().getAsDouble();
//                            }else{
//                                return -1;
//                            }
                            double re = Arrays.stream(minimalDis).min().getAsDouble();
                            if (OUTPUTMATCHMIN) {
                                System.err.println(reverseIDNameMap.get(s)+"\t" + decimalFormat.format(re));
                            }
                            return re; //min distance to any target gene.
                        })
                        .average();
    }
    
    /**
     * Compute sorted distance to targetGenes and corresponding targetGene ID.
     * @param testGeneID
     * @param targetGeneIDs
     * @param G
     * @return list of ranked gene distance and geneIDs, ranked by distance to testGeneID.
     */
    private static List sortedDistanceAndID(int testGeneID, int [] targetGeneIDs, EdgeWeightedGraph G){
        DijkstraUndirectedSP sp = new DijkstraUndirectedSP(G, testGeneID);
        
        //get the distance to all targetGenes.
        double[] distance = new double[targetGeneIDs.length];
        for (int i = 0; i < targetGeneIDs.length; i++) {
            int targetGeneID = targetGeneIDs[i];
            if (sp.hasPathTo(targetGeneID)) {
                distance[i] = sp.distTo(targetGeneID);
            }else{
                distance[i] = Double.MAX_VALUE;
            }
        }
        
        //find the closest one.
        double[][] results = new double[distance.length][2]; //array of (distance, targetGeneID)
        for (int i = 0; i < results.length; i++) {
            double[] r = new  double[2];
            r[0] = distance[i];
            r[1] = targetGeneIDs[i];
            results[i] = r;
        }
        
        //rank by distance.
        List<double[]> reList = Arrays.asList(results);
        Collections.sort(reList,new Comparator<double[]>(){
                    public int compare(double[] i, double[] j){
                        double t =  i[0] - j[0];
                        if(t==0){
                            return 0;
                        }else if (t>0) {
                            return 1;
                        }else{
                            return -1;
                        }
                    }
                });
//        for (double[] ds : reList) {
//            System.out.println(Arrays.toString(ds));
//        }
        return reList;
    }
    
    
    /**
     * Compute the proportion of times for each known which were ranked as topN closest genes 
     * for each query gene among all query genes. 
     * @param geneIDs
     * @param topN 
     */
    private static void outputProportion(int[] geneIDs, final int topN, EdgeWeightedGraph G){
        final HashMap<Integer, Integer> knownGeneCountMap = new HashMap<>(knownGenes.length);
        Arrays.stream(knownGenes)
              .mapToInt(k->geneNameMap.get(k))
              .forEach(k->{
                  knownGeneCountMap.put(k, 0);
              });
                            
        int[] knownGeneIDs = Arrays.stream(knownGenes).mapToInt(k->geneNameMap.get(k)).toArray();

        Arrays.stream(geneIDs)
                .forEach(i->{
                   List<double[]> re = sortedDistanceAndID(i, knownGeneIDs, G);
                   
                   for (int j = 0; j < topN; j++) {
                       int id = (int) re.get(j)[1];
                       knownGeneCountMap.put(id, knownGeneCountMap.get(id) + 1);
                   }
                });
        Arrays.stream(geneIDs).mapToObj(s->reverseIDNameMap.get(s)).forEach(s->System.out.println(s));
        //System.out.println(knownGeneCountMap);
        //output resutls.
        StringBuilder sb = new  StringBuilder();
        for (Map.Entry<Integer, Integer> entry : knownGeneCountMap.entrySet()) {
            Integer key = entry.getKey();
            Integer value = entry.getValue();
            sb.setLength(0);
            sb.append(reverseIDNameMap.get(key)).append("\t");
            sb.append(decimalFormat.format(value * 1.0 / geneIDs.length));
            System.out.println(sb.toString());
        }
    }

    public static void main(String[] args) {
            Map<String, Object> opts =
                 new Docopt(DOC).withVersion("1.8").parse(args);
//		     System.err.println(opts);
            if(opts.get("-t") != null){
                    System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", (String) opts.get("-t"));
            }
            
            if(opts.get("-r") != null){
                RandomPickGenes = Integer.parseInt((String) opts.get("-r"));
            }
            if(opts.get("-b") != null){
                bootTimes = Integer.parseInt((String) opts.get("-b"));
            }
            if(opts.get("--task") != null){
                TASK = (String) opts.get("--task");
            }
            //set the penalty if no connection between two genes.
            if(opts.get("--penalty") != null){
                PENALTY4NOCONNECTION = Double.parseDouble((String)opts.get("--penalty"));
            }
            //output the minvalue of a gene with known genes.
            if((boolean)opts.get("--omatchmin") == true){
                OUTPUTMATCHMIN = true;
            }
            
            //Make the distance with steps.
            if((boolean)opts.get("--step") == true){
                withSteps = true;
            }
            
            
            
      
            boolean T_OUTPUTBOOTGENELIST = false;
            if((boolean)opts.get("-n") == true){
                T_OUTPUTBOOTGENELIST = true;
            }
            final boolean OUTPUTBOOTGENELIST = T_OUTPUTBOOTGENELIST;
            
//            System.out.println(Arrays.toString(topGenes));
//            System.out.println(Arrays.toString(knownGenes));

            //read network from stdin.
            BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
            final LinkedList<String[]> tempEdgesList = new LinkedList();
            in.lines()
              .filter(s -> s.length() > 0)
              .sequential() //*** Very Important **** Here must to be sequential, because we will update index.
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
            
            //read knwon genes.
            if(opts.get("-k") != null){
                try {
                   knownGenes = Files.lines(Paths.get((String) opts.get("-k")))
                        .filter(s -> s.length() > 0)
                        .toArray(String[]::new);
                   
                   //checking whether in network.
                   boolean error = false;
                   for (String string: knownGenes) {
                        if (!geneNameMap.containsKey(string)) {
                            System.err.println("Can not find this knownGenes in graph: " + string);
                            error = true;
                        }
                   }
                    if (error) {
                        System.exit(-1);
                    }
                   
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
                   
                   boolean error = false;
                   for (String string: topGenes) {
                        if (!geneNameMap.containsKey(string)) {
                            System.err.println("Can not find this topGenes in graph: " + string);
                            error = true;
                        }
                   }
                   if (error) {
                        System.exit(-1);
                   }
                   
                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
            }
            //set this after -c.
            if(opts.get("-c") != null){
                outputClosestNKnownGenes = Integer.parseInt((String) opts.get("-c"));
                if(topGenes == null){
                    System.err.println("ERROR: -c option must be used in combination with -i option!");
                    System.exit(-1);
                }
            }
            //read gene length annotation. 
            if(opts.get("-l") != null){
                try {
                   Files.lines(Paths.get((String) opts.get("-l")))
                        .filter(s -> s.length() > 0)
                        .forEach(s->{
                            String[] ss = s.split("\\s+");
                            
                            GENELENGTH_MAP.put(ss[0], Integer.valueOf(ss[1]));
                        });
                    
                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
            }
            
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
            
            for (Map.Entry<String, Integer> entry : geneNameMap.entrySet()) {
                String key = entry.getKey();
                Integer value = entry.getValue();
                reverseIDNameMap.put(value, key);
            }
            
            //perform different tasks from here.
            if (TASK.equalsIgnoreCase("path")) {
                String[] genes = null;
                if(opts.get("-g") != null){
                    genes = ((String) opts.get("-g")).split(",");
                }else{
                    System.out.println("ERROR: TASK 'path' should be work with option '-g!");
                    System.exit(-1);
                }
                
                boolean chxError = false;
                for (String gene : genes) {
                    if (!geneNameMap.containsKey(gene)) {
                        System.err.println("ERROR: can not find this node in network: " + gene);
                        chxError = true;
                    }
                }
                if (chxError) {
                    System.exit(-1);
                }
               
                int[] gIDs = Arrays.stream(genes)
                        .mapToInt(s->{
                            return geneNameMap.get(s);
                        })
                        .toArray();
                
                DijkstraUndirectedSP sp = new DijkstraUndirectedSP(G, gIDs[0]);
                StringBuilder sb = new StringBuilder();
                if (sp.hasPathTo(gIDs[1])) {
                    Iterable<Edge> edges = sp.pathTo(gIDs[1]);
                    for (Edge edge : edges) {
                        sb.setLength(0);
                        sb.append(reverseIDNameMap.get(edge.either())).append("\t");
                        sb.append(reverseIDNameMap.get(edge.other(edge.either()))).append("\t");
                        sb.append(decimalFormat.format(edge.weight()));
                        System.out.println(sb.toString());
                    }
                }else{
                    System.err.println("There is no path between:" + Arrays.toString(genes));
                }
                System.exit(-1);
            }
            
            //Compute the proportion of times for each knwon gene, "
            //   which can be ranked to the top -p closest genes amoung all the query genes."
            if(TASK.equalsIgnoreCase("proportion")){
                int temp = -1;
                if(opts.get("-p") != null){
                    temp = Integer.parseInt((String) opts.get("-p"));
                }else{
                    System.out.println("ERROR: TASK 'proportion' should be work with option '-p!");
                    System.exit(-1);
                }
                final int topN = temp;
                if (knownGenes == null) {
                    System.out.println("ERROR: TASK 'proportion' should be work with option '-k!");
                    System.exit(-1);
                }
                if (bootTimes <=0) {
                    System.out.println("ERROR: TASK 'proportion' should be work with option '-b!");
                    System.exit(-1);
                }
                if (RandomPickGenes <= 0) {
                    System.out.println("ERROR: TASK 'proportion' should be work with option '-r!");
                    System.exit(-1);
                }
                
                //compute bootstrapping values.
                //candidate gene list, do not including known genes.
                HashSet<String> knownSet = new HashSet(Arrays.asList(knownGenes));
                String[] candidateGenes = geneNameMap.keySet().stream()
                        .map(e->e.toString())
                        .filter(s->!knownSet.contains(s))
                        .toArray(String[]::new);
                
                //bootstrap for compute proportion:
                IntStream.range(0, bootTimes)
                        .parallel()
                        //.sequential()
                        .forEach(s->{
                            int[] randomGeneIDs = FisherYatesArrayShuffle.Shuffle(candidateGenes, RandomPickGenes)
                                    .stream()
                                    .mapToInt(e->geneNameMap.get(e))
                                    .toArray();
                            //System.out.println(Arrays.toString(randomGeneIDs));
                            //output porportion.
                            outputProportion(randomGeneIDs, topN, G);
                        });
                System.exit(-1);
            }
            
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

                //output the distance between top genes and known genes.
                if (outputClosestNKnownGenes <= 0) {
                    
                    OptionalDouble results = averageMindistance(topGeneIDs, knwonGeneIDs, G);
                    //System.out.println(decimalFormat.format(results));
                    //collect total gene length.
                    if (GENELENGTH_MAP.size()>0) {
                        int totalLen = Arrays.stream(topGeneIDs)
                            .mapToObj(x->reverseIDNameMap.get(x))
                            .mapToInt(x->GENELENGTH_MAP.get(x))
                            .sum();
                        System.out.print(Integer.toString(totalLen) + "\t");
                    }
                    
                    
                    if (results.isPresent()) {
                        System.out.println(decimalFormat.format(results.getAsDouble()));
                    }else{
                        System.out.println("NA");
                    }

                   
                }else{ //output the closest top n known genes for each top genes.
                    
//                    System.out.println(outputClosestNKnownGenes);
//                    System.out.println("shortestPath.edu.princeton.cs.algs4.MyShortestPath.main()");
                    Arrays.stream(topGeneIDs)
                            .parallel()
                            .forEach(s->{
                                List<double[]> re = sortedDistanceAndID(s, knwonGeneIDs, G);
                                StringBuilder sb = new StringBuilder();
                                DecimalFormat formater = new  DecimalFormat("#.####");
//                                System.out.println("shortestPath.edu.princeton.cs.algs4.MyShortestPath.main()");
                                int len = Math.min(outputClosestNKnownGenes, re.size());
                                for (int i = 0; i < len; i++) {
                                    sb.setLength(0);
                                    sb.append(i+1).append("\t");
                                    sb.append(reverseIDNameMap.get(s)).append("\t");
                                    double[] r = re.get(i);
                                    if(r[0] == Double.MAX_VALUE){
                                        sb.append("NA").append("\t");
                                        sb.append("NA");
                                    }else{
                                        sb.append(reverseIDNameMap.get((int)r[1])).append("\t");
                                        sb.append(formater.format(r[0]));
                                    }
                                    System.out.println(sb.toString());
                                }
                            });
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
            
            // For random bootstap.
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
                            
                            if (OUTPUTBOOTGENELIST){
                                for (int rg : randomGenes) {
                                    System.out.print(reverseIDNameMap.get(rg) + "\t");
                                }
                            }
                            
                            //System.out.println(decimalFormat.format(results));
                            if (GENELENGTH_MAP.size()>0) {
                                int totalLen = Arrays.stream(randomGenes)
                                    .mapToObj(x->reverseIDNameMap.get(x))
                                    .mapToInt(x->{
                                        if(GENELENGTH_MAP.get(x) == null){
                                            System.err.println("ERROR: Can not find length annotation for gene: " + x);
                                            System.exit(-1);
                                            return -1;
                                        }else{
                                            return GENELENGTH_MAP.get(x);
                                        }
                                    })
                                    .sum();
                                System.out.print(Integer.toString(totalLen) + "\t");
                            }
                            
                            if (results.isPresent()) {
                                System.out.println(decimalFormat.format(results.getAsDouble()));
                            }else{
                                System.out.println("NA");
                            }
                        });
            }
    }
}