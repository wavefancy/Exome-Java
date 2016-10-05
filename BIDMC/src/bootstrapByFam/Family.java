/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bootstrapByFam;

import java.util.Arrays;
import java.util.LinkedList;
import org.apache.commons.math3.util.Pair;

/**
 *
 * @author wallace
 */
public class Family {
    String famName = "";
    //geneName, meta information for this gene.
    //like: ZFYVE9, chr1    52703445        A       T 
    //Pair<String, String>[] geneListPairs = null;
    //LinkedList<Pair<String, String>> geneList = new LinkedList<>();
    LinkedList<String> geneList = new LinkedList<>();
    private String[] geneArray = null;

    public Family(String name) {
        this.famName = name;
    }

    public void addGene(String gName){
        //geneList.add(Pair.create(gName, metaString));
        geneList.add(gName);
    }
    
    public void listToArray(){
        geneArray = geneList.toArray(new String[0]);
        geneList = null;
    }
    
    /**
     * Make a deep copy of gene array.
     * @return 
     */
    public String[] getGeneArrayCopy(){
        return Arrays.copyOf(geneArray, geneArray.length);
    }
}
