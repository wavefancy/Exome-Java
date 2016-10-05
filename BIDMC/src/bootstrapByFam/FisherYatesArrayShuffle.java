/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bootstrapByFam;

import java.util.concurrent.ThreadLocalRandom;

/**
 * Shuffle array by Fisher-Yates Algorithm.
 * https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 * http://www.programming-algorithms.net/article/43676/Fisher-Yates-shuffle
 * -- To shuffle an array a of n elements (indices 0..n-1):
 * for i from n−1 downto 1 do
 *    j ← random integer such that 0 ≤ j ≤ i
 *    exchange a[j] and a[i]
 * @author wallace
 */
public class FisherYatesArrayShuffle {
    
    /**
     * Shuffle Array in place.
     * @param inArr 
     */
    public static void Shuffle(String[] inArr){
        ThreadLocalRandom random = ThreadLocalRandom.current();
        for (int i = inArr.length -1 ; i > 0; i--) {
            int index = random.nextInt(i+1);
            String t = inArr[i];
            inArr[i] = inArr[index];
            inArr[index] = t;
            
        }
    }
    
    public static void main(String[] args) {
        String[] temp = {"1","2","3","4"};
        FisherYatesArrayShuffle.Shuffle(temp);
        for (String string : temp) {
            System.out.println(string);
        }
    }
}
