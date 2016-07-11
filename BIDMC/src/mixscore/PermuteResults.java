package mixscore;

import java.text.DecimalFormat;

/**
 * Contain permutation results.
 * @author wallace
 *
 */
public class PermuteResults {
	int smallerNum = 0;
	int largerNum = 0;
	double original = 0.0;
	
	//check results when add 100 points. If significance less than checkThreshold, stop permutation.
	int checkCycle = 100; 
	double checkThreshold = 0.2;
	int cTemp = 0;
	
	DecimalFormat formater = new DecimalFormat("#.####");
	
	/**
	 * Add current points, and check whether to continue permutation.
	 * @param value
	 * @return true: continue permutation, else stop permutation.
	 */
	public synchronized boolean addAndCheckLarger(double value) {
//		System.err.println(value);
		cTemp ++;
		if (value >= original) {
			largerNum ++;
		}else {
			smallerNum ++;
		}
		if (cTemp == 100) {
			cTemp = 0;
			if (largerNum*1.0 /(largerNum + smallerNum) >= checkThreshold) {
				return false;
			}
		}
		return true;
	}
	
	public void setOriginal(double ori) {
//		System.err.println("ori: "+  ori);
		this.original = ori;
	}
	
	@Override
	public String toString() {
		return formater.format(original) + "\t" + smallerNum + "\t" + largerNum;
	}
}
