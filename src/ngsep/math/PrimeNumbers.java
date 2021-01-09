package ngsep.math;

public class PrimeNumbers {
	private int capacity = 1000000;
	
	private int [] nextPrimes;
	private boolean [] prime;
	
	public PrimeNumbers (int capacity) {
		this.capacity = capacity;
		initialize();
		
	}
	public PrimeNumbers () {
		initialize();
	}
	private void initialize() {
		nextPrimes = new int [capacity+1];
		prime = new boolean [capacity+1];
		int p = 2;
		for(int i=0;i<=capacity;i++) {
			prime[i] = (i == p);
			if(prime[i]) p = findNextPrime (p);
			nextPrimes[i] = p;
		}
	}
	private int findNextPrime(int n) {
		int j=n+1;
		if(j%2==0) j++;
		while (true) {
			boolean prime = checkPrime(j);
			if(prime) break;
			j+=2;
		}
		return j;
	}
	private boolean checkPrime(int number) {
		if(number <2) return false;
		if(number ==2) return true;
		if(number%2==0) return false;
		double q = Math.sqrt(number);
		for(int k=3;k<q+1;k+=2) {
			if(number%k==0) return false;
		}
		return true;
	}
	public int getNextPrime (int number) {
		if(number > capacity) findNextPrime(number);
		return nextPrimes[number];
	}
	public boolean isPrime (int number) {
		if(number > capacity) return checkPrime(number);
		return prime[number];
	}
	
	public int getCapacity() {
		return capacity;
	}
	public static void main(String[] args) {
		long t1 = System.currentTimeMillis();
		new PrimeNumbers(10000000);
		long t2 = System.currentTimeMillis();
		long sec = (t2-t1)/1000;
		System.out.println("Prime numbers created. Time: "+sec);
	}
}
