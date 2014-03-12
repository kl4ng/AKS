/* Kevin Lang
 * Evan Dorundo
 */

import java.math.BigInteger;


public class AKS
{
	public static boolean isPrime(int num)
	{
		return isPrime(BigInteger.valueOf(num));
	}
	
	public static boolean isPrime(BigInteger n)
	{
		// (1)see if n = a^b 
		if(isPowerOfInteger(n))
			return false;
		
		//(2) find smallest r s.t. Or(n) > (log n)^2
		BigInteger r = multOrder(n);
		
		//(3) see if 1 < gcd(a,n) < n for some a<=r
		BigInteger a = BigInteger.valueOf(2);
		BigInteger one = BigInteger.valueOf(1);
		
		for(; a.compareTo(r) <= 0; a.add(one))
		{
			BigInteger gcd = gcd(a, r);
			if(gcd.compareTo(one) > 0 && gcd.compareTo(n) < 0)
				return false;
		}
		
		//(4) if n <= r, we know its prime
		if(n.compareTo(r) <= 0)
			return true;
		
		//(5) 
		
		//(6) we know its not composite by this point
		return true;
	}
	
	private static boolean isPowerOfInteger(BigInteger num)
	{
		
	}
	
	private static BigInteger multOrder(BigInteger num)
	{
		
	}
	
	//mod euclidian gcd algorithm
	private static BigInteger gcd(BigInteger a, BigInteger b)
	{
		if(b.equals(BigInteger.valueOf(0)))
			return a;
			
		return gcd(b, a.mod(b));
	}
	
	public static void main(String[] args)
	{
		
	}
}
