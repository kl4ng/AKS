/* Kevin Lang
 * Evan Dorundo
 */

import java.math.*;


public class AKS
{
	public static final BigDecimal EPSILON = new BigDecimal(1E-9);
	
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
	
	private static boolean isPowerOfInteger(BigInteger n)
	{
		BigInteger maxB = BigInteger.valueOf(n.bitLength());
		for(BigInteger b = BigInteger.valueOf(2); b.compareTo(maxB) < 0; b = b.add(BigInteger.ONE))
		{
			BigInteger lo = BigInteger.valueOf(2);
			BigInteger hi = n.add(BigInteger.ZERO);
			
			while(lo.compareTo(hi) < 0)
			{
				BigInteger a = hi.add(lo).divide(BigInteger.valueOf(2));
				BigInteger tmp = pow(a,b);
				
				if(tmp.compareTo(n) == 0)
					return true;
				else if(tmp.compareTo(n) < 0)
					lo = a.add(BigInteger.ONE);
				else
					hi = a;
			}
			
		}
		
		return false;
	}
	
	private static BigInteger multOrder(BigInteger n)
	{
		BigInteger maxK = log(n, 2).multiply(log(n, 2)).add(EPSILON).toBigInteger();
		
		for(BigInteger r = new BigInteger("2"); ; r = r.add(BigInteger.ONE))
		{
			BigInteger a = BigInteger.ONE;
			boolean valid = true;
			for(BigInteger k = BigInteger.ONE; k.compareTo(maxK) <= 0; k = k.add(BigInteger.ONE))
			{
				a = a.multiply(n).mod(r);
				
				if(a.compareTo(BigInteger.ONE) <= 0)
				{
					valid = false;
					break;
				}
			}
			
			if(valid)
			{
				return r;
			}
		}
	}
	
	//mod euclidian gcd algorithm
	private static BigInteger gcd(BigInteger a, BigInteger b)
	{
		if(b.equals(BigInteger.valueOf(0)))
			return a;
			
		return gcd(b, a.mod(b));
	}
	
	private static BigDecimal log(BigInteger a, int b)
	{
		
	}
	
	private static BigInteger pow(BigInteger a, BigInteger b)
	{
		if(b.compareTo(BigInteger.ONE) == 0)
			return a;
		
		BigInteger tmp = pow(a, b.shiftRight(1));
		tmp = tmp.multiply(tmp);
		
		//check if b is odd
		if(b.and(BigInteger.ONE).compareTo(BigInteger.ONE) == 0)
			tmp = tmp.multiply(a);
		
		return tmp;
	}
	
	public static void main(String[] args)
	{
		//System.out.println(pow(new BigInteger("7"), new BigInteger("3")));
		System.out.println(isPowerOfInteger(new BigInteger("6436343")));
	}
}
