/* Kevin Lang
 * Evan Dorundo
 */

import java.math.*;
import java.util.*;

/*
 * (1) BigIntegers
 * (2) longs for r
 * (3) longs for a
 * (4) cast n to long
 * (5) longs for r, BigInteger for n
 */

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
		long r = multOrder(n);
		
		//(3) see if 1 < gcd(a,n) < n for some a<=r
		if(!checkGCD(n, r))
		{
			return false;
		}
		
		//(4) if n <= r, we know its prime
		if(n.compareTo(BigInteger.valueOf(r)) <= 0)
			return true;
		
		//(5) 
		checkCondition(n, r);
		
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
	
	private static long multOrder(BigInteger n)
	{
		BigInteger maxK = log(n, 2).multiply(log(n, 2)).add(EPSILON).toBigInteger();
		
		for(long r = 2; ; r++)
		{
			BigInteger a = BigInteger.ONE;
			boolean valid = true;
			for(BigInteger k = BigInteger.ONE; k.compareTo(maxK) <= 0; k = k.add(BigInteger.ONE))
			{
				a = a.multiply(n).mod(BigInteger.valueOf(r));
				
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
	
	private static boolean checkGCD(BigInteger n, long r)
	{
		for(long a = 2; a <= r; a++)
		{
			BigInteger gcd = gcd(BigInteger.valueOf(a), n);
			if(gcd.compareTo(BigInteger.ONE) > 0 && gcd.compareTo(n) < 0)
				return false;
		}
		
		return true;
	}
	
	private static boolean checkCondition(BigInteger n, long r)
	{
		BigInteger maxA = Math.sqrt(totient(r)).multiply(log(r, 2)).toBigInteger();
		for(BigInteger a = BigInteger.ONE; a.compareTo(maxA) <= 0; a = a.add(BigInteger.ONE))
		{
			
		}
	}
	
	//mod euclidian gcd algorithm
	private static BigInteger gcd(BigInteger a, BigInteger b)
	{
		if(b.equals(BigInteger.valueOf(0)))
			return a;
			
		return gcd(b, a.mod(b));
	}
	
	private static double log(BigInteger a, int b)
	{
		String s = a.toString();
		int l = Math.min(18, s.length());
		double d = Double.parseDouble(s.substring(0,l)) / Math.pow(10, l-1);
		return (Math.log10(d) + (s.length()-1)) / Math.log10(b);
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
	
	private static long totient(long r)
	{
		long result = r;
		for(long a = 2; a * a <= r; a++)
		{
			if(r % a == 0)
			{
				result /= a;
				result *= a - 1;
				
				while(r % a == 0)
				{
					r /= a;
				}
			}
		}
		
		if(r > 1)
		{
			result /= r;
			result *= r - 1;
		}
		
		return result;
	}
	
	private class Polynomial
	{
		ArrayList<BigInteger> coef;
		BigInteger r;
		BigInteger n;
		
		public Polynomial(BigInteger R, BigInteger N, BigInteger a)
		{
			coef = new ArrayList<BigInteger>();
			r = R;
			n = N;
			
			coef.add(BigInteger.ONE);
			coef.add(a);
		}
		
		private Polynomial pow(BigInteger e)
		{
			if(e.compareTo(BigInteger.ONE) == 0)
			{
				return this;
			}
			
			Polynomial temp = pow(e.shiftRight(1));
			temp = temp.multiply(temp);
			
			if(e.and(BigInteger.ONE).compareTo(BigInteger.ONE) == 0)
			{
				temp = temp.multiply(this);
			}
			
			return temp;
		}
		
		private Polynomial multiply(Polynomial p)
		{
			
		}
	}
	
	public static void main(String[] args)
	{
		//System.out.println(pow(new BigInteger("7"), new BigInteger("3")));
		System.out.println(log(new BigInteger("6"), 2));
	}
}
