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

// TEST CASES: 27644437

public class AKS
{
	public static final BigDecimal EPSILON = new BigDecimal(1E-9);
	
	public static boolean isPrime(int n)
	{
		return isPrime(BigInteger.valueOf(n));
	}
	
	public static boolean isPrime(long n)
	{
		return isPrime(BigInteger.valueOf(n));
	}
	
	public static boolean isPrime(BigInteger n)
	{
		System.out.println(Double.MIN_NORMAL);
		
		// (1)see if n = a^b 
		if(isPowerOfInteger(n))
			return false;
		
		//(2) find smallest r s.t. Or(n) > (log n)^2
		long r = multOrder(n);
		
		//(3) see if 1 < gcd(a,n) < n for some a<=r
		if(!checkGCD(n, r))
			return false;
		
		//(4) if n <= r, we know its prime
		if(n.compareTo(BigInteger.valueOf(r)) <= 0)
			return true;
		
		System.out.println(r);
		
		//(5) 
		if(!checkCondition(n, r))
			return false;
		
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
		long maxK = (long)Math.floor(log(n, 2) * log(n, 2));
		
//		BigInteger maxK = log(n, 2).multiply(log(n, 2)).add(EPSILON).toBigInteger();
		
		for(long r = 2; ; r++)
		{
			BigInteger a = BigInteger.ONE;
			boolean valid = true;
			for(long k = 1; k <= maxK; k++)
			{
				a = a.multiply(n).mod(BigInteger.valueOf(r));
				
				if(a.compareTo(BigInteger.ONE) <= 0)
				{
					valid = false;
					break;
				}
			}
			
//			for(BigInteger k = BigInteger.ONE; k.compareTo(maxK) <= 0; k = k.add(BigInteger.ONE))
//			{
//				a = a.multiply(n).mod(BigInteger.valueOf(r));
//				
//				if(a.compareTo(BigInteger.ONE) <= 0)
//				{
//					valid = false;
//					break;
//				}
//			}
			
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
		long maxA = (long)Math.floor(Math.sqrt(totient(r)) * log(n, 2));
		for(long a = 1; a <= maxA; a++)
		{
			System.out.println(a);
			
			Polynomial rhs = new Polynomial(r, n);
			rhs.coef[(int)(n.mod(BigInteger.valueOf(r)).longValue())] = BigInteger.ONE;
			rhs.coef[0] = rhs.coef[0].add(BigInteger.valueOf(a));
			
			if(!(new Polynomial(r, n, a).pow(n).equals(rhs)))
			{
				return false;
			}
		}
		
		return true;
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
	
	public static BigComplex[] fft(BigInteger[] a, int m, BigComplex w)
	{
		BigComplex[] newCoef = new BigComplex[a.length];
		for(int i = 0; i < a.length; i++)
		{
			newCoef[i] = new BigComplex(new BigDecimal(a[i]), BigDecimal.ZERO);
		}
		
		return fft(newCoef, m, w);
	}
	
	public static BigComplex[] fft(BigComplex[] a, int m, BigComplex w)
	{
		if(m == 1)
		{
			BigComplex[] f = new BigComplex[1];
			f[0] = a[0];
			
			return f;
		}
		else
		{
			BigComplex[] aEven = new BigComplex[m >> 1];
			BigComplex[] aOdd = new BigComplex[m >> 1];
			
			for(int i = 0; i < a.length; i++)
			{
				if((i & 1) == 0)
				{
					aEven[i >> 1] = a[i];
				}
				else
				{
					aOdd[i >> 1] = a[i];
				}
			}
			
			BigComplex[] fEven = fft(aEven, m >> 1, w.multiply(w));
			BigComplex[] fOdd = fft(aOdd, m >> 1, w.multiply(w));
			
			BigComplex[] f = new BigComplex[m];
			
			BigComplex x = new BigComplex(BigDecimal.ONE, BigDecimal.ZERO);
			for(int j = 0; j < m >> 1; j++)
			{
				f[j] = fEven[j].add(x.multiply(fOdd[j]));
				
				f[j + (m >> 1)] = fEven[j].subtract(x.multiply(fOdd[j]));
				
				x = x.multiply(w);
			}
			
			return f;
		}
	}
	
	public static Complex[] fft(BigInteger[] a, int m, Complex w)
	{
		Complex[] newCoef = new Complex[a.length];
		for(int i = 0; i < a.length; i++)
		{
			newCoef[i] = new Complex(a[i].doubleValue(), 0);
		}
		
		return fft(newCoef, m, w);
	}
	
	public static Complex[] fft(Complex[] a, int m, Complex w)
	{
		if(m == 1)
		{
			Complex[] f = new Complex[1];
			f[0] = a[0];
			
			return f;
		}
		else
		{
			Complex[] aEven = new Complex[m >> 1];
			Complex[] aOdd = new Complex[m >> 1];
			
			for(int i = 0; i < a.length; i++)
			{
				if((i & 1) == 0)
				{
					aEven[i >> 1] = a[i];
				}
				else
				{
					aOdd[i >> 1] = a[i];
				}
			}
			
			Complex[] fEven = fft(aEven, m >> 1, w.multiply(w));
			Complex[] fOdd = fft(aOdd, m >> 1, w.multiply(w));
			
			Complex[] f = new Complex[m];
			
			Complex x = new Complex(1, 0);
			for(int j = 0; j < m >> 1; j++)
			{
				f[j] = fEven[j].add(x.multiply(fOdd[j]));
				
				f[j + (m >> 1)] = fEven[j].subtract(x.multiply(fOdd[j]));
				
				x = x.multiply(w);
			}
			
			return f;
		}
	}
	
	private static class Polynomial
	{
		BigInteger[] coef;
		long r;
		BigInteger n;
		
		public Polynomial(long R, BigInteger N)
		{
			coef = new BigInteger[(int)R];
			for(int i = 0; i < coef.length; i++)
			{
				coef[i] = BigInteger.ZERO;
			}
			
			r = R;
			n = N;
		}
		
		public Polynomial(long R, BigInteger N, long a)
		{
			coef = new BigInteger[(int)R];
			for(int i = 0; i < coef.length; i++)
			{
				coef[i] = BigInteger.ZERO;
			}
			
			r = R;
			n = N;
			
			coef[1] = BigInteger.ONE;
			coef[0] = BigInteger.valueOf(a);
		}
		
		public boolean equals(Polynomial p)
		{
			if(coef.length != p.coef.length)
			{
				return false;
			}
			
			for(int i = 0; i < coef.length; i++)
			{
				if(!coef[i].equals(p.coef[i]))
				{
					return false;
				}
			}
			
			return true;
		}
		
		public Polynomial pow(BigInteger e)
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
		
		private Polynomial expand()
		{
			int tmp = (int)((r << 1) - 1);
			int pow2 = 0;
			while(tmp != 0)
			{
				tmp >>= 1;
				pow2++;
			}
			
			if(Long.bitCount((r << 1) - 1) == 1)
			{
				pow2--;
			}
			
			Polynomial newPoly = new Polynomial(1 << pow2, n);
			for(int i = 0; i < newPoly.coef.length; i++)
			{
				if(i < r)
					newPoly.coef[i] = coef[i];
				else
					newPoly.coef[i] = BigInteger.ZERO;
			}
			
			return newPoly;
		}
		
		// Assumes both polynomials are the same degree
		private Polynomial multiply(Polynomial p)
		{
			Polynomial q = expand();
			p = p.expand();
			int m = (int) Math.max(q.r, p.r);
			
			System.out.println(Math.sin(2.0 * Math.PI / (double)m));
			
			Complex root = new Complex(Math.cos(2.0 * Math.PI / (double)m), Math.sin(2.0 * Math.PI / (double)m));
			
			Complex[] f1 = fft(q.coef, m, root);
			Complex[] f2 = fft(p.coef, m, root);
			
			Complex[] f3 = new Complex[m];
			for(int i = 0; i < m; i++)
			{
				f3[i] = f1[i].multiply(f2[i]);
			}
			
			Complex[] c = fft(f3, m, new Complex(1, 0).divide(root));
			
			Polynomial result = new Polynomial(c.length, n);
			for(int j = 0; j < c.length; j++)
			{
				result.coef[j] = BigInteger.valueOf(Math.round(c[j].real + 1E-9));
//				System.out.println(c[j].imag);
			}
			
			return result.mod((int)r);
			
//			BigComplex root = new BigComplex(Math.cos(2 * Math.PI / m), Math.sin(2 * Math.PI / m));
//			
//			BigComplex[] f1 = fft(q.coef, m, root);
//			BigComplex[] f2 = fft(p.coef, m, root);
//			
//			BigComplex[] f3 = new BigComplex[m];
//			for(int i = 0; i < m; i++)
//			{
//				f3[i] = f1[i].multiply(f2[i]);
//			}
//			
//			BigComplex[] c = fft(f3, m, new BigComplex(BigDecimal.ONE, BigDecimal.ZERO).divide(root));
//			
//			Polynomial result = new Polynomial(c.length, n);
//			for(int j = 0; j < c.length; j++)
//			{
//				result.coef[j] = c[j].real.add(EPSILON).toBigInteger();
//				System.out.println(c[j].imag);
//			}
//			
//			return result.mod((int)r);
		}
		
		public Polynomial mod(int m)
		{
			Polynomial p = new Polynomial(m, n);
			for(int i = 0; i < coef.length; i++)
			{
				p.coef[i % m] = p.coef[i % m].add(coef[i]).mod(n);
			}
			
			return p;
		}
		
		// We do the mod during the multiplication
//		private Polynomial multiply(Polynomial p)
//		{
//			Polynomial ret = new Polynomial(r, n);
//			for(int i = 0; i < coef.length; i++)
//			{
//				for(int j = 0; j < p.coef.length; j++)
//				{
//					BigInteger c = coef[i].multiply(p.coef[j]).mod(n);
//					
//					ret.coef[(i + j) % (int)r] = ret.coef[(i + j) % (int)r].add(c).mod(n);
//				}
//			}
//			
//			return ret;
//		}
	}
	
	private static class Complex
	{
		double real;
		double imag;
		
		public Complex(double r, double i)
		{
			real = r;
			imag = i;
			
//			System.out.println("Here");
			
//			if(real <= 1E-307 || imag <= 1E-307)
//			{
//				System.out.println("Underflow!");
//			}
		}
		
		public Complex add(Complex c)
		{
			return new Complex(real + c.real, imag + c.imag);
		}
		
		public Complex subtract(Complex c)
		{
			return new Complex(real - c.real, imag - c.imag);
		}
		
		public Complex multiply(Complex c)
		{
			return new Complex(real * c.real - imag * c.imag, real * c.imag + imag * c.real);
		}
		
		public Complex divide(Complex c)
		{
			double denom = c.real * c.real + c.imag * c.imag;
			return new Complex((real * c.real + imag * c.imag) / denom, (imag * c.real - real * c.imag) / denom);
		}
	}
	
	private static class BigComplex
	{
		BigDecimal real;
		BigDecimal imag;
		
		public BigComplex(double r, double i)
		{
			real = new BigDecimal(r);
			imag = new BigDecimal(i);
		}
		
		public BigComplex(BigDecimal r, BigDecimal i)
		{
			real = r;
			imag = i;
		}
		
		public BigComplex add(BigComplex c)
		{
			return new BigComplex(real.add(c.real), imag.add(c.imag));
		}
		
		public BigComplex subtract(BigComplex c)
		{
			return new BigComplex(real.subtract(c.real), imag.subtract(c.imag));
		}
		
		public BigComplex multiply(BigComplex c)
		{
			return new BigComplex(real.multiply(c.real).subtract(imag.multiply(c.imag)),
					              real.multiply(c.imag).add(imag.multiply(c.real)));
		}
		
		public BigComplex divide(BigComplex c)
		{
			BigDecimal denom = c.real.multiply(c.real).add(c.imag.multiply(c.imag));
			return new BigComplex(real.multiply(c.real).add(imag.multiply(c.imag)).divide(denom, RoundingMode.HALF_EVEN),
								  imag.multiply(c.real).subtract(real.multiply(c.imag)).divide(denom, RoundingMode.HALF_EVEN));
		}
	}
	
	public static void main(String[] args)
	{
		Scanner in = new Scanner(System.in);
		
		while(true)
		{
			System.out.println(isPrime(new BigInteger(in.next())));
		}
	}
}
