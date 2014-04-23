/* Kevin Lang
 * Evan Dorundo
 */

import java.util.*;
import java.io.*;
import java.math.*;

public class AKS
{
	/*
	 * Values for MODE:
	 * 0: Slow multiplication
	 * 1: Fast fourier transform with BigDecimals
	 * 2: Fast fourier transform with doubles
	 */
	public static final int MODE = 2;
	
	public static final MathContext PRECISION = MathContext.DECIMAL64;
	public static final BigDecimal EPSILON = new BigDecimal(1E-3, PRECISION);
	
	//wrapper function for testing integers
	public static boolean isPrime(int n)
	{
		return isPrime(BigInteger.valueOf(n));
	}
	
	//wrapper function for testing longs
	public static boolean isPrime(long n)
	{
		return isPrime(BigInteger.valueOf(n));
	}
	
	public static boolean isPrime(BigInteger n)
	{
		//(1) see if n = a^b 
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
		
		//(5) 
		if(!checkCondition(n, r))
			return false;
		
		//(6) we know its not composite by this point
		return true;
	}
	
	// Determine if n is a power of an integer
	private static boolean isPowerOfInteger(BigInteger n)
	{
		// Iterate through possible exponents
		BigInteger maxB = BigInteger.valueOf(n.bitLength());
		for(BigInteger b = BigInteger.valueOf(2); b.compareTo(maxB) < 0; b = b.add(BigInteger.ONE))
		{
			BigInteger lo = BigInteger.valueOf(2);
			BigInteger hi = n.add(BigInteger.ZERO);
			
			// Binary search for value of a (the base)
			while(lo.compareTo(hi) < 0)
			{
				BigInteger a = hi.add(lo).divide(BigInteger.valueOf(2));
				
				// Calculate a^b
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
	
	// Find the smallest r such that or(n) > (log n)^2
	private static long multOrder(BigInteger n)
	{
		// Calculate the maximum power of n that shouldn't be 1 (mod r)
		long maxK = (long)Math.floor(log(n, 2) * log(n, 2));
		
		// Iterate through values of r
		for(long r = 2; ; r++)
		{
			BigInteger a = BigInteger.ONE;
			boolean valid = true;
			
			// Compute powers of n (mod r)
			for(long k = 1; k <= maxK; k++)
			{
				a = a.multiply(n).mod(BigInteger.valueOf(r));
				
				// Check if the value if 1
				if(a.compareTo(BigInteger.ONE) <= 0)
				{
					valid = false;
					break;
				}
			}
			
			// Check if we found a valid value for r
			if(valid)
			{
				return r;
			}
		}
	}
	
	//see if 1 < gcd(a,n) < n for some a<=r
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
	
	//the main part of AKS. check is (x+a)^n == x + a^n (mod n)
	//if the above condition does not hold, we know that n is composite!
	private static boolean checkCondition(BigInteger n, long r)
	{
		//upper bound of for loop that we need to check to in order to know that there is
		//no chance of a false positive prime
		long maxA = (long)Math.floor(Math.sqrt(totient(r)) * log(n, 2));
		for(long a = 1; a <= maxA; a++)
		{
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
	
	//quick and dirty logarithm function for BigInteger
	private static double log(BigInteger a, int b)
	{
		String s = a.toString();
		int l = Math.min(18, s.length());
		double d = Double.parseDouble(s.substring(0,l)) / Math.pow(10, l-1);
		return (Math.log10(d) + (s.length()-1)) / Math.log10(b);
	}
	
	//repeated squaring recursive pow function
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
	
	//totient function that iterates through all the divisors of r
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
	
	//wrapper function for fft for BigIntegers where we convert BigInteger to BigComplex first
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
			//the base case is if the polynomial is degree 1
			BigComplex[] f = new BigComplex[1];
			f[0] = a[0];
			
			return f;
		}
		else
		{
			BigComplex[] aEven = new BigComplex[m >> 1];
			BigComplex[] aOdd = new BigComplex[m >> 1];
			
			//divide the a's into the even and odd arrays	
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
			
			//compute FFT on each 'half'
			BigComplex[] fEven = fft(aEven, m >> 1, w.multiply(w));
			BigComplex[] fOdd = fft(aOdd, m >> 1, w.multiply(w));
			
			BigComplex[] f = new BigComplex[m];
			
			//build the answer using the two halves and the roots of unity
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
	
	//wrapper function, BigInteger->Complex before doing FFT
	public static Complex[] fft(BigInteger[] a, int m, Complex w)
	{
		Complex[] newCoef = new Complex[a.length];
		for(int i = 0; i < a.length; i++)
		{
			newCoef[i] = new Complex(a[i].doubleValue(), 0);
		}
		
		return fft(newCoef, m, w);
	}
	
	//same code, but uses Complex (which uses doubles to store the values)
	public static Complex[] fft(Complex[] a, int m, Complex w)
	{
		if(m == 1)
		{
			//the base case is if the polynomial is degree 1
			Complex[] f = new Complex[1];
			f[0] = a[0];
			
			return f;
		}
		else
		{
			Complex[] aEven = new Complex[m >> 1];
			Complex[] aOdd = new Complex[m >> 1];
			
			//divide the a's into the even and odd arrays
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

			//compute FFT on each 'half'
			Complex[] fEven = fft(aEven, m >> 1, w.multiply(w));
			Complex[] fOdd = fft(aOdd, m >> 1, w.multiply(w));
			
			Complex[] f = new Complex[m];
			
			//build the answer using the two halves and the roots of unity
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
		
		// Constructor for a polynomial mod x^r - 1, n
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
		
		// Constructor for a polynomial of the form x + a (mod x^r - 1, n)
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
		
		// Determine if two polynomials are equal
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
		
		//fast exponentiation of polynomials
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
		
		//multiply wrapper that uses a different multiply function based on the MODE
		private Polynomial multiply(Polynomial p)
		{
			switch(MODE)
			{
			case 0:
				return multiplySlow(p);
			default:
				return multiplyFFT(p);
			}
		}
		
		// basic O(n^2) elementary multiplication
		private Polynomial multiplySlow(Polynomial p)
		{
			Polynomial ret = new Polynomial(r, n);
			
			// We do the mod during the multiplication
			for(int i = 0; i < coef.length; i++)
			{
				for(int j = 0; j < p.coef.length; j++)
				{
					BigInteger c = coef[i].multiply(p.coef[j]).mod(n);
					
					ret.coef[(i + j) % (int)r] = ret.coef[(i + j) % (int)r].add(c).mod(n);
				}
			}
			
			return ret;
		}
		
		// Multiply two polynomials using FFT
		// Assumes both polynomials are the same degree
		@SuppressWarnings("unused")
		private Polynomial multiplyFFT(Polynomial p)
		{
			// Expand polynomials so the number of coefficients is a power of 2
			Polynomial q = expand();
			p = p.expand();
			int m = (int) Math.max(q.r, p.r);
			
			if(MODE == 1)
			{
				// mth root of unity
				BigComplex root = new BigComplex(Math.cos(2 * Math.PI / m), Math.sin(2 * Math.PI / m));
				
				// Compute fft of both polynomials
				BigComplex[] f1 = fft(q.coef, m, root);
				BigComplex[] f2 = fft(p.coef, m, root);
				
				// Multiply result of fft for each x-value
				BigComplex[] f3 = new BigComplex[m];
				for(int i = 0; i < m; i++)
				{
					f3[i] = f1[i].multiply(f2[i]);
				}
				
				// Convert to polynomial using inverse fft
				BigComplex[] c = fft(f3, m, new BigComplex(1, 0).divide(root));
				
				// Extract real coefficients into polynomial form
				Polynomial result = new Polynomial(c.length, n);
				for(int j = 0; j < c.length; j++)
				{
					BigComplex temp = c[j].divide(new BigComplex(m, 0));
					
					result.coef[j] = temp.real.add(EPSILON).toBigInteger();
				}
				
				// Perform mod
				return result.mod((int)r);
			}
			else
			{
				// mth root of unity
				Complex root = new Complex(Math.cos(2.0 * Math.PI / (double)m), Math.sin(2.0 * Math.PI / (double)m));
				
				// Compute fft of both polynomials
				Complex[] f1 = fft(q.coef, m, root);
				Complex[] f2 = fft(p.coef, m, root);
				
				// Multiply result of fft for each x-value
				Complex[] f3 = new Complex[m];
				for(int i = 0; i < m; i++)
				{
					f3[i] = f1[i].multiply(f2[i]);
				}
				
				// Convert to polynomial using inverse fft
				Complex[] c = fft(f3, m, new Complex(1, 0).divide(root));
				
				// Extract real coefficients into polynomial form
				Polynomial result = new Polynomial(c.length, n);
				for(int j = 0; j < c.length; j++)
				{
					Complex temp = c[j].divide(new Complex(m, 0));
					
					result.coef[j] = new BigDecimal(temp.real + 0.5 + 1E-3).toBigInteger();
				}
				
				// Perform mod
				return result.mod((int)r);
			}
		}
		
		//we need to expand the polynomials to some power of 2
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
			
			//fill expanded coefficients with 0 so we dont actually change the polynomial itself
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
		
		//Compute mod x^r - 1
		public Polynomial mod(int m)
		{
			//we know that cx^a mod n = cx^(a mod n)
			Polynomial p = new Polynomial(m, n);
			for(int i = 0; i < coef.length; i++)
			{
				p.coef[i % m] = p.coef[i % m].add(coef[i]).mod(n);
			}
			
			return p;
		}
	}
	
	private static class Complex
	{
		double real;
		double imag;
		
		// Contructor for Complex numbers
		public Complex(double r, double i)
		{
			real = r;
			imag = i;
		}
		
		// Add two complex numbers
		public Complex add(Complex c)
		{
			return new Complex(real + c.real, imag + c.imag);
		}
		
		// Subtract two complex numbers
		public Complex subtract(Complex c)
		{
			return new Complex(real - c.real, imag - c.imag);
		}
		
		// Multiply two complex numbers
		public Complex multiply(Complex c)
		{
			return new Complex(real * c.real - imag * c.imag, real * c.imag + imag * c.real);
		}
		
		// Divide two complex numbers
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
		
		// Constructor for complex number using doubles
		public BigComplex(double r, double i)
		{
			real = new BigDecimal(r, PRECISION);
			imag = new BigDecimal(i, PRECISION);
		}
		
		// Constructor for complex number using BigDecimals
		public BigComplex(BigDecimal r, BigDecimal i)
		{
			real = r;
			imag = i;
		}
		
		// Add two complex numbers
		public BigComplex add(BigComplex c)
		{
			return new BigComplex(real.add(c.real), imag.add(c.imag));
		}
		
		// Subtract two complex numbers
		public BigComplex subtract(BigComplex c)
		{
			return new BigComplex(real.subtract(c.real), imag.subtract(c.imag));
		}
		
		// Multiply two complex numbers
		public BigComplex multiply(BigComplex c)
		{
			return new BigComplex(real.multiply(c.real).subtract(imag.multiply(c.imag)),
					              real.multiply(c.imag).add(imag.multiply(c.real)));
		}
		
		// Divide two complex numbers
		public BigComplex divide(BigComplex c)
		{
			BigDecimal denom = c.real.multiply(c.real).add(c.imag.multiply(c.imag));
			return new BigComplex(real.multiply(c.real).add(imag.multiply(c.imag)).divide(denom, RoundingMode.HALF_EVEN),
								  imag.multiply(c.real).subtract(real.multiply(c.imag)).divide(denom, RoundingMode.HALF_EVEN));
		}
	}
	
	// Main method for testing
	public static void main(String[] args) throws IOException
	{
		Scanner in = new Scanner(System.in);
		
		while(true)
		{
			System.out.print("Enter a number >= 2 that you would like to test for primality (or 0 to exit): ");
			
			String num = in.next();
			
			if(num.equals("0"))
			{
				break;
			}
			
			boolean result = isPrime(new BigInteger(num));
			
			System.out.print(num + " is ");
			if(result)
			{
				System.out.println("prime.");
			}
			else
			{
				System.out.println("not prime.");
			}
			System.out.println();
		}
	}
}
