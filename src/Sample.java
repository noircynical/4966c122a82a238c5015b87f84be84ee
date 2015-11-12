import java.io.IOException;
import java.math.BigInteger;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Random;

public class Sample {
	protected static final BigInteger
    	_0 = BigInteger.valueOf(0L),
    	_1 = BigInteger.valueOf(1L),
    	_2 = BigInteger.valueOf(2L),
    	_3 = BigInteger.valueOf(3L),
    	_4 = BigInteger.valueOf(4L),
    	_5 = BigInteger.valueOf(5L),
    	_6 = BigInteger.valueOf(6L);
	
	public static int iteration= 10;
	
	public static byte[] randSeed= null;
	public static SecureRandom rnd= null;
	public static BNParams bn= null;
	public static BNCurve E= null;
	public static BNCurve2 Et= null;
	public static BNPairing pair= null;
	
	
	public static void main(String[] args) throws IOException, NoSuchAlgorithmException {
		if(randSeed == null) randSeed = new byte[20];
        (new Random()).nextBytes(randSeed);
        for (int i = 0; i < randSeed.length; i++) randSeed[i] = (byte)i;
        if(rnd == null) rnd = new SecureRandom(randSeed);
        
        for (int i = 0; i < BNParams.validBitsRange.length; i++) {
            System.out.println("\n======== bits: " + BNParams.validBitsRange[i]);
            if(bn == null) bn = new BNParams(BNParams.validBitsRange[i]);
            System.out.println("u = " + bn.u + " (fam 1)");
            if(E == null) E = new BNCurve(bn); System.out.println("E: " + E);
            E.G.getSerializedTable();
            if(Et == null) Et = new BNCurve2(E); System.out.println("Et: "+Et);
            if(pair == null) pair = new BNPairing(Et);
            test(true);
        }
	}
	
	public static void test(boolean verbose){
		BNPoint w, x, y, z, ecZero;
        BigInteger m, n;
        int numBits = 256;
        
        long totalElapsed = -System.currentTimeMillis();
        for (int i = 0; i < iteration; i++) {
            if (verbose) System.out.print("test #" + i);
            long elapsed = -System.currentTimeMillis();
            x = E.G.randomize(rnd);
            y = E.G.randomize(rnd);
            z = E.G.randomize(rnd);
            ecZero = E.G.E.infinity;
            m = new BigInteger(numBits, rnd);
            n = new BigInteger(numBits, rnd);

// ------------------------------------------------------------------------------------------------
            if (iteration == 1) System.out.print("\nchecking cloning/comparison/pertinence");
            if (!x.equals(x)) throw new RuntimeException("Comparison failure");
            if (verbose) System.out.print(".");
            if (!x.isOnSameCurve(x)) throw new RuntimeException("Inconsistent pertinence self-comparison");
            if (verbose) System.out.print(".");
            if (!x.E.contains(x)) throw new RuntimeException("Inconsistent curve pertinence");
            if (verbose) System.out.print(".");
            if (iteration == 1) System.out.print(" done.\nchecking addition properties");
            if (!x.add(y).equals(y.add(x))) throw new RuntimeException("x + y != y + x");
            if (verbose) System.out.print(".");
            if (!x.add(ecZero).equals(x)) throw new RuntimeException("x + 0 != x");
            if (verbose) System.out.print(".");
            if (!x.add(x.negate()).isZero()) throw new RuntimeException("x + (-x) != 0");
            if (verbose) System.out.print(".");
            if (!x.add(y).add(z).equals(x.add(y.add(z)))) throw new RuntimeException("(x + y) + z != x + (y + z)");
            if (verbose) System.out.print(".");
            if (!x.negate().negate().equals(x)) throw new RuntimeException("-(-x) != x");
            if (iteration == 1) System.out.print(" done.\nchecking scalar multiplication properties");
            if (!x.multiply(BigInteger.valueOf(0L)).equals(ecZero)) throw new RuntimeException("0*x != 0");
            if (verbose) System.out.print(".");
            if (!x.multiply(BigInteger.valueOf(1L)).equals(x)) throw new RuntimeException("1*x != x");
            if (verbose) System.out.print(".");
            if (!x.multiply(BigInteger.valueOf(2L)).equals(x.twice(1))) throw new RuntimeException("2*x != twice x");
            if (verbose) System.out.print(".");
            if (!x.multiply(BigInteger.valueOf(2L)).equals(x.add(x))) throw new RuntimeException("2*x != x + x");
            if (verbose) System.out.print(".");
            if (!x.multiply(BigInteger.valueOf(-1L)).equals(x.negate())) throw new RuntimeException("(-1)*x != -x");
            if (verbose) System.out.print(".");
            if (!x.multiply(m.negate()).equals(x.negate().multiply(m))) throw new RuntimeException("(-m)*x != m*(-x)");
            if (verbose) System.out.print(".");
            if (!x.multiply(m.negate()).equals(x.multiply(m).negate())) throw new RuntimeException("(-m)*x != -(m*x)");
            if (verbose) System.out.print(".");
            if (!x.multiply(m.add(n)).equals(x.multiply(m).add(x.multiply(n)))) throw new RuntimeException("(m + n)*x != m*x + n*x");
            if (verbose) System.out.print(".");
            w = x.multiply(n).multiply(m);
            if (!w.equals(x.multiply(m).multiply(n))) throw new RuntimeException("m*(n*x) != n*(m*x)");
            if (verbose) System.out.print(".");
            if (!w.equals(x.multiply(m.multiply(n)))) throw new RuntimeException("m*(n*x) != (m*n)*x");
// ------------------------------------------------------------------------------------------------            
            
            BNField12 g, a, b, c;
            
            System.out.println("Testing Tate pairing");
            g = pair.tate(E.G, Et.Gt);
            System.out.println("g = " + g);
            System.out.println("g^n = " + g.exp(bn.n));
            if (g.isZero()) throw new RuntimeException("degeneracy error!");
            if (!g.exp(bn.n).isOne()) throw new RuntimeException("G_T order error!");
            a = pair.tate(E.G.twice(1).normalize(), Et.Gt);
            b = pair.tate(E.G, Et.Gt.twice(1).normalize());
            c = g.square();
            if(!(a.equals(b) && b.equals(c))) System.out.println("bilinear? false");
            if (!(a.equals(b) && b.equals(c)) || a.isOne()) {
                System.out.println(">>>> a = " + a);
                System.out.println(">>>> b = " + b);
                System.out.println(">>>> c = " + c);
                throw new RuntimeException("Bilinearity error!");
            }
            for (int j = 0; j < 10; j++) {
                m = new BigInteger(BNParams.validBitsRange_oneValue, rnd);
                a = pair.tate(E.G.multiply(m), Et.Gt);
                b = pair.tate(E.G, Et.Gt.multiply(m));
                c = g.exp(m);
                if(!(a.equals(b) && b.equals(c))) System.out.println("bilinear? false");
                if (!(a.equals(b) && b.equals(c)) || a.isOne()) {
                    System.out.println(">>>> a = " + a);
                    System.out.println(">>>> b = " + b);
                    System.out.println(">>>> c = " + c);
                    throw new RuntimeException("Bilinearity error!");
                }
            }
        }
        totalElapsed += System.currentTimeMillis();
        System.out.println(" OK; all " + iteration + " tests done in " + (float)totalElapsed/1000 + " s.");
	}
}
