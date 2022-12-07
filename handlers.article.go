package main

import (
	"fmt"
	"math/big"
	"net/http"
	"strconv"

	"github.com/gin-gonic/gin" // for the web page
)

func showIndexPage(c *gin.Context) {
	// Call the render function with the name of the template to render
	render(c, gin.H{
		"title": "AKS Algorithm",
	}, "index.html")

}

func showAnswer(c *gin.Context) {
	numString := c.PostForm("number")
	print(numString)
	n, err := strconv.Atoi(numString)
	if err != nil {
		// issue w string -> int
		panic(err)
	}
	//println(n)
	if aks(*big.NewInt(int64(n))) { // if prime
		render(c, gin.H{
			"title": "Prime Found"}, "index.html")
		print("Number: ", n, " is Prime") // for the "server" to check

	} else { // if non-prime
		render(c, gin.H{
			"title": "Prime NOT Found"}, "index.html")
		print("Number: ", n, " is Composite") // for the "server" to check
	}

}

// Render one of HTML, JSON or CSV based on the 'Accept' header of the request
// If the header doesn't specify this, HTML is rendered, provided that
// the template name is present
func render(c *gin.Context, data gin.H, templateName string) {
	// for database usage if needed
	switch c.Request.Header.Get("Accept") {
	case "application/json":
		// Respond with JSON
		c.JSON(http.StatusOK, data["payload"])
	case "application/xml":
		// Respond with XML
		c.XML(http.StatusOK, data["payload"])
	default:
		// Respond with HTML
		c.HTML(http.StatusOK, templateName, data)
	}

}

func aks(n big.Int) bool {

	var jobs int = 10
	// Set jobs to the number of goroutines to use when testing n (proccessing stuff)

	r := AKSMod(&n)
	M := AKSUpperBound(&n, r)
	a := GetAKSWitness(&n, r, &big.Int{}, M, jobs)

	if a != nil {
		// n is composite
		return false
	} else {
		// n is prime
		return true
	}
}

func isAKSWitness(n, a big.Int, tmp1, tmp2, tmp3 *bigIntPoly) bool {
	// Left-hand side: (X + a)^n mod (n, X^r - 1)
	tmp1.Set(a, *big.NewInt(1), n)
	tmp1.Pow(n, tmp2, tmp3)

	// Right-hand side: (X^n + a) mod (n, X^r - 1)
	tmp2.Set(a, n, n)

	isWitness := !tmp1.Eq(tmp2)
	return isWitness
}

// Returns the first AKS witness of n with the parameters r and M, or
// nil if there isn't one
func getFirstAKSWitness(n, r, M *big.Int) *big.Int {
	tmp1 := newBigIntPoly(*n, *r)
	tmp2 := newBigIntPoly(*n, *r)
	tmp3 := newBigIntPoly(*n, *r)

	for a := big.NewInt(1); a.Cmp(M) < 0; a.Add(a, big.NewInt(1)) {
		fmt.Printf("Testing %d (M = %d)...\n", a, M)
		isWitness := isAKSWitness(*n, *a, tmp1, tmp2, tmp3)
		if isWitness {
			return a
		}
	}
	return nil
}

// Holds the result of an AKS witness test
type witnessResult struct {
	a         *big.Int
	isWitness bool
}

// Tests all numbers received on numberCh if they are witnesses of n with parameter r
func testAKSWitnesses(
	n, r *big.Int,
	numberCh chan *big.Int,
	resultCh chan witnessResult) {
	tmp1 := newBigIntPoly(*n, *r)
	tmp2 := newBigIntPoly(*n, *r)
	tmp3 := newBigIntPoly(*n, *r)

	for a := range numberCh {
		fmt.Printf("Testing %v...\n", a)
		isWitness := isAKSWitness(*n, *a, tmp1, tmp2, tmp3)
		fmt.Printf("Finished testing %v (isWitness=%t)\n",
			a, isWitness)
		resultCh <- witnessResult{a, isWitness}
	}
}
func GetAKSWitness(
	n, r, start, end *big.Int,
	maxOutstanding int) *big.Int {
	numberCh := make(chan *big.Int, maxOutstanding)
	defer close(numberCh)
	resultCh := make(chan witnessResult, maxOutstanding)
	for i := 0; i < maxOutstanding; i++ {
		go testAKSWitnesses(n, r, numberCh, resultCh)
	}

	// Send off all numbers for testing (counted by i), draining
	// any results that come in (counted by j) while we're doing
	// so
	var i, j big.Int
	i.Set(start)
	j.Set(start)
	logResult := func(result witnessResult) {
		fmt.Printf("%v isWitness=%t\n", result.a, result.isWitness)
	}
	for i.Cmp(end) < 0 {
		select {
		case result := <-resultCh:
			j.Add(&j, big.NewInt(1))
			logResult(result)
			if result.isWitness {
				return result.a
			}
		default:
			var a big.Int
			a.Set(&i)
			numberCh <- &a
			i.Add(&i, big.NewInt(1))
		}
	}

	// Get rid of anything else while working
	for j.Cmp(end) < 0 {
		result := <-resultCh
		j.Add(&j, big.NewInt(1))
		logResult(result)
		if result.isWitness {
			return result.a
		}
	}

	return nil
}

// Upper bound for r s.t. o_r(n) > ceil(lg(n))^2 that is polylog in n
func AKSModUpperBound(n *big.Int) *big.Int {
	two := big.NewInt(2)
	three := big.NewInt(3)
	five := big.NewInt(5)
	eight := big.NewInt(8)

	//max(ceil(lg(n))^5, 3).
	ceilLgN := big.NewInt(int64(n.BitLen()))
	rUpperBound := &big.Int{}
	rUpperBound.Exp(ceilLgN, five, nil)
	rUpperBound = max(rUpperBound, three)

	var nMod8 big.Int
	nMod8.Mod(n, eight)
	if (nMod8.Cmp(three) == 0) || (nMod8.Cmp(five) == 0) {
		//8*ceil(lg(n))^2.
		var rUpperBound2 big.Int
		rUpperBound2.Exp(ceilLgN, two, nil)
		rUpperBound2.Mul(&rUpperBound2, eight)
		rUpperBound = min(rUpperBound, &rUpperBound2)
	}
	return rUpperBound
}

// Returns the least r s.t. o_r(n) > ceil(lg(n))^2 >= ceil(lg(n)^2)
func AKSMod(n *big.Int) *big.Int {
	one := big.NewInt(1)
	two := big.NewInt(2)

	ceilingMumbo := big.NewInt(int64(n.BitLen()))
	ceilingMumbo.Mul(ceilingMumbo, ceilingMumbo)
	var r big.Int
	r.Add(ceilingMumbo, two)
	rUpperBound := AKSModUpperBound(n)
	for ; r.Cmp(rUpperBound) < 0; r.Add(&r, one) {
		var gcd big.Int
		gcd.GCD(nil, nil, n, &r)
		if gcd.Cmp(one) != 0 {
			continue
		}
		o := calcMultOrder(n, &r)
		if o.Cmp(ceilingMumbo) > 0 {
			return &r
		}
	}

	panic("Could not find AKS modulus")
}

// Returns floor(sqrt(Phi(r)))*ceil(lg(n)) + 1 > floor(sqrt(Phi(r)))* lg(n)
func AKSUpperBound(n, r *big.Int) *big.Int {
	one := big.NewInt(1)
	two := big.NewInt(2)

	M := calculateEulerPhi(r)
	M = floorRoot(M, two)

	M.Mul(M, big.NewInt(int64(n.BitLen())))
	M.Add(M, one)
	return M
}

// Returns the first factor of n less than M
func GetFirstFactorBelow(n, M *big.Int) *big.Int {
	var factor *big.Int
	var mMinusOne big.Int
	mMinusOne.Sub(M, big.NewInt(1))
	trialDivide(n, func(q, e *big.Int) bool {
		if q.Cmp(M) < 0 && q.Cmp(n) < 0 {
			factor = q
		}
		return false
	}, &mMinusOne)
	return factor
}

// Returns the smaller of x and y- pointer magic
func min(x, y *big.Int) *big.Int {
	if x.Cmp(y) < 0 {
		return x
	}
	return y
}

// Returns the larger of x and y - pointer magic
func max(x, y *big.Int) *big.Int {
	if x.Cmp(y) > 0 {
		return x
	}
	return y
}

// Returns the greatest number y such that y^k <= x. x must be non-negative and k must be positive
func floorRoot(x, k *big.Int) *big.Int {
	if x.Sign() < 0 {
		panic("negative radicand")
	}
	if k.Sign() <= 0 {
		panic("non-negative index")
	}
	if x.Sign() == 0 {
		return &big.Int{}
	}
	one := big.NewInt(1)
	var kMinusOne big.Int
	kMinusOne.Sub(k, one)

	// Calculate p = ceil((floor(lg(x)) + 1)/k)
	var p, r big.Int
	p.DivMod(big.NewInt(int64(x.BitLen())), k, &r)
	if r.Sign() > 0 {
		p.Add(&p, one)
	}

	y := &big.Int{}
	y.Exp(big.NewInt(2), &p, nil)
	for y.Cmp(one) > 0 {
		// z = floor(((k-1)y + floor(x/y^{k-1}))/k)
		var z1 big.Int
		z1.Mul(&kMinusOne, y)

		var z2 big.Int
		var yPowKMinusOne big.Int
		yPowKMinusOne.Exp(y, &kMinusOne, nil)
		z2.Div(x, &yPowKMinusOne)

		var z big.Int
		z.Add(&z1, &z2)
		z.Div(&z, k)

		if z.Cmp(y) >= 0 {
			return y
		}
		y = &z
	}
	return one
}

// Assuming p is prime, calculates and returns Phi(p^k) quickly
func calcEulerwPrime(p, k *big.Int) *big.Int {
	var pMinusOne, kMinusOne big.Int
	pMinusOne.Sub(p, big.NewInt(1))
	kMinusOne.Sub(k, big.NewInt(1))
	var phi big.Int
	phi.Exp(p, &kMinusOne, nil)
	phi.Mul(&phi, &pMinusOne)
	return &phi
}

// A factorFunction takes a prime and its multiplicity and returns
// whether or not to continue trying to find more factors.
type factorFunction func(p, m *big.Int) bool

// Does trial division to find factors of n and passes them to the
// given factorFunction until it indicates otherwise. If upperBound is
// not nil, only factors less than or equal to it will be tried
func trialDivide(n *big.Int, factorFn factorFunction, upperBound *big.Int) {
	one := big.NewInt(1)
	two := big.NewInt(2)
	three := big.NewInt(3)
	four := big.NewInt(4)
	five := big.NewInt(5)
	six := big.NewInt(6)
	seven := big.NewInt(7)
	eleven := big.NewInt(11)

	if n.Sign() < 0 {
		panic("negative n")
	}
	if n.Sign() == 0 {
		return
	}

	if upperBound == nil {
		upperBound = floorRoot(n, two)
	}

	t := &big.Int{}
	t.Set(n)

	// Factors out d from t as much as possible and calls factorFn if d divides t.
	factorOut := func(d *big.Int) bool {
		var m big.Int
		for {
			var q, r big.Int
			q.QuoRem(t, d, &r)
			if r.Sign() != 0 {
				break
			}
			t = &q
			upperBound = min(upperBound, t)
			m.Add(&m, one)
		}
		if m.Sign() != 0 {
			if !factorFn(d, &m) {
				return false
			}
		}
		return true
	}

	// Try small primes first.
	if two.Cmp(upperBound) <= 0 && !factorOut(two) {
		return
	}

	if three.Cmp(upperBound) <= 0 && !factorOut(three) {
		return
	}

	if five.Cmp(upperBound) <= 0 && !factorOut(five) {
		return
	}

	if seven.Cmp(upperBound) <= 0 && !factorOut(seven) {
		return
	}

	//Run through a mod-30 wheel(it sucks >:( but it super cuts down proccessing time)
	mod30Wheel := []*big.Int{four, two, four, two, four, six, two, six}
	for i, d := 1, eleven; d.Cmp(upperBound) <= 0; {
		if !factorOut(d) {
			return
		}
		d.Add(d, mod30Wheel[i])
		i = (i + 1) % len(mod30Wheel)
	}
	if t.Cmp(one) != 0 {
		factorFn(t, one)
	}
}

// Assuming that p is prime and a and p^k are coprime, returns the smallest power e of a such that a^e = 1 (mod p^k).
func calcMultOrderPrime(a, p, k *big.Int) *big.Int {
	var n big.Int
	n.Exp(p, k, nil)
	t := calcEulerwPrime(p, k)

	o := big.NewInt(1)
	one := big.NewInt(1)
	processPrimeFactor := func(q, e *big.Int) bool {
		// Calculate x = a^(t/q^e) (mod n).
		var x big.Int
		x.Exp(q, e, nil)
		x.Div(t, &x)
		x.Exp(a, &x, &n)
		for x.Cmp(one) != 0 {
			o.Mul(o, q)
			x.Exp(&x, q, &n)
		}
		return true
	}

	if k.Cmp(one) > 0 {
		var kMinusOne big.Int
		kMinusOne.Sub(k, one)
		processPrimeFactor(p, &kMinusOne)
	}

	var pMinusOne big.Int
	pMinusOne.Sub(p, one)
	trialDivide(&pMinusOne, processPrimeFactor, nil)
	return o
}

// Assuming that a and n are coprime, returns the smallest power e s.t. a^e = 1 (mod n).
func calcMultOrder(a, n *big.Int) *big.Int {
	o := big.NewInt(1)
	trialDivide(n, func(q, e *big.Int) bool {
		oq := calcMultOrderPrime(a, q, e)
		// Set o to lcm(o, oq).
		var gcd big.Int
		gcd.GCD(nil, nil, o, oq)
		o.Div(o, &gcd)
		o.Mul(o, oq)
		return true
	}, nil)
	return o
}

// Calculate Phi(n) by factorizing it.
func calculateEulerPhi(n *big.Int) *big.Int {
	phi := big.NewInt(1)
	trialDivide(n, func(q, e *big.Int) bool {
		phi.Mul(phi, calcEulerwPrime(q, e))
		return true
	}, nil)
	return phi
}

// A bigIntPoly represents a polynomial with big.Int coefficients mod
// some (N, X^R - 1).
//
// The zero value for a bigIntPoly represents the zero polynomial.
type bigIntPoly struct {
	R int
	// k is the number of big.Words required to hold a coefficient
	// in calculations without overflowing.
	k int
	// If p(x) is the polynomial as a function, phi is
	// p(2^{k*bitsize(big.Word)}). Since a coefficient fits into k
	// big.Words, this is a lossless transformation; that is, one
	// can recover all coefficients of p(x) from phi.
	//
	// phi is set to have capacity for the largest possible
	// (intermediate) polynomial. No assumptions can be made about
	// the bytes in the unused capacity except for that the unused
	// bytes for the leading coefficient (if any) is guaranteed to
	// be zeroed out.
	phi big.Int
}

// Only polynomials built with the same value of N and R may be used
// together in one of the functions below.

// Builds a new bigIntPoly representing the zero polynomial
// mod (N, X^R - 1). R must fit into an int.
func newBigIntPoly(N, R big.Int) *bigIntPoly {
	// A coefficient can be up to R*(N - 1)^2 in intermediate
	// calculations.
	var maxCoefficient big.Int
	maxCoefficient.Sub(&N, big.NewInt(1))
	maxCoefficient.Mul(&maxCoefficient, &maxCoefficient)
	maxCoefficient.Mul(&maxCoefficient, &R)

	var phi big.Int
	rInt := int(R.Int64())
	k := len(maxCoefficient.Bits())
	// Up to 2*R coefficients may be needed in intermediate
	// calculations.
	maxWordCount := 2 * rInt * k
	phi.SetBits(make([]big.Word, maxWordCount))
	return &bigIntPoly{rInt, k, phi}
}

// Returns 1 + the degree of this polynomial, or 0 if the polynomial
// is the zero polynomial.
func (p *bigIntPoly) getCoefficientCount() int {
	l := len(p.phi.Bits())
	if l == 0 {
		return 0
	}
	coefficientCount := l / p.k
	if l%p.k != 0 {
		coefficientCount++
	}
	return coefficientCount
}

// Sets the coefficient count to the given number, which must be at
// most p.R.
func (p *bigIntPoly) setCoefficientCount(coefficientCount int) {
	p.phi.SetBits(p.phi.Bits()[0 : coefficientCount*p.k])
}

// Returns the ith coefficient of this polynomial. i must be less than
// p.getCoefficientCount().
func (p *bigIntPoly) getCoefficient(i int) big.Int {
	var c big.Int
	start := i * p.k
	end := (i + 1) * p.k
	// Since the unused data for the leading coefficient is
	// guaranteed to be zeroed out, this is okay.
	c.SetBits(p.phi.Bits()[start:end])
	return c
}

// Clears the unused bytes of the given coefficient.
func (p *bigIntPoly) commitCoefficient(c big.Int) {
	cBits := c.Bits()
	unusedBits := cBits[len(cBits):p.k]
	for j := 0; j < len(unusedBits); j++ {
		unusedBits[j] = 0
	}
}

// Sets p to X^k + a mod (N, X^R - 1).
func (p *bigIntPoly) Set(a, k, N big.Int) {
	c0 := p.getCoefficient(0)
	c0.Mod(&a, &N)
	p.commitCoefficient(c0)

	R := big.NewInt(int64(p.R))
	var kModRBig big.Int
	kModRBig.Mod(&k, R)
	kModR := int(kModRBig.Int64())

	for i := 1; i <= kModR; i++ {
		c := p.getCoefficient(i)
		c.Set(&big.Int{})
		p.commitCoefficient(c)
	}

	cKModR := p.getCoefficient(kModR)
	cKModR.Set(big.NewInt(1))
	p.commitCoefficient(cKModR)

	p.setCoefficientCount(kModR + 1)
}

// Returns whether p has the same coefficients as q
func (p *bigIntPoly) Eq(q *bigIntPoly) bool {
	return p.phi.Cmp(&q.phi) == 0
}

// Sets p to the product of p and q mod (N, X^R - 1). Assumes R >=2
func (p *bigIntPoly) mul(q *bigIntPoly, N big.Int, tmp *bigIntPoly) {
	tmp.phi.Mul(&p.phi, &q.phi)
	p.phi, tmp.phi = tmp.phi, p.phi

	// Mod p by X^R - 1.
	mid := p.R * p.k
	pBits := p.phi.Bits()
	if len(pBits) > mid {
		var lo, hi big.Int
		lo.SetBits(pBits[:mid])
		hi.SetBits(pBits[mid:])
		p.phi.Add(&lo, &hi)
		pBits = p.phi.Bits()
	}

	// Clear the unused bits of the leading coefficient if necessary
	if len(pBits)%p.k != 0 {
		start := len(pBits)
		end := start + p.k - start%p.k
		unusedBits := pBits[start:end]
		for i := 0; i < len(unusedBits); i++ {
			unusedBits[i] = 0
		}
	}
	// Commit the leading coefficient before we access it
	oldCoefficientCount := p.getCoefficientCount()
	if oldCoefficientCount > 0 {
		p.commitCoefficient(p.getCoefficient(oldCoefficientCount - 1))
	}

	// Mod p by N
	newCoefficientCount := 0
	tmp2 := tmp.getCoefficient(0)
	tmp3 := tmp.getCoefficient(1)
	for i := 0; i < oldCoefficientCount; i++ {
		c := p.getCoefficient(i)
		if c.Cmp(&N) >= 0 {
			// Mod c by N. Use big.Int.QuoRem() instead of
			// big.Int.Mod() since the latter allocates an
			// extra big.Int
			tmp2.QuoRem(&c, &N, &tmp3)
			c.Set(&tmp3)
			p.commitCoefficient(c)
		}
		if c.Sign() != 0 {
			newCoefficientCount = i + 1
		}
	}
	p.setCoefficientCount(newCoefficientCount)
}

// Sets p to p^N mod (N, X^R - 1), where R is the size of p
func (p *bigIntPoly) Pow(N big.Int, tmp1, tmp2 *bigIntPoly) {
	tmp1.phi.Set(&p.phi)

	for i := N.BitLen() - 2; i >= 0; i-- {
		tmp1.mul(tmp1, N, tmp2)
		if N.Bit(i) != 0 {
			tmp1.mul(p, N, tmp2)
		}
	}

	p.phi, tmp1.phi = tmp1.phi, p.phi
}

// fmt.Formatter implementation
func (p *bigIntPoly) Format(f fmt.State, c rune) {
	if p.phi.Sign() == 0 {
		fmt.Fprint(f, "0")
		return
	}

	// Formats coeff*x^deg
	formatNumMonomial := func(
		f fmt.State, c rune,
		coeff big.Int, deg int) {
		if coeff.Cmp(big.NewInt(1)) != 0 || deg == 0 {
			fmt.Fprint(f, &coeff)
		}
		if deg != 0 {
			fmt.Fprint(f, "x")
			if deg > 1 {
				fmt.Fprint(f, "^", deg)
			}
		}
	}

	i := p.getCoefficientCount() - 1
	formatNumMonomial(f, c, p.getCoefficient(i), i)

	for i--; i >= 0; i-- {
		coeff := p.getCoefficient(i)
		if coeff.Sign() != 0 {
			fmt.Fprint(f, " + ")
			formatNumMonomial(f, c, coeff, i)
		}
	}
}
