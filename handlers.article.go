// handlers.article.go

package main

import (
	"net/http"
	"strconv"

	//"strconv"

	"github.com/gin-gonic/gin"
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

	if aks(n) { // if prime
		render(c, gin.H{
			"title": "Prime Found"}, "index.html")
		print("Prime")

	} else { // if non-prime
		render(c, gin.H{
			"title": "Prime NOT Found"}, "index.html")
		//c.HTML(http.StatusBadRequest, "login.html", gin.H{
		//"ErrorTitle":   "Login Failed",
		//"ErrorMessage": "Invalid credentials provided"})
		print("Nonprime")
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

/*
// prime if ( x – 1 )^n – ( x^n – 1) is divisible by n, x being some value, n being input value
// true if prime, false if nonprime
func AKSalgProcess(num string) bool {
	n, err := strconv.Atoi(num)
	if err != nil {
		// issue w string -> int
		panic(err)
	}
	coef := allCoef(n)
	// now check to see if ALL coef are divisable by n
	i := n
	for (i > 0) && (coef[i]%n == 0) { // go until there's no more coef to check or coef is not divisable
		i-- // check the next coef value
	}

	return i < 0 // if loop broke on non divisable, return false
}

// returns all coefficients of ( x – 1 )^n – ( x^n – 1), using pascal's triangle
func allCoef(n int) []int {
	var coef = make([]int, n+1)
	for i := 0; i < n; i++ { // for each ^n
		coef[1+i] = 1 // set pascal triangle position

		for j := i; j > 0; j-- {
			coef[j] = coef[j-1] - coef[j]
		}
		coef[0] = -coef[0] // swap the first val (same position as i++)

	}
	coef[0]++ // the "-(x^n -1)" one
	coef[n]-- // the (x-1)^n one
	return coef
}
*/

func bc(p int) []int64 {
	c := make([]int64, p+1)
	r := int64(1)
	for i, half := 0, p/2; i <= half; i++ {
		c[i] = r
		c[p-i] = r
		r = r * int64(p-i) / int64(i+1)
	}
	for i := p - 1; i >= 0; i -= 2 {
		c[i] = -c[i]
	}
	return c
}

func aks(p int) bool {
	c := bc(p)
	c[p]-- // for the (x^p-1)
	c[0]++ // for the (-1)^p
	for _, d := range c {
		if d%int64(p) != 0 { // making sure it's all divisable by num
			return false
		}
	}
	return true
}
