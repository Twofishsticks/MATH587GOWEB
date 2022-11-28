// handlers.article.go

package main

import (
	"net/http"

	"github.com/gin-gonic/gin"
)

func showIndexPage(c *gin.Context) {
	// Call the render function with the name of the template to render
	render(c, gin.H{
		"title": "AKS Algorithm",
	}, "index.html")

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

// prime if ( x – 1 )^n – ( x^n – 1) is divisible by n, x being some value, n being input value
// true if prime, false if nonprime
func AKSalgProcess(n int) bool {
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
	var coef []int
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
