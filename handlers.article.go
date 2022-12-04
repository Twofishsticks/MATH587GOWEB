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

func fullequation(p int) []int64 { // using int 64 to make it as big as possible w/o float
	c := make([]int64, p+1) // initialize array of right size
	r := int64(1)
	for i, half := 0, p/2; i <= half; i++ { // for all possible nums
		c[i] = r
		c[p-i] = r
		r = r * int64(p-i) / int64(i+1)
	}
	for i := p - 1; i >= 0; i -= 2 { // for all the odd exponents
		c[i] = -c[i]
	}
	return c
}

func aks(p int) bool {
	c := fullequation(p)
	c[p]-- // for the (x^p-1)
	c[0]++ // for the (-1)^p
	for _, d := range c {
		if d%int64(p) != 0 { // making sure it's all divisable by num
			return false
		}
	}
	return true
}
