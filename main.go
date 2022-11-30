package main

import (
	"net/http"

	"github.com/gin-gonic/gin"
)

var router *gin.Engine

func main() {

	gin.SetMode(gin.ReleaseMode)
	router = gin.Default()
	router.LoadHTMLGlob("templates/*")

	// Handle Index
	router.GET("/", func(c *gin.Context) { // the / means on initial connection (welcome page)
		c.HTML(
			http.StatusOK,
			"index.html",
			gin.H{
				"title": "AKS Algorithm",
			},
		)
	})
	router.POST("/init", showAnswer) // at initial login to base page (localhost8080)

	router.GET("/answer", showAnswer) // for performing AKS

	router.Run()

}
