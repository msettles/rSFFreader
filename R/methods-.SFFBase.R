

setMethod(show,
          signature=signature(object=".SFFBase"),
          function(object) {
              cat("class: ", class(object), "\n", sep="")
          })

setMethod(detail,
          signature=signature(x=".SFFBase"),
          function(x, ...) {
              cat("class: ", class(x), "\n", sep="")
          })
