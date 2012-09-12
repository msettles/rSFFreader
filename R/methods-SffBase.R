

setMethod(show,
          signature=signature(object=".SffBase"),
          function(object) {
              cat("class: ", class(object), "\n", sep="")
          })

setMethod(detail,
          signature=signature(x=".SffBase"),
          function(x, ...) {
              cat("class: ", class(x), "\n", sep="")
          })
