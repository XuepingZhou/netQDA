softmax <- function(data){
  data <- as.data.frame(exp(data))
  data$sum <- apply(data, 1, sum)
  data$x1 <- data[, 1] /  data$sum
  data$x2 <- data[, 2] /  data$sum

  return(data)
}

