#modelMatrixFunction
#takes a matrix X, and the number of a column for the modifier
#note that this is not going to be the final funciton, however I need to test my understanding of 
#pushing, pulling, commiting, and stuff to gitHub 
#so hopefully I will be able to upload this my github branch
#This function takes in a data frame, and a column name 
#returns a model matrix without the specified column or intercept
modelMatrixFunction <- function(X, modifier_name){
  X_full <- X[, !(names(X) %in% modifier_name), drop = FALSE]
  modelMat <- model.matrix(~., data = X_full)
  modelMat[, -1]
}