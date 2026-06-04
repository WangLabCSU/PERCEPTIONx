# Opposite of standard %in% function
'%!in%' <- function(x,y)!('%in%'(x,y))


# Compute a colMax
colMax <- function (colData) {
  apply(colData, MARGIN=c(2), max)
}


# Compute a colMedian
colMedian <- function (colData) {
  apply(colData, MARGIN=c(2), median)
}


# Compute a colMin
colMin <- function (colData) {
  apply(colData, MARGIN=c(2), min)
}


# Compute a RowMax
rowMax <- function (colData) {
  apply(colData, MARGIN=c(1), max)
}


# Compute a RowMin
rowMin <- function (colData) {
  apply(colData, MARGIN=c(1), min)
}


# <!-- Capping a matrix to value 1 -->
capping_at1 <- function(x){
  x[x>1]=1
  x
}

capping_at0 <- function(x){
  x[x<0]=0
  x
}


# RowMeans functions considering the boundary case where mat only has one row
rowMeans_if_one_row <- function(mat){
  if(ncol(mat)>1){
    return(rowMeans(mat))
  } else {
    return(mat)
  }
}


# Count the number of NAs in each row of a matrix
count_row_NAs <- function(df){
  apply(df, 1, function(x) sum(is.na(x)))
}


# Handing a error when a function is run and returning NA in case of error
# instead of stopping the task.
err_handle <- function(x){ tryCatch(x, error=function(e){NA}) }


# Custom Head which only first returns 5 rows and columns
myhead <- function(mat){
  mat[1:min(5, nrow(mat)), 1:min(5, ncol(mat))]
}


# Subsetting a set of columns
colSubset <- function(mat, column_Names){
  mat[,na.omit(match(column_Names, colnames(mat)))]
}


# Subsetting a set of rows
rowSubset <- function(mat, row_Names){
  mat[na.omit(match(row_Names, rownames(mat))),]
}


# Taking vector as an input, it converts the factor vector into numeric
factor2numeric <- function(x){
  as.numeric(as.character(x))
}


# Strip all non-char and non-numeric and make lower case
# this is primarily to facilitate inconsistent naming comparison (eg. drugNames)
stripall2match <- function(x){
  tolower(gsub('[^A-z0-9]', '', x) )
}


# Split a string and return the item of interest (retreving_onject_id)
strsplit_customv0 <- function(infunc_list=pred_viab$cellLines_mapping$cellLine_ID,
                              infunc_split_by='_',
                              retreving_onject_id=1){
  sapply(strsplit(infunc_list, split = infunc_split_by), function(x) x[retreving_onject_id])
}

