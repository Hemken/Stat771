Admit <- UCBAdmissions
# Admit[1,,]
# Admit[2,,]
# Admit[,1,]
# Admit[,2,]

df <- matrix(NA, ncol=4, nrow=24)
l <- 0
for (i in 1:2) {
  admit <- 2-i
  for (j in 1:2) {
    female <- j - 1
    for (k in 1:6) {
      dept <- k
      l <- l + 1
      df[l,] <- c(admit,female,dept,Admit[i,j,k])
      
    }
  }
}

df <- as.data.frame(df)
names(df) <- c("admit", "female", "dept", "count")

df$dept <- factor(df$dept, labels=LETTERS[1:6])
# df$deptA <- as.integer(df$dept=="A")
# df$deptB <- as.integer(df$dept=="B")
# df$deptC <- as.integer(df$dept=="C")
# df$deptD <- as.integer(df$dept=="D")
# df$deptE <- as.integer(df$dept=="E")
# df$deptF <- as.integer(df$dept=="F")
# df$Cons <- 1

fit <- glm(admit~female+dept, weight=count, data=df, family="binomial")
summary(fit)

X <- model.matrix(~female+dept, data=df)
X

B <- rep(1,7)
Y <- df$admit
counts <- df$count

admitll <- function(B,Y,X,counts){
  sum(counts*(Y*(X%*%B)- log(1 + exp(X%*%B))))
}

admitll(B,Y,X,counts)

admitll(coef(fit),Y,X,counts)
admitll(coef(fit)+.01,Y,X,counts)
admitll(coef(fit)-.01,Y,X,counts)

admitgr <- function(B,Y,X,counts){
  J <- matrix(NA, ncol=1,nrow=24)
  for (i in 1:24){
    J[i,] <- (Y[i]-1/(1+exp(-X[i,]%*%B)))
  }
  return(t(X)%*%J)
}

alpha <- 0.001
B <- rep(1,7)
#B <- coef(fit)
B.inc <- admitgr(B,Y,X,counts)
B <- B + alpha*B.inc/norm(as.matrix(B.inc))
B

admitgr(coef(fit),Y,X,counts)

