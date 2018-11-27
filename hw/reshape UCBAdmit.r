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
# X

B <- rep(1,7)
Y <- df$admit
counts <- df$count

admitll <- function(B,Y,X,counts){
  sum(counts*(Y*(X%*%B)- log(1 + exp(X%*%B))))
}

admitll(B,Y,X,counts)

# Check the log-likelihood around the solution
admitll(coef(fit),Y,X,counts)
logLik(fit) # check
admitll(coef(fit)+.01,Y,X,counts)
admitll(coef(fit)-.01,Y,X,counts)

admitgr <- function(B,Y,X,counts){
  J <- (Y-1/(1+exp(-X%*%B)))*counts
  return((t(X)%*%J)/sum(counts))
}

B <- rep(1,7)
admitgr(B,Y,X,counts)
# solve(admitgr(B,Y,X,counts), B)
# B <- coef(fit)
alpha <- 1
n <- 3000
i <- 0
ll <- rep(NA, n)
ngr <- rep(NA, n)
while (norm(admitgr(B,Y,X,counts), "F") > 1e-8 && i < n) {
#for (i in 1:n) {
  i <- i+1
  B.inc <- admitgr(B,Y,X,counts)
  ngr[i] <- norm(as.matrix(B.inc), "F")
  B <- B + alpha*B.inc #/norm(as.matrix(B.inc))
  ll[i] <- admitll(B,Y,X,counts)
}

B; i
norm(admitgr(B,Y,X,counts), "F")
admitll(B,Y,X,counts); logLik(fit)
plot(ll)
plot(ngr)
coef(fit)

# admitgr(coef(fit),Y,X,counts)
norm(admitgr(coef(fit),Y,X,counts), "F")

B <- rep(1,7)
H <- t(X)%*%diag(as.vector(1/(1+exp(-X%*%B))*(1-(1/(1+exp(-X%*%B)))))*counts)%*%X/sum(counts)
H

B <- rep(1,7)
n <- 25
i <- 0
ll <- rep(NA, n)
ngr <- rep(NA, n)
while (norm(admitgr(B,Y,X,counts), "F") > 1e-8 && i < n) {
  i <- i+1
  B.inc <- solve(H)%*%admitgr(B,Y,X,counts)
  ngr[i] <- norm(as.matrix(B.inc), "F")
  B <- B + alpha*B.inc
  ll[i] <- admitll(B,Y,X,counts)
}
cbind(B, coef(fit));i

plot(ll)
plot(ngr)
