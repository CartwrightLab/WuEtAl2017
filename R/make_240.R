source("dm.R")

## Note the file I used.
dat <- read.delim("base_counts_hets_240_byref_flag_filtered.txt",header=T)
load("search_em2_2_88005_id240.RData")
load("search_em2_3_88015_id240.RData")
load("search_em2_4_88025_id240.RData")
load("search_em2_5_88035_id240.RData")
load("search_em2_6_88045_id240.RData")
load("search_em_2_87930_id240.RData")
load("search_em_3_87961_id240.RData")
load("search_em_4_87973_id240.RData")
load("search_em_5_87983_id240.RData")
load("search_em_6_87993_id240.RData")

x <- cbind(dat$refs,dat$alts,dat$e1s+dat$e2s)
x <- data.matrix(x)
row.names(x) <- dat$pos
## This is the Clean Data
x <- x[dat$callby == 2 & dat$snp == 1 & dat$snpdif == 0,]

## This is the Dirty Data including "g" a grouping factor
#x <- x[dat$callby == 2,]
#g <- (dat$snp == 1 & dat$snpdif == 0)[dat$callby==2]

n <- rowSums(x)
oo <- n > 10 & n < 150
x <- x[oo,]
n <- n[oo]
#g <- ifelse(g[oo],1,2)

y <- x/n
nn <- sum(n)
xx <- colSums(x)

pdf(file="qqplots_240.pdf", width=12, height=4)
par(mai=c(0.6,0.7,0.2,0.1),mfrow=c(1,3),cex.main=1.2^4,cex.lab=1.2^2)
mains = c("Reference Allele", "Alternate Allele", "Error")

#simple multinomial (no bias)
cat("\n**** Multinomial ****\n")
p <- c((xx[1]+xx[2])/2,(xx[1]+xx[2])/2,xx[3])/nn
ll <- sum(log(p)*xx)
cat(sprintf("  ll = %0.16g\n", ll))
print(p)
b <- rmultinomial(length(n), n, p)
z <- b/rowSums(b)
for(i in 1:3) {
	qqplot(z[,i],y[,i],xlim=c(0,1),ylim=c(0,1),xlab="Estimated Frequency",ylab="Observed Frequency",main=mains[i])
	abline(0,1)
	print(ad.stat.k(y[,i],z[,i]))
}

#multinomial (with ref bias)
cat("\n**** Biased Multinomial ****\n")
p <- xx/nn
ll <- sum(log(p)*xx)
cat(sprintf("  ll = %0.16g\n", ll))
print(p)
b <- rmultinomial(length(n), n, p)
z <- b/rowSums(b)
for(i in 1:3) {
	qqplot(z[,i],y[,i],xlim=c(0,1),ylim=c(0,1),xlab="Estimated Frequency",ylab="Observed Frequency",main=mains[i])
	abline(0,1)
	print(ad.stat.k(y[,i],z[,i]))
}

#dirichlet-multinomial
cat("\n**** Dirichlet-Multinomial ****\n")
m1 <- fitdm.mle(x,eps=1e-15)
cat(sprintf("  ll = %0.16g\n", m1$ll))
print(m1$param)
b <- rdm(length(n), n, m1$param)
z <- b/rowSums(b)
for(i in 1:3) {
	qqplot(z[,i],y[,i],xlim=c(0,1),ylim=c(0,1),xlab="Estimated Frequency",ylab="Observed Frequency",main=mains[i])
	abline(0,1)
	print(ad.stat.k(y[,i],z[,i]))
}

#dirichlet-multinomial mixture
cat("\n**** Dirichlet-Multinomial (2) ****\n")
# m2 <- res2[[1]]
# cat(sprintf("  ll = %0.16g\n", m2$ll))
# print(m2$param)
# b <- rmdm(length(n),n,m2$param.p,m2$param.a)
b <- rmdm(length(n),n,m2$f, m2$param[,-1])
z <- b/rowSums(b)
for(i in 1:3) {
	qqplot(z[,i],y[,i],xlim=c(0,1),ylim=c(0,1),xlab="Estimated Frequency",ylab="Observed Frequency",main=mains[i])
	abline(0,1)
# 	print(ad.stat.k(y[,i],z[,i]))
}



#dirichlet-multinomial mixture
cat("\n**** Dirichlet-Multinomial (3) ****\n")
m3 <- res3[[1]]
cat(sprintf("  ll = %0.16g\n", m3$ll))
print(m3$param)
# b <- rmdm(length(n),n,m3$param.p,m3$param.a)
b <- rmdm(length(n),n,m3$f,m3$param)
#NOTE: difference between $param $parmas
# mdmParams should be $params but $params don't have class
 class(m3$params)<- "mdmParams"

b <- rmdm(length(n),n,m3$f,m3$params)

z <- b/rowSums(b)
for(i in 1:3) {
	qqplot(z[,i],y[,i],xlim=c(0,1),ylim=c(0,1),xlab="Estimated Frequency",ylab="Observed Frequency",main=mains[i])
	abline(0,1)
# 	print(ad.stat.k(y[,i],z[,i]))
}



#dirichlet-multinomial mixture
cat("\n**** Dirichlet-Multinomial (4) ****\n")
m4 <- res4[[1]]
cat(sprintf("  ll = %0.16g\n", m4$ll))
print(m4$param)
b <- rmdm(length(n),n,m4$param.p,m4$param.a)
z <- b/rowSums(b)
for(i in 1:3) {
	qqplot(z[,i],y[,i],xlim=c(0,1),ylim=c(0,1),xlab="Estimated Frequency",ylab="Observed Frequency",main=mains[i])
	abline(0,1)
	print(ad.stat.k(y[,i],z[,i]))
}

#dirichlet-multinomial mixture
cat("\n**** Dirichlet-Multinomial (5) ****\n")
m5 <- res5[[1]]
cat(sprintf("  ll = %0.16g\n", m5$ll))
print(m5$param)
b <- rmdm(length(n),n,m5$param.p,m5$param.a)
z <- b/rowSums(b)
for(i in 1:3) {
	qqplot(z[,i],y[,i],xlim=c(0,1),ylim=c(0,1),xlab="Estimated Frequency",ylab="Observed Frequency",main=mains[i])
	abline(0,1)
	print(ad.stat.k(y[,i],z[,i]))
}

#dirichlet-multinomial mixture
cat("\n**** Dirichlet-Multinomial (6) ****\n")
m6 <- res6[[1]]
cat(sprintf("  ll = %0.16g\n", m6$ll))
print(m6$param)
b <- rmdm(length(n),n,m6$param.p,m6$param.a)
z <- b/rowSums(b)
for(i in 1:3) {
	qqplot(z[,i],y[,i],xlim=c(0,1),ylim=c(0,1),xlab="Estimated Frequency",ylab="Observed Frequency",main=mains[i])
	abline(0,1)
	print(ad.stat.k(y[,i],z[,i]))
}

graphics.off()

#save data
#save(dat,x,n,m1,m2,m3,m3,m4,m5,v2,v3,v4,v5,file="make_results.RData")

#   ll2 <- mdm.ll(x,m2z$param.p,m2z$param.a)
#   s2 <- split(as.data.frame(ll2$p.row),g)


