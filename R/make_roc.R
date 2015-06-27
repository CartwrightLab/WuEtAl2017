source("mdm.R")

## Note the file I used.
dat <- read.delim("base_counts_hets_878_byref_flag_filtered.txt",header=T)
load("search_em.RData")

x <- cbind(dat$refs,dat$alts,dat$e1s+dat$e2s)
x <- data.matrix(x)
row.names(x) <- dat$pos
## This is the Clean Data
#x <- x[dat$callby == 2 & dat$snp == 1 & dat$snpdif == 0,]

## This is the Dirty Data including "g" a grouping factor
x <- x[dat$callby == 2,]
g <- (dat$snp == 1 & dat$snpdif == 0)[dat$callby==2]

n <- rowSums(x)
oo <- n > 20 & n < 110
x <- x[oo,]
n <- n[oo]
g <- ifelse(g[oo],1,2)

y <- x/n
nn <- sum(n)
xx <- colSums(x)

m2<-res2.dirty[[1]]
ll<-mdm.ll(x,m2$param.p,m2$param.a)
p<-ll$p.row[,2]
p1<-p[g==1]
p2<-p[g==2]