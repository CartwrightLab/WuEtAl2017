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
oo <- n > 10 & n < 150
x <- x[oo,]
n <- n[oo]
g <- ifelse(g[oo],1,2)

y <- x/n
nn <- sum(n)
xx <- colSums(x)
m2

m2<-res2.dirty[[1]]
ll<-mdm.ll(x,m2$param.p,m2$param.a)


m2<-maxModel[[6]] 
ll<-mdm.ll(x,m2$f, mdmAlphas(m2$params))
p<-ll$p.row[,1]
p1<-p[g==1]
p2<-p[g==2]

# g <- ifelse(g[oo],1,2)
pred<- prediction(p, g)
perf <- performance(pred,"tpr","fpr")
performance(pred,"auc")
plot(perf)










m2<- maxModel[[6]]

x<- dataRefDirty
source("/home/steven/Postdoc2/Project_MDM/RachelCode/dm.R")
ll<-mdm.ll(x,m2$f, mdmAlphas(m2$params))
#
ml<- maxLikelihoodTable[[6]]
pp<- m2$f

#
pr<- ml
p<- pp
prm <- apply(pr,1,max)
prm <- 0
pr <- pr-prm
pr <- t(t(exp(pr))*p)
rs <- rowSums(pr)
rss <- log(rs)+prm
all.equal(sum(rss), ll$ll)
all.equal(rss, ll$ll.row)
all.equal(pr/rs, ll$p.row)

#

eml<- t(apply(ml,1,function(x){ exp(x)*p }))
#eml<- t(t(exp(ml))*p)
peml<- eml/rowSums(eml)
all.equal(peml, ll$p.row)


#
dat<- dataFull
g <- (dat$snp == 1 & dat$snpdif == 0)[dat$callby==2]


p<-ll$p.row[,1]
p1<-p[g==1]
p2<-p[g==2]

# g <- ifelse(g[oo],1,2)
require("ROCR")
x <- cbind(dat$refs,dat$alts,dat$e1s+dat$e2s)
x <- data.matrix(x)
row.names(x) <- dat$pos
x <- x[dat$callby == 2,]

ll<-mdm.ll(x, m2$f, mdmAlphas(m2$params))
pred<- prediction(ll$p.row[,3], g)
perf <- performance(pred,"tpr","fpr")
performance(pred,"auc")
plot(perf)
