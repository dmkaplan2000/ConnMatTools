library(ConnMatTools)
data(chile.loco)

num <- prod(dim(chile.loco)) / sum(chile.loco)
betas <- betasVectorDefault(n=num,steps=20)
chile.loco.split <- optimalSplitConnMat(chile.loco,normalize.cols=FALSE,
                                        betas=betas)

pops <- subpopsVectorToList(chile.loco.split$subpops[,chile.loco.split$best.splits$"3"$index])

reduce.loco <- reducedConnMat(pops,chile.loco)

sr <- selfRecruitment(reduce.loco)
lr <- localRetention(reduce.loco)
rlr <- relativeLocalRetention(reduce.loco)
