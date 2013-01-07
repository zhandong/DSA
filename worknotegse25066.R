gse25066 = extractgse("GSE25066")



ERMarkers = c("205225_at", "205009_at","205862_at")
ERpatients = rownames(gse25066$featureData[gse25066$featureData[,4]=="er_status_ihc: P",])


splist=list()
splist[[1]]=ERMarkers
ERW=weight_estimate_missing(2^gse25066$expData[,ERpatients],splist)

ERW = t(t(ERW)/rowSums(t(ERW)))

est_ER  =log2(deconvoltion_general(2^gse25066$expData[,ERpatients],t(ERW), method="QP_LM"))

pred2=est_ER[,1]-est_ER[,2]



pred2 = prediction(pred2, reflabel)
perf2 = performance(pred2,"tpr","fpr")
plot(perf2,col="red")