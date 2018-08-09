library("FNN")

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("At least two arguments (k-neighbor value and kdTreeMetrics file)  must be supplied", call.=FALSE)
}


metrics<- read.table(args[2], header=TRUE)
#metrics<- read.table("ALL_kdTreeMetrics.txt", header=TRUE)

# Get k-nearest neighbors for each sample ex: k.param <- 32
k.param <- as.numeric(args[1])

#get output path
outdir <- args[3]

#treat all exomes in batch mode (from scratch)
batch <- args[4]
 


for (i in 2:ncol(metrics)) {
	mini <- min(metrics[,i])
	maxi <- max(metrics[,i])
	metrics[,i] <- apply(metrics, 1, function(row) { 
                if (maxi == 1) {
                  row[[i]] <- (as.numeric(row[[i]]) - mini)
                }
                else {
                  row[[i]] <- (as.numeric(row[[i]]) - mini) / (maxi - mini)
                }
              } 
	)
}



knns <- get.knn(metrics[,c(seq(2,ncol(metrics)))],k=k.param,algorithm="kd_tree")


# Generate a single file for each sample listing its k-nearest neighbor sample IDs
if (is.na(batch)){
	for (i in 1:1) {
		fname <- paste(outdir,metrics$SAMPLE[i], ".", k.param, "nns.txt", sep="")
		nn.sampleids <- metrics$SAMPLE[ knns$nn.index[i,] ]
		write.table(nn.sampleids, fname, quote=F, row.names=F, col.names=F)
	}
}else{
	for (i in 1:nrow(metrics)) {
		fname <- paste(outdir,metrics$SAMPLE[i], ".", k.param, "nns.txt", sep="")
		nn.sampleids <- metrics$SAMPLE[ knns$nn.index[i,] ]
		write.table(nn.sampleids, fname, quote=F, row.names=F, col.names=F)
	}
}


# To check how well each sample's kNNs fit, compute the distance to its kNN cluster mean
metrics$DistanceToClusterMean <- sapply(1:nrow(metrics),  function(x) {
	this.knns <- knns$nn.index[x,];
	center <- colMeans(metrics[this.knns, 2:ncol(metrics)]);
	return(as.numeric(dist(rbind(as.numeric(metrics[x, 2:ncol(metrics)]), as.numeric(center)))))
	}
)

# Plot distance distribution
#plot(ecdf(metrics$DistanceToClusterMean))


