#WGCNA smallRNAs
###LFMP0

screen -S WGCNAsmall


cd /home/lfmp/Carabica/WGCNA

#####################
export ALLOW_WGCNA_THREADS=8

R
library("WGCNA")
library(edgeR)
load("DEsmallRNAseq.RData")
enableWGCNAThreads(8)

PHASdata <-read.delim("TE.andSK.Results.txt",header=T)

cpm.matrix <- cpm(analysis_matrix,normalized.lib.sizes=FALSE)
cpm.matrix <- log(cpm.matrix+1,2)
cpm.matrix <- t(cpm.matrix)
#colnames(cpm.matrix) <- #gsub("NC_039898.1","1",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039899.1","2",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039900.1","3",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039901.1","4",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039902.1","5",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039903.1","6",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039904.1","7",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039905.1","8",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039906.1","9",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039907.1","10",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039908.1","11",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039909.1","12",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039910.1","13",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039911.1","14",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039912.1","15",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039913.1","16",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039914.1","17",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039915.1","18",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039916.1","19",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039917.1","20",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039918.1","21",colnames(cpm.matrix))

#colnames(cpm.matrix) <- gsub("NC_039919.1","22",colnames(cpm.matrix))



#We first check for genes and samples with too many missing values:

gsg = goodSamplesGenes(cpm.matrix, verbose = 3);
gsg$allOK




#If the last statement returns TRUE , all genes have passed the cuts. If not, we remove the offending genes and samplesfrom the data:
if (!gsg$allOK)
{
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(cpm.matrix)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(cpm.matrix)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
cpm.matrix0 = cpm.matrix[gsg$goodSamples, gsg$goodGenes]
}

#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.

sampleTree = hclust(dist(cpm.matrix), method = "average");


# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 15, col = "red")
dev.off()


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(cpm.matrix, powerVector = powers, verbose = 5)

# Plot the results:
pdf("pick_soft_treshold.pdf", width = 12, height = 9)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


softPower = 6;

adjacency = adjacency(cpm.matrix, power = softPower);

#Topological Overlap Matrix (TOM)
#To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

#Clustering using TOM

#We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes. Note that we use the function hclust that provides a much faster hierarchical clustering routine than the standard hclust function.

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf("Gene_clustering_on_TOM-based_dissimilarity.pdf", width = 12, height = 9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04);
dev.off()


#The clustering dendrogram plotted by the last command is shown in Figure 2. In the clustering tree (dendrogram), each leaf, that is a short vertical line, corresponds to a gene. Branches of the dendrogram group together densely interconnected, highly co-expressed genes. Module identification amounts to the identification of individual branches (”cutting the branches off the dendrogram”). There are several methods for branch cutting; our standard method is the Dynamic Tree Cut from the package dynamicTreeCut. The next snippet of code illustrates its use.

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
#I choosed to set this parameter to 100 to reduce the total number of modules - we are reporting more than 20 here
#minModuleSize = 100;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
library(RColorBrewer)
dynamicColors = labels2colors(dynamicMods,colorSeq=brewer.pal(8,"Paired"))
table(dynamicColors)
# Plot the dendrogram and colors underneath

pdf("Gene_clustering_on_TOM-based_dissimilarity_colors.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
dev.off()

dynamicColorsNames <- dynamicColors

dynamicColorsNames[dynamicColorsNames =="#1F78B4" ] <- "Mariner"


dynamicColorsNames[dynamicColorsNames =="#33A02C" ] <- "LimeGreen"


dynamicColorsNames[dynamicColorsNames =="#A6CEE3" ] <- "BlueJeans"


dynamicColorsNames[dynamicColorsNames =="#B2DF8A" ] <- "YellowGreen"


#dynamicColorsNames[dynamicColorsNames =="#E31A1C" ] <- "RubyRed"


dynamicColorsNames[dynamicColorsNames =="#FB9A99" ] <- "LightSalmon"


#dynamicColorsNames[dynamicColorsNames =="#FDBF6F" ] <- "Peach"

plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function

moduleColors=dynamicColors
png(file = "network_heatmap_4000.png", width = 4000, height = 4000);
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes", margin = c(50,50))
dev.off()


#terrainColors
png(file = "network_heatmap_4000terrainCorlors.png", width = 4000, height = 4000);
TOMplot(plotTOM, geneTree, moduleColors,terrainColors=T,  main = "Network heatmap plot, all genes", margin = c(50,50))
dev.off()


#ColorCodes <-c("#1F78B4","#33A02C","#A6CEE3","#B2DF8A","#E31A1C","#FB9A99","#FDBF6F")
#ColorNames <- c("Mariner","LimeGreen","JeansBlue","YellowGreen","RubyRed","LightSalmon","Peach")

ColorCodes <-c("#1F78B4","#33A02C","#A6CEE3","#B2DF8A","#FB9A99")
ColorNames <- c("Mariner","LimeGreen","BlueJeans","YellowGreen","LightSalmon")


equivalenceTable <- cbind(ColorCodes,ColorNames)



###############
probes = colnames(cpm.matrix)

#PHAS21 <- PHASdata$X.Locus[PHASdata$DicerCall==21]

#PHAS21<- PHAS21[!is.na(PHAS21)]

#for (loci in PHAS21){

#  probes[grep(loci,probes)] <- gsub("PHAS","21PHAS",probes[grep(loci,probes)])
#}


#PHAS24 <- PHASdata$X.Locus[PHASdata$DicerCall==24]

#PHAS24<- PHAS24[!is.na(PHAS24)]

#for (loci in PHAS24){
#  probes[grep(loci,probes)] <- gsub("PHAS","24PHAS",probes[grep(loci,probes)])
}

#probes[grep("PHAS",probes)]

#Mariner
inModule = is.finite(match(moduleColors, "#1F78B4"));

#probes = colnames(cpm.matrix)

MarinermodProbes = probes[inModule];

write.table(MarinermodProbes,file="Mariner.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)


#LimeGreen
inModule = is.finite(match(moduleColors, "#33A02C"));

#probes = colnames(cpm.matrix)

LimeGreenmodProbes = probes[inModule];

write.table(LimeGreenmodProbes,file="LimeGreen.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#JeansBlue
inModule = is.finite(match(moduleColors, "#A6CEE3"));

#probes = colnames(cpm.matrix)

JeansBluemodProbes = probes[inModule];

write.table(JeansBluemodProbes,file="BlueJeans.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)


#YellowGreen
inModule = is.finite(match(moduleColors, "#B2DF8A"));

#probes = colnames(cpm.matrix)

YellowGreenmodProbes = probes[inModule];

write.table(YellowGreenmodProbes,file="YellowGreen.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)


##RubyRed
#inModule = is.finite(match(moduleColors, "#E31A1C"));

#probes = colnames(cpm.matrix)

#RubyRedmodProbes = probes[inModule];

#write.table(RubyRedmodProbes,file="RubyRed.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)


#LightSalmon
inModule = is.finite(match(moduleColors, "#FB9A99"));

#probes = colnames(cpm.matrix)

LightSalmonmodProbes = probes[inModule];

write.table(LightSalmonmodProbes,file="LightSalmon.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)


#Peach
#inModule = is.finite(match(moduleColors, "#FDBF6F"));

#probes = colnames(cpm.matrix)

#PeachmodProbes = probes[inModule];

#write.table(PeachmodProbes,file="Peach.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)


################################################
pastelColors <- c("#B97C7C",
"#698595",
"#9BA9CE",
"#C0C8CE",
"#CAB393",
"#BF7F58",
"#5F6D54",
"#BCDDE3"

)
################################################
waffle <- function(x, rows, cols = seq_along(x), ...) {
  xx <- rep(cols, times = x)
  lx <- length(xx)
  m <- matrix(nrow = rows, ncol = (lx %/% rows) + (lx %% rows != 0))
  m[1:length(xx)] <- xx

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  par(list(...))
  plot.new()
  o <- cbind(c(row(m)), c(col(m))) + 1
  plot.window(xlim = c(0, max(o[, 2]) + 1), ylim = c(0, max(o[, 1]) + 1),
              asp = 1, xaxs = 'i', yaxs = 'i')
  rect(o[, 2], o[, 1], o[, 2] + .85, o[, 1] + .85, col = c(m), border = NA)

  invisible(list(m = m, o = o))
}


cols <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
pdf("test waffle.pdf")
waffle(c(80, 30, 20, 10), rows = 8, cols = cols, mar = c(0,0,0,7),bg = 'cornsilk')
legend('right', legend = LETTERS[1:4], pch = 15, col = cols, pt.cex = 2, bty = 'n')
dev.off()

################################################

modules_Mariner = "#1F78B4"


inModule = is.finite(match(moduleColors, modules_Mariner));
modProbes <- probes[inModule]


#get the number of each type of small RNAs in the module

SNRNA <- length(grep("SNRNA", modProbes))

SNORNA <- length(grep("SNORNA", modProbes))

TRNA <- length(grep("TRNA", modProbes))

PHAS21 <- length(grep("21PHAS", modProbes))

PHAS24 <- length(grep("24PHAS", modProbes))

L24P_like <- length(grep("L24P", modProbes))

UNK <- length(grep("UNK", modProbes))

miRNA <- length(grep("miR", modProbes))



allTypes <- c(SNRNA, SNORNA, TRNA, PHAS21, PHAS24, L24P_like ,UNK, miRNA)
names(allTypes) <- c("sn-RNA","sno-RNA","t-RNA","21-PHAS","24-PHAS","L24P-like","UNK","mi-RNA")

#pie
pdf("All.Mariner.pie.pdf")
par(mar=c(2,2,2,2))
pie(allTypes,col=pastelColors,labels="")
dev.off()


pdf("All.Mariner.waffle.pdf")
waffle(allTypes, rows = 20, cols = pastelColors, mar = c(1,1,1,1),bg = 'white')
#legend('right', legend = LETTERS[1:4], pch = 15, col = cols, pt.cex = 2, bty = 'n')
dev.off()



adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

summary(adj[adj != 1])

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.35] = 1
adj[adj != 1] = 0


library(igraph)
(MarinerNetwork <- graph.adjacency(adj))
(MarinerNetwork <- simplify(MarinerNetwork))  # removes self-loops
(MarinerNetwork <- delete.vertices(MarinerNetwork, degree(MarinerNetwork)==0))


ColorTable <- cbind(probes,moduleColors)

V(MarinerNetwork)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(MarinerNetwork)))),2]

# remove unconnected nodes
MarinerNetwork <- delete.vertices(MarinerNetwork, degree(MarinerNetwork)==0)

length(names(V(MarinerNetwork)))

topConnectedMariner <- names(V(MarinerNetwork))

write.table(topConnectedMariner,file="topConnectedMariner.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the MarinerNetwork. It can be thought of as how critical the vertex is to the flow of information through a MarinerNetwork. Individuals with high betweenness are key bridges between different parts of a MarinerNetwork. In our measles transmission MarinerNetwork, vertices with high betweenness are those children who were central to passing on the disease to other parts of the MarinerNetwork. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the MarinerNetwork adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(MarinerNetwork, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweeness_Mariner.pdf")
hist(g.b, breaks = 80)
dev.off()

png(file = "nicely.MarinerNetwork_Mariner.png", width = 4000, height = 4000);
plot(MarinerNetwork,layout=layout_nicely(MarinerNetwork), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.MarinerNetwork_Mariner.png", width = 4000, height = 4000);
plot(MarinerNetwork,layout=layout_with_fr(MarinerNetwork), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.MarinerNetwork_Mariner.png", width = 4000, height = 4000);
plot(MarinerNetwork,layout=layout_as_tree(MarinerNetwork), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.MarinerNetwork_Mariner.png", width = 4000, height = 4000);
plot(MarinerNetwork,layout=layout_in_circle(MarinerNetwork), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the MarinerNetwork graph.
gd <- edge_density(MarinerNetwork)

#Another measure of how interconnected a MarinerNetwork is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the MarinerNetwork. The longest path length between any pair of vertices is called the diameter of the MarinerNetwork graph.

# Get the diameter of the graph g
diameter(MarinerNetwork, directed = FALSE)
get_diameter(MarinerNetwork, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(MarinerNetwork)
# Shows the path sequence between two furthest apart vertices.
get_diameter(MarinerNetwork)

# Get the average path length of the graph g
g.apl <- mean_distance(MarinerNetwork, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an important method for investigating how likely or unlikely other MarinerNetwork metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph
Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(MarinerNetwork)

gorder(MarinerNetwork)

g.random <- erdos.renyi.game(n = gorder(MarinerNetwork), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random
mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(MarinerNetwork), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthMariner.pdf")
hist(gl.apls,xlim=c(1.73,1.77),col="#1F78B4",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.003, y= 150,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,2),col="red",line=.46,cex=.9)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length is lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the MarinerNetwork.

g.b <- betweenness(MarinerNetwork, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
#pdf("histVertexBetweenessMariner.pdf")
#hist(g.b,breaks = g.b[which.max(g.b)], col ="#1F78B4", main="Betweenness Histogram",xlab="Betweenness")
#dev.off()

#Distances between vertices

#The inter-connectivity of a MarinerNetwork can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(MarinerNetwork, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

#pdf ("histDegreeMariner.pdf")
#hist(g.outd, breaks = g.outd[which.max(g.outd)]
#, col ="#1F78B4", main="Degree #Histogram",xlab="Degree")
#dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(MarinerNetwork)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our MarinerNetwork comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoMarinerNetwork <- make_ego_graph(MarinerNetwork, 1, "TRNA_11:3621510-3621582", mode=c("all"))[[1]]


png(file = "nicely.MarinerNetwork_ego_Mariner.png", width = 4000, height = 4000);
plot(egoMarinerNetwork,layout=layout_nicely(egoMarinerNetwork), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

equivalenceTable

################################################


modules_LimeGreen = "#33A02C"


inModule = is.finite(match(moduleColors, modules_LimeGreen));
modProbes <- probes[inModule]


#get the number of each type of small RNAs in the module

SNRNA <- length(grep("SNRNA", modProbes))

SNORNA <- length(grep("SNORNA", modProbes))

TRNA <- length(grep("TRNA", modProbes))

PHAS21 <- length(grep("21PHAS", modProbes))

PHAS24 <- length(grep("24PHAS", modProbes))

L24P_like <- length(grep("L24P", modProbes))

UNK <- length(grep("UNK", modProbes))

miRNA <- length(grep("miR", modProbes))



allTypes <- c(SNRNA, SNORNA, TRNA, PHAS21, PHAS24, L24P_like ,UNK, miRNA)
names(allTypes) <- c("sn-RNA","sno-RNA","t-RNA","21-PHAS","24-PHAS","L24P-like","UNK","mi-RNA")


pdf("All.LimeGreen.pie.pdf")
par(mar=c(2,2,2,2))
pie(allTypes,col=pastelColors,labels='')
dev.off()


pdf("All.LimeGreen.waffle.pdf")
waffle(allTypes, rows = 20, cols = pastelColors, mar = c(1,1,1,1),bg = 'white')
#legend('right', legend = LETTERS[1:4], pch = 15, col = cols, pt.cex = 2, bty = 'n')
dev.off()



adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

summary(adj[adj != 1])

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.58] = 1
adj[adj != 1] = 0


library(igraph)
(LimeGreenNetwork <- graph.adjacency(adj))
(LimeGreenNetwork <- simplify(LimeGreenNetwork))  # removes self-loops
(LimeGreenNetwork <- delete.vertices(LimeGreenNetwork, degree(LimeGreenNetwork)==0))


ColorTable <- cbind(probes,moduleColors)

V(LimeGreenNetwork)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(LimeGreenNetwork)))),2]

# remove unconnected nodes
LimeGreenNetwork <- delete.vertices(LimeGreenNetwork, degree(LimeGreenNetwork)==0)

length(names(V(LimeGreenNetwork)))

topConnectedLimeGreen <- names(V(LimeGreenNetwork))

write.table(topConnectedLimeGreen,file="topConnectedLimeGreen.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the LimeGreenNetwork. It can be thought of as how critical the vertex is to the flow of information through a LimeGreenNetwork. Individuals with high betweenness are key bridges between different parts of a LimeGreenNetwork. In our measles transmission LimeGreenNetwork, vertices with high betweenness are those children who were central to passing on the disease to other parts of the LimeGreenNetwork. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the LimeGreenNetwork adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(LimeGreenNetwork, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweeness_LimeGreen.pdf")
hist(g.b, breaks = 80)
dev.off()

png(file = "nicely.LimeGreenNetwork_LimeGreen.png", width = 4000, height = 4000);
plot(LimeGreenNetwork,layout=layout_nicely(LimeGreenNetwork), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.LimeGreenNetwork_LimeGreen.png", width = 4000, height = 4000);
plot(LimeGreenNetwork,layout=layout_with_fr(LimeGreenNetwork), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.LimeGreenNetwork_LimeGreen.png", width = 4000, height = 4000);
plot(LimeGreenNetwork,layout=layout_as_tree(LimeGreenNetwork), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.LimeGreenNetwork_LimeGreen.png", width = 4000, height = 4000);
plot(LimeGreenNetwork,layout=layout_in_circle(LimeGreenNetwork), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the LimeGreenNetwork graph.
gd <- edge_density(LimeGreenNetwork)

#Another measure of how interconnected a LimeGreenNetwork is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the LimeGreenNetwork. The longest path length between any pair of vertices is called the diameter of the LimeGreenNetwork graph.

# Get the diameter of the graph g
diameter(LimeGreenNetwork, directed = FALSE)
get_diameter(LimeGreenNetwork, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(LimeGreenNetwork)
# Shows the path sequence between two furthest apart vertices.
get_diameter(LimeGreenNetwork)

# Get the average path length of the graph g
g.apl <- mean_distance(LimeGreenNetwork, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an important method for investigating how likely or unlikely other LimeGreenNetwork metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph
Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(LimeGreenNetwork)

gorder(LimeGreenNetwork)

g.random <- erdos.renyi.game(n = gorder(LimeGreenNetwork), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random
mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

#Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(LimeGreenNetwork), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthLimeGreen.pdf")
hist(gl.apls,xlim=c(1.73,1.77),col="#1F78B4",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.003, y= 150,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,2),col="red",line=.46,cex=.9)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length is lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the LimeGreenNetwork.

g.b <- betweenness(LimeGreenNetwork, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
#pdf("histVertexBetweenessLimeGreen.pdf")
#hist(g.b,breaks = g.b[which.max(g.b)], col ="#1F78B4", main="Betweenness Histogram",xlab="Betweenness")
#dev.off()

#Distances between vertices

#The inter-connectivity of a LimeGreenNetwork can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(LimeGreenNetwork, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

#pdf ("histDegreeLimeGreen.pdf")
#hist(g.outd, breaks = g.outd[which.max(g.outd)]
#, col ="#1F78B4", main="Degree #Histogram",xlab="Degree")
#dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(LimeGreenNetwork)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our LimeGreenNetwork comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoLimeGreenNetwork <- make_ego_graph(LimeGreenNetwork, 1, "UNK_5:40072703-40082576", mode=c("all"))[[1]]


png(file = "nicely.LimeGreenNetwork_ego_LimeGreen.png", width = 4000, height = 4000);
plot(egoLimeGreenNetwork,layout=layout_nicely(egoLimeGreenNetwork), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

equivalenceTable

################################################


modules_BlueJeans = "#A6CEE3"


inModule = is.finite(match(moduleColors, modules_BlueJeans));
modProbes <- probes[inModule]


#get the number of each type of small RNAs in the module

SNRNA <- length(grep("SNRNA", modProbes))

SNORNA <- length(grep("SNORNA", modProbes))

TRNA <- length(grep("TRNA", modProbes))

PHAS21 <- length(grep("21PHAS", modProbes))

PHAS24 <- length(grep("24PHAS", modProbes))

L24P_like <- length(grep("L24P", modProbes))

UNK <- length(grep("UNK", modProbes))

miRNA <- length(grep("miR", modProbes))



allTypes <- c(SNRNA, SNORNA, TRNA, PHAS21, PHAS24, L24P_like ,UNK, miRNA)
names(allTypes) <- c("sn-RNA","sno-RNA","t-RNA","21-PHAS","24-PHAS","L24P-like","UNK","mi-RNA")


pdf("All.BlueJeans.pie.pdf")
par(mar=c(2,2,2,2))
pie(allTypes,col=pastelColors,labels='')
dev.off()


pdf("All.BlueJeans.waffle.pdf")
waffle(allTypes, rows = 25, cols = pastelColors, mar = c(1,1,1,1),bg = 'white')
#legend('right', legend = LETTERS[1:4], pch = 15, col = cols, pt.cex = 2, bty = 'n')
dev.off()

adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

summary(adj[adj != 1])

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.53] = 1
adj[adj != 1] = 0


library(igraph)
(BlueJeansNetwork <- graph.adjacency(adj))
(BlueJeansNetwork <- simplify(BlueJeansNetwork))  # removes self-loops
(BlueJeansNetwork <- delete.vertices(BlueJeansNetwork, degree(BlueJeansNetwork)==0))


ColorTable <- cbind(probes,moduleColors)

V(BlueJeansNetwork)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(BlueJeansNetwork)))),2]

# remove unconnected nodes
BlueJeansNetwork <- delete.vertices(BlueJeansNetwork, degree(BlueJeansNetwork)==0)

length(names(V(BlueJeansNetwork)))

topConnectedBlueJeans <- names(V(BlueJeansNetwork))

write.table(topConnectedBlueJeans,file="topConnectedBlueJeans.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the BlueJeansNetwork. It can be thought of as how critical the vertex is to the flow of information through a BlueJeansNetwork. Individuals with high betweenness are key bridges between different parts of a BlueJeansNetwork. In our measles transmission BlueJeansNetwork, vertices with high betweenness are those children who were central to passing on the disease to other parts of the BlueJeansNetwork. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the BlueJeansNetwork adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(BlueJeansNetwork, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweeness_BlueJeans.pdf")
hist(g.b, breaks = 80)
dev.off()

png(file = "nicely.BlueJeansNetwork_BlueJeans.png", width = 4000, height = 4000);
plot(BlueJeansNetwork,layout=layout_nicely(BlueJeansNetwork), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.BlueJeansNetwork_BlueJeans.png", width = 4000, height = 4000);
plot(BlueJeansNetwork,layout=layout_with_fr(BlueJeansNetwork), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.BlueJeansNetwork_BlueJeans.png", width = 4000, height = 4000);
plot(BlueJeansNetwork,layout=layout_as_tree(BlueJeansNetwork), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.BlueJeansNetwork_BlueJeans.png", width = 4000, height = 4000);
plot(BlueJeansNetwork,layout=layout_in_circle(BlueJeansNetwork), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the BlueJeansNetwork graph.
gd <- edge_density(BlueJeansNetwork)

#Another measure of how interconnected a BlueJeansNetwork is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the BlueJeansNetwork. The longest path length between any pair of vertices is called the diameter of the BlueJeansNetwork graph.

# Get the diameter of the graph g
diameter(BlueJeansNetwork, directed = FALSE)
get_diameter(BlueJeansNetwork, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(BlueJeansNetwork)
# Shows the path sequence between two furthest apart vertices.
get_diameter(BlueJeansNetwork)

# Get the average path length of the graph g
g.apl <- mean_distance(BlueJeansNetwork, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an important method for investigating how likely or unlikely other BlueJeansNetwork metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph
Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(BlueJeansNetwork)

gorder(BlueJeansNetwork)

g.random <- erdos.renyi.game(n = gorder(BlueJeansNetwork), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random
mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

#Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(BlueJeansNetwork), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthBlueJeans.pdf")
hist(gl.apls,xlim=c(1.73,1.77),col="#1F78B4",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.003, y= 150,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,2),col="red",line=.46,cex=.9)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length is lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the BlueJeansNetwork.

g.b <- betweenness(BlueJeansNetwork, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
#pdf("histVertexBetweenessBlueJeans.pdf")
#hist(g.b,breaks = g.b[which.max(g.b)], col ="#1F78B4", main="Betweenness Histogram",xlab="Betweenness")
#dev.off()

#Distances between vertices

#The inter-connectivity of a BlueJeansNetwork can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(BlueJeansNetwork, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

#pdf ("histDegreeBlueJeans.pdf")
#hist(g.outd, breaks = g.outd[which.max(g.outd)]
#, col ="#1F78B4", main="Degree #Histogram",xlab="Degree")
#dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(BlueJeansNetwork)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our BlueJeansNetwork comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoBlueJeansNetwork <- make_ego_graph(BlueJeansNetwork, 1, "L24P.like_15:5429180-5452342", mode=c("all"))[[1]]


png(file = "nicely.BlueJeansNetwork_ego_BlueJeans.png", width = 4000, height = 4000);
plot(egoBlueJeansNetwork,layout=layout_nicely(egoBlueJeansNetwork), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

equivalenceTable

################################################


modules_YellowGreen = "#B2DF8A"


inModule = is.finite(match(moduleColors, modules_YellowGreen));
modProbes <- probes[inModule]


#get the number of each type of small RNAs in the module

SNRNA <- length(grep("SNRNA", modProbes))

SNORNA <- length(grep("SNORNA", modProbes))

TRNA <- length(grep("TRNA", modProbes))

PHAS21 <- length(grep("21PHAS", modProbes))

PHAS24 <- length(grep("24PHAS", modProbes))

L24P_like <- length(grep("L24P", modProbes))

UNK <- length(grep("UNK", modProbes))

miRNA <- length(grep("miR", modProbes))



allTypes <- c(SNRNA, SNORNA, TRNA, PHAS21, PHAS24, L24P_like ,UNK, miRNA)
names(allTypes) <- c("sn-RNA","sno-RNA","t-RNA","21-PHAS","24-PHAS","L24P-like","UNK","mi-RNA")


pdf("All.YellowGreen.pie.pdf")
par(mar=c(2,2,2,2))
pie(allTypes,col=pastelColors,labels='')
dev.off()

pdf("All.YellowGreen.waffle.pdf")
waffle(allTypes, rows = 20, cols = pastelColors, mar = c(1,1,1,1),bg = 'white')
#legend('right', legend = LETTERS[1:4], pch = 15, col = cols, pt.cex = 2, bty = 'n')
dev.off()


adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

summary(adj[adj != 1])

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.27] = 1
adj[adj != 1] = 0


library(igraph)
(YellowGreenNetwork <- graph.adjacency(adj))
(YellowGreenNetwork <- simplify(YellowGreenNetwork))  # removes self-loops
(YellowGreenNetwork <- delete.vertices(YellowGreenNetwork, degree(YellowGreenNetwork)==0))


ColorTable <- cbind(probes,moduleColors)

V(YellowGreenNetwork)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(YellowGreenNetwork)))),2]

# remove unconnected nodes
YellowGreenNetwork <- delete.vertices(YellowGreenNetwork, degree(YellowGreenNetwork)==0)

length(names(V(YellowGreenNetwork)))

topConnectedYellowGreen <- names(V(YellowGreenNetwork))

write.table(topConnectedYellowGreen,file="topConnectedYellowGreen.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the YellowGreenNetwork. It can be thought of as how critical the vertex is to the flow of information through a YellowGreenNetwork. Individuals with high betweenness are key bridges between different parts of a YellowGreenNetwork. In our measles transmission YellowGreenNetwork, vertices with high betweenness are those children who were central to passing on the disease to other parts of the YellowGreenNetwork. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the YellowGreenNetwork adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(YellowGreenNetwork, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweeness_YellowGreen.pdf")
hist(g.b, breaks = 80)
dev.off()

png(file = "nicely.YellowGreenNetwork_YellowGreen.png", width = 4000, height = 4000);
plot(YellowGreenNetwork,layout=layout_nicely(YellowGreenNetwork), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.YellowGreenNetwork_YellowGreen.png", width = 4000, height = 4000);
plot(YellowGreenNetwork,layout=layout_with_fr(YellowGreenNetwork), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.YellowGreenNetwork_YellowGreen.png", width = 4000, height = 4000);
plot(YellowGreenNetwork,layout=layout_as_tree(YellowGreenNetwork), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.YellowGreenNetwork_YellowGreen.png", width = 4000, height = 4000);
plot(YellowGreenNetwork,layout=layout_in_circle(YellowGreenNetwork), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the YellowGreenNetwork graph.
gd <- edge_density(YellowGreenNetwork)

#Another measure of how interconnected a YellowGreenNetwork is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the YellowGreenNetwork. The longest path length between any pair of vertices is called the diameter of the YellowGreenNetwork graph.

# Get the diameter of the graph g
diameter(YellowGreenNetwork, directed = FALSE)
get_diameter(YellowGreenNetwork, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(YellowGreenNetwork)
# Shows the path sequence between two furthest apart vertices.
get_diameter(YellowGreenNetwork)

# Get the average path length of the graph g
g.apl <- mean_distance(YellowGreenNetwork, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an important method for investigating how likely or unlikely other YellowGreenNetwork metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph
Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(YellowGreenNetwork)

gorder(YellowGreenNetwork)

g.random <- erdos.renyi.game(n = gorder(YellowGreenNetwork), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random
mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

#Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(YellowGreenNetwork), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthYellowGreen.pdf")
hist(gl.apls,xlim=c(1.73,1.77),col="#1F78B4",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.003, y= 150,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,2),col="red",line=.46,cex=.9)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length is lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the YellowGreenNetwork.

g.b <- betweenness(YellowGreenNetwork, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
#pdf("histVertexBetweenessYellowGreen.pdf")
#hist(g.b,breaks = g.b[which.max(g.b)], col ="#1F78B4", main="Betweenness Histogram",xlab="Betweenness")
#dev.off()

#Distances between vertices

#The inter-connectivity of a YellowGreenNetwork can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(YellowGreenNetwork, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

#pdf ("histDegreeYellowGreen.pdf")
#hist(g.outd, breaks = g.outd[which.max(g.outd)]
#, col ="#1F78B4", main="Degree #Histogram",xlab="Degree")
#dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(YellowGreenNetwork)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our YellowGreenNetwork comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoYellowGreenNetwork <- make_ego_graph(YellowGreenNetwork, 1, "21PHAS_8:9616284-9619710", mode=c("all"))[[1]]


png(file = "nicely.YellowGreenNetwork_ego_YellowGreen.png", width = 4000, height = 4000);
plot(egoYellowGreenNetwork,layout=layout_nicely(egoYellowGreenNetwork), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

equivalenceTable

################################################


modules_LightSalmon = "#FB9A99"


inModule = is.finite(match(moduleColors, modules_LightSalmon));
modProbes <- probes[inModule]


#get the number of each type of small RNAs in the module

SNRNA <- length(grep("SNRNA", modProbes))

SNORNA <- length(grep("SNORNA", modProbes))

TRNA <- length(grep("TRNA", modProbes))

PHAS21 <- length(grep("21PHAS", modProbes))

PHAS24 <- length(grep("24PHAS", modProbes))

L24P_like <- length(grep("L24P", modProbes))

UNK <- length(grep("UNK", modProbes))

miRNA <- length(grep("miR", modProbes))



allTypes <- c(SNRNA, SNORNA, TRNA, PHAS21, PHAS24, L24P_like ,UNK, miRNA)
names(allTypes) <- c("sn-RNA","sno-RNA","t-RNA","21-PHAS","24-PHAS","L24P-like","UNK","mi-RNA")


pdf("All.LightSalmon.pie.pdf")
par(mar=c(2,2,2,2))
pie(allTypes,col=pastelColors,labels='')
dev.off()


pdf("All.LightSalmon.waffle.pdf")
waffle(allTypes, rows = 20, cols = pastelColors, mar = c(1,1,1,1),bg = 'white')
#legend('right', legend = LETTERS[1:4], pch = 15, col = cols, pt.cex = 2, bty = 'n')
dev.off()


adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

summary(adj[adj != 1])

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.18] = 1
adj[adj != 1] = 0


library(igraph)
(LightSalmonNetwork <- graph.adjacency(adj))
(LightSalmonNetwork <- simplify(LightSalmonNetwork))  # removes self-loops
(LightSalmonNetwork <- delete.vertices(LightSalmonNetwork, degree(LightSalmonNetwork)==0))


ColorTable <- cbind(probes,moduleColors)

V(LightSalmonNetwork)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(LightSalmonNetwork)))),2]

# remove unconnected nodes
LightSalmonNetwork <- delete.vertices(LightSalmonNetwork, degree(LightSalmonNetwork)==0)

length(names(V(LightSalmonNetwork)))

topConnectedLightSalmon <- names(V(LightSalmonNetwork))

write.table(topConnectedLightSalmon,file="topConnectedLightSalmon.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the LightSalmonNetwork. It can be thought of as how critical the vertex is to the flow of information through a LightSalmonNetwork. Individuals with high betweenness are key bridges between different parts of a LightSalmonNetwork. In our measles transmission LightSalmonNetwork, vertices with high betweenness are those children who were central to passing on the disease to other parts of the LightSalmonNetwork. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the LightSalmonNetwork adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(LightSalmonNetwork, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweeness_LightSalmon.pdf")
hist(g.b, breaks = 80)
dev.off()

png(file = "nicely.LightSalmonNetwork_LightSalmon.png", width = 4000, height = 4000);
plot(LightSalmonNetwork,layout=layout_nicely(LightSalmonNetwork), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.LightSalmonNetwork_LightSalmon.png", width = 4000, height = 4000);
plot(LightSalmonNetwork,layout=layout_with_fr(LightSalmonNetwork), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.LightSalmonNetwork_LightSalmon.png", width = 4000, height = 4000);
plot(LightSalmonNetwork,layout=layout_as_tree(LightSalmonNetwork), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.LightSalmonNetwork_LightSalmon.png", width = 4000, height = 4000);
plot(LightSalmonNetwork,layout=layout_in_circle(LightSalmonNetwork), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the LightSalmonNetwork graph.
gd <- edge_density(LightSalmonNetwork)

#Another measure of how interconnected a LightSalmonNetwork is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the LightSalmonNetwork. The longest path length between any pair of vertices is called the diameter of the LightSalmonNetwork graph.

# Get the diameter of the graph g
diameter(LightSalmonNetwork, directed = FALSE)
get_diameter(LightSalmonNetwork, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(LightSalmonNetwork)
# Shows the path sequence between two furthest apart vertices.
get_diameter(LightSalmonNetwork)

# Get the average path length of the graph g
g.apl <- mean_distance(LightSalmonNetwork, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an important method for investigating how likely or unlikely other LightSalmonNetwork metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph
Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(LightSalmonNetwork)

gorder(LightSalmonNetwork)

g.random <- erdos.renyi.game(n = gorder(LightSalmonNetwork), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random
mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

#Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(LightSalmonNetwork), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthLightSalmon.pdf")
hist(gl.apls,xlim=c(1.73,1.77),col="#1F78B4",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.003, y= 150,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,2),col="red",line=.46,cex=.9)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length is lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the LightSalmonNetwork.

g.b <- betweenness(LightSalmonNetwork, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
#pdf("histVertexBetweenessLightSalmon.pdf")
#hist(g.b,breaks = g.b[which.max(g.b)], col ="#1F78B4", main="Betweenness Histogram",xlab="Betweenness")
#dev.off()

#Distances between vertices

#The inter-connectivity of a LightSalmonNetwork can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(LightSalmonNetwork, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

#pdf ("histDegreeLightSalmon.pdf")
#hist(g.outd, breaks = g.outd[which.max(g.outd)]
#, col ="#1F78B4", main="Degree #Histogram",xlab="Degree")
#dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(LightSalmonNetwork)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our LightSalmonNetwork comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoLightSalmonNetwork <- make_ego_graph(LightSalmonNetwork, 1, "SNORNA_13:12799762-12799886", mode=c("all"))[[1]]


png(file = "nicely.LightSalmonNetwork_ego_LightSalmon.png", width = 4000, height = 4000);
plot(egoLightSalmonNetwork,layout=layout_nicely(egoLightSalmonNetwork), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

equivalenceTable

################################################
