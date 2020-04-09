
# Relating relative prey contribution and dietary diversity metrics

This contents are related with the article "Harmonising relative prey contribution and dietary diversity metrics in DNA metabarcoding studies" that will be soon sent to a peer-review journal.

### Four basic predator-prey matrices

All relative prey contribution and dietary diversity metrics are computed based on the following four predator-prey matrices. Note that predator can refer to a replicate, a sample, an individual or any other unit defined by researchers, and prey can refer to unique DNA haplotypes, ASVs, OTUs, species or another type defined by researchers.

````R
#Create the original "Absolute Prey Abundance Matrix" (APAM) - Figure 1A
APAM <- cbind(Predator1=c(10,0,20,30),Predator2=c(10,0,10,0),Predator3=c(20,25,20,5))
rownames(APAM) <- c("Prey1","Prey2","Prey3","Prey4")
APAM

#Transform APAM to "Absolute Prey Occurrence Matrix" (APOM) - Figure 1B
APOM <- APAM
APOM[APOM != 0] <- 1

#Transform APAM to "Normalised Prey Abundance Matrix" (NPAM) - Figure 1C
NPAM <- apply(APAM,2,function(x){x/sum(x)})
NPAM

#Transform APOM to "Normalised Prey Occurrence Matrix" (NPOM) - Figure 1D
NPOM <- apply(APOM,2,function(x){x/sum(x)})
NPOM
````

--ADD Figure 1--

### Relative prey contribution metrics

The different relative prey contribution metrics are computed by dividing the predator sum values by the total sum. 
````R
#Relative Prey Abundance (RPA) - Figure 2A
APAMsums <- rowSums(APAM)
RPA <- APAMsums / sum(APAM)
#or
RPA <- APAMsums / sum(APAMsums)
RPA

#Relative Prey Occurrence (RPO) - Figure 2B
APOMsums <- rowSums(APOM)
RPO <- APOMsums / sum(APOM)
#or
RPO <- APOMsums / sum(APOMsums)
RPO

#Normalised Relative Prey Abundance (nRPA) - Figure 2C
NPAMsums <- rowSums(NPAM)
nRPA <- NPAMsums / sum(NPAM)
#or
nRPA <- NPAMsums / sum(NPAMsums)
nRPA

#Normalised Relative Prey Occurrence (nRPO) - Figure 2D
NPOMsums <- rowSums(NPOM)
nRPO <- NPOMsums / sum(NPOM)
#or
nRPO <- NPOMsums / sum(NPOMsums)
nRPO
````

--ADD Figure 2--

### Diversity metrics

--ADD Figure 3--

#### Abundance-based diversity of a single predator - Figure 3A

If q=0, the computed metric is prey richness. If q=1, the computed metric is Shannon diversity (exponential of Shannon index). If q=2, the computed metric is Simpson diversity (the multiplicative inverse of Simpson index).

````R
predator="Predator1" #Can be changed to "Predator2" or "Predator3"
q=1 #Can be changed to 0, 2 or any other positive value
if(q == 1){q=0.999999999} # Because the function is not defined for the unity
pi <- NPAM[,predator]
pi.q <- pi^q
basicsum <- sum(pi.q)
basicsum^(1/(1-q))
````

Luckily, functions to easily compute individual diversity metrics exist.
````R
library(hilldiv)
hill_div(NPAM,qvalue=1)
````

####  Abundance-based diversity of the entire predator system (gamma diversity) - Figure 3B

Using even predator weights.
````R
evenweight=c(1/3,1/3,1/3)
q=2 #Can be changed to 0, 2 or any other positive value
if(q == 1){q=0.999999999} # Because the function is not defined for the unity
pi <- NPAM
pi.w <- sweep(pi, 2, evenweight, "*")
pi.w.sum <- rowSums(pi.w)
pi.w.sum.q <- pi.w.sum^q
basicsum <- sum(pi.w.sum.q)
basicsum^(1/(1-q))
````

Using predator weights leveled to sequencing depth
````R
seqdepthweight <- colSums(APAM) / sum(colSums(APAM))
q=1 #Can be changed to 0, 2 or any other positive value
if(q == 1){q=0.999999999} # Because the function is not defined for the unity
pi <- NPAM
pi.w <- sweep(pi, 2, seqdepthweight, "*")
pi.w.sum <- rowSums(pi.w) # Note that this is identical to RPA
pi.w.sum.q <- pi.w.sum^q
basicsum <- sum(pi.w.sum.q)
basicsum^(1/(1-q))
````

Luckily, functions to easily compute individual diversity metrics exist
````R
library(hilldiv)
gamma_div(NPAM,qvalue=1)
gamma_div(NPAM,qvalue=1,weight=evenweight)
gamma_div(NPAM,qvalue=1,weight=seqdepthweight)
````

#### Incidence-based diversity of the entire predator system (gamma diversity) - Figure 3C
````R
q=2 #Can be changed to 0, 2 or any other positive value
if(q == 1){q=0.999999999} # Because the function is not defined for the unity
pi <- RPO
pi.q <- pi^q
basicsum <- sum(pi.q)
basicsum^(1/(1-q))
````
Luckily, functions to easily compute this exist
````R
library(hilldiv)
hill_div(to.incidence(APAM),qvalue=1)
hill_div(RPO,qvalue=1)
````
