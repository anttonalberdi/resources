
# Relating relative prey contribution and dietary diversity metrics

This contents are related with the article "Harmonising relative prey contribution and dietary diversity metrics in DNA metabarcoding studies" that will be soon sent to a peer-review journal.

````R
#Create the original predator-prey matrix
matrix <- cbind(Predator1=c(10,0,20,30),Predator2=c(10,0,10,0),Predator3=c(20,25,20,5))
rownames(matrix) <- c("Prey1","Prey2","Prey3","Prey4")

#Visualise it
matrix
````
