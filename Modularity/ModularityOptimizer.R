#Run ModularityOptimizer in R
#Requires Java JDK 7

#w: Adjacency matrix from doSNN
#edge: Name of the input edge file
#output: Name of the output file
#modularity_function:	Modularity function (1 = standard; 2 = alternative)
#resolution_parameter: Value of the resolution parameter; 
#Use a value of 1.0 for standard modularity-based community detection. Use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#optimization_algorithm: Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm)
#n_start: Number of random starts
#n_iter: Maximal number of iterations per random start
#random_seed: Seed of the random number generator
#print_output: Whether or not to print output to the console (0 = no; 1 = yes)

doModularity_Clust=function(object, w=matrix(), output, modularity_function=1, resolution_param=1.0, algorithm=1, n_start=1000, n_iter=10, random_seed=0, print_output=1, set.ident=TRUE){
  diag(w)=0
  edge=cbind((which(w!=0,arr.ind = TRUE)-1),w[which(w!=0,arr.ind = TRUE)])
  rownames(edge)=NULL; colnames(edge)=NULL
  write.table(x = edge,file = "edge.txt",sep = "\t",row.names = FALSE,col.names = FALSE)
  
  command=paste("java -jar ModularityOptimizer.jar", "edge.txt", output, modularity_function, resolution_param, algorithm, n_start, n_iter, random_seed, print_output, sep = " ")
  system(command, wait=TRUE)
  
  if (set.ident){
    ident=read.table(file = output,header = FALSE,sep = "\t")[,1]
    object=set.ident(object,ident.use=ident)
    return(object)
  }
}

