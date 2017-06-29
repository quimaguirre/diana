library(RCurl)
library(rjson)

association.file.location='/home/quim/project/diana/toolbox/funcassociate_my_associations.tsv'  ### the content of this file defines each Attribute (GO) and their Entities (Genes) 
  #### note 1: there is a limit in the numberof character per row in the association file <=32000 
  #### note 2: the lined endigs must be LF (Unix style); including the last line
  #### note 3: see the acompanying association file for an example of format (3 columns, no comments ...)

genespace.file.location  ='/home/quim/project/diana/toolbox/funcassociate_my_genespace_weighted.tsv' ### the content of this file defines Entity Space (Genes) and their weights (weights must be non zero integers) 

funcassociate.result.destination='/home/quim/project/diana/toolbox/funcassociate_result.xls'  ### optional we store the result as a tab delimited file

####### Example : using an association file defined by user (bypassing GO) and a user defined GeneSpace with weights
query="HEPACAM AAK1 ABL1 ABL2 ACTR2 ADCK1 ADCK2 ADCK3 ADCK4 ADCK5 AKT1 AKT2 AKT3 ALK ARAF ATM ATR AXL BCKDK BCR BLK BMPR1A BMPR1B BMPR2 BMX BRAF BRD2 BRD3 BRD4 BRDT BRSK1 BRSK2 BTK BUB1 CAMK1G CAMK2A CAMK2B CAMK2D CAMK2G CAMK4 CAMKK1 CASK CDC7 CDK10 CDK2 CDK3 CDK4 CDK5 CDK6 CDK7 CDK8 CDK9 CDKL1 CDKL2" ### large number of Kinases
associations.file=fileUpload(filename='funcassociate_my_associations.tsv'  ### an arbitrary name that is transmitted to the Funcassociate Web Server
                            ,contents=paste(paste(readLines(association.file.location),collapse="\n"),"\n",sep=''));  ### read line by line and separate by Unix LineFeed 
genespace.file   =fileUpload(filename='funcassociate_my_genespace_weighted.tsv'  ### an arbitrary name that is transmitted to the Funcassociate Web Server
                            ,contents=paste(paste(readLines(genespace.file.location),collapse="\n"),"\n",sep=''));  ### read line by line and separate by Unix LineFeed 



#### define the other params needed by FuncAssociate Client
reps = 1000;
mode = 'unordered'
which = 'over'
cutoff= 0.05

##### obtain the results as a JSON text after we post the form to the FuncAssociate Server
my.text<-postForm('http://llama.mshri.on.ca/cgi/funcassociate/serv'
                  ,method='functionate'
                  ,query =query
                  ,`associations-file`=associations.file
                  ,`genespace-file`   =genespace.file
                  ,reps  =reps
                  ,mode  =mode
                  ,which =which
                  ,cutoff=cutoff)

##### parse the JSON response into an R structure
my.str<-fromJSON(my.text) 

##### look at the warnings
print(my.str$result$entity_attribute_table$warnings)

##### assembly the table with results
result.table<-as.matrix(Reduce(rbind,lapply(my.str$result$over,unlist)))
colnames(result.table)<-c('nmb entities in query for row','nmb entities in query','nmb entities in TOTAL for row','log10 ods','pValue','pValue adjusted by resampling','GO term','GO term description')
for (j in 1:6) result.table[,j]<-as.numeric(result.table[,j])

##### write the result as an excel tab delimited file (this is optional)
write.table(result.table, funcassociate.result.destination, sep="\t",quote=FALSE,row.names=FALSE)
