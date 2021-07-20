#' @title Converts MAF into mutation matrix.
#' @description The function `get_mut_matrix` converts mutation annotation file (MAF) format data into a mutations matrix.Then use the fisher exact test to select the geneset with higher mutation frequency in alive sample group.Finally return the higher mutation frequency matrix.
#' @param maffile Input mutation annotation file (MAF) format data. It must be an absolute path or the name relatived to the current working directory.
#' @param is.TCGA 	Is input MAF file from TCGA source. If TRUE uses only first 15 characters from Tumor_Sample_Barcode.
#' @param mut_fre A threshold value(zero as the default value). The genes with a given mutation frequency equal or greater than the threshold value are retained for the following analysis.
#' @param nonsynonymous Logical,tell if extract the non-synonymous somatic mutations (nonsense mutation, missense mutation, frame-shif indels, splice site, nonstop mutation, translation start site, inframe indels).
#' @param cut_fisher.pval The significant cut_off pvalue for the fisher test.
#' @param cut_oddsRatio The cut_off oddsRatio for the fisher test,uses to select the geneset with higher mutation frequency in alive sample group.
#' @param sur A nx2 data frame of samples' survival data,the first line is samples' survival event and the second line is samples' overall survival.
#' @return The survival-related mutations matrix.
#' @importFrom utils read.delim
#' @importFrom maftools read.maf
#' @importFrom stats fisher.test
#' @export
#' @examples
#' #get the path of the mutation annotation file and samples' survival data
#' maf<-system.file("extdata","data_mutations_extended.txt",package = "pathwayTMB")
#' sur_path<-system.file("extdata","sur.csv",package = "pathwayTMB")
#' sur<-read.csv(sur_path,header=TRUE,row.names = 1)
#' #perform the function 'get_mut_matrix'
#' mut_matrix<-get_mut_matrix(maffile=maf,mut_fre = 0.01,is.TCGA=FALSE,sur=sur)
get_mut_matrix<-function(maffile,is.TCGA=TRUE,mut_fre=0,nonsynonymous = TRUE,cut_fisher.pval=1,cut_oddsRatio=1,sur){
  if(nonsynonymous){
    data<-as.data.frame(read.maf(maffile)@data,stringsAsFactors = default.stringsAsFactors())
  }else{data<-read.delim( maffile, header = TRUE, comment.char = '#', stringsAsFactors = FALSE )}
  if(is.TCGA){
    data_9<-data[,c(1,2,9,10,16)]
    for(i in 1:dim(data_9)[2]){
      data_9[,i]<-as.character(data_9[,i])
    }
    data_9[,5]<-substr(gsub(pattern = "\\-",replacement = "\\.",x =data_9[,5]))
  }else{
    data_9<-data[,c(1,2,10,11,17)]
  }
  samples<-unique(data_9[,5])
  genes<-unique(data_9[,1])
  freq_matrix<-matrix(0,length(genes),length(samples))
  rownames(freq_matrix)<-genes
  colnames(freq_matrix)<-samples
  for(i in 1:length(samples)){
    s1<-table(data_9[which(data_9[,5]==samples[i]),1])
    freq_matrix[names(s1),i]<-s1
  }
  a<-apply(freq_matrix,1,function(x){(length(which(x!=0))/length(x))})
  freq_matrix<-freq_matrix[which(a >= mut_fre),]
  #fisher.test
  mut_matrix<-freq_matrix
  inter<-intersect(colnames(mut_matrix),rownames(sur))
  mut_matrix<-mut_matrix[,inter]
  sur<-sur[inter,]
  mut_matrix[mut_matrix>0]<-1
  pvalue<-c()
  oddsRatio<-c()
  for(i in 1:length(mut_matrix[,1])){
    a<-length(intersect(which(mut_matrix[i,]==1),which(sur[,1]==0)))
    b<-length(intersect(which(mut_matrix[i,]==1),which(sur[,1]==1)))
    c<-length(intersect(which(mut_matrix[i,]==0),which(sur[,1]==0)))
    d<-length(intersect(which(mut_matrix[i,]==0),which(sur[,1]==1)))
    x<-matrix(c(a,c,b,d),ncol=2,nrow=2)
    p<-fisher.test(x)$p.value
    pvalue<-c(pvalue,p)
    OR<-fisher.test(x)$estimate
    oddsRatio<-c(oddsRatio,OR)
  }
  names(pvalue)<-rownames(mut_matrix)
  names(oddsRatio)<-rownames(mut_matrix)
  gene<-intersect(names(which(oddsRatio>cut_oddsRatio)),names(which(pvalue<cut_fisher.pval)))
  freq_matrix<-freq_matrix[gene,]
  return(freq_matrix)
}

#' @title Extract coding genes' length from gene transfer format(GTF) file.
#' @description Extract coding genes' length from gene transfer format(GTF) file.
#' @param filepath Input gene transfer format file. It must be an absolute path or the name relatived to the current working directory.
#' @return Return a list of genes' length.
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicFeatures exonsBy
#' @importFrom purrr reduce
#' @importFrom BiocGenerics width
#' @importFrom clusterProfiler bitr
#' @export
get_gene_length<-function(filepath){
  txdb <-makeTxDbFromGFF(filepath,format="gtf")
  exon<-exonsBy(txdb, "gene")
  exonic.gene.sizes <- lapply(exon,function(x){sum(width(reduce(x)))})
  a<-names(exonic.gene.sizes)
  c<-c()
  for(i in 1:length(a)){
    b<-strsplit(a[i],split = "\\.")[[1]][1]
    c<-c(c,b)
  }
  exonic.gene.sizes<-exonic.gene.sizes[-which(duplicated(c))]
  uni_c<-unique(c)
  eg<-bitr(uni_c,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
  b<-c()
  for(i in 1:length(eg[,1])){
    a<-which(grepl(eg[i,1],names(exonic.gene.sizes),ignore.case = T))
    b<-c(b,a)
  }
  genesmbol<-exonic.gene.sizes[b]
  names(genesmbol)<-eg[,2]
  genesmbol<-genesmbol[-which(duplicated(names(genesmbol)))]
  return(genesmbol)}

#' @title Calculate the Pathway-based Tumor Mutational Burden.
#' @description The function `get_PTMB` uses to calculate the Pathway-based Tumor Mutational Burden (PTMB). PTMB is defined as pathway-based tumor mutational burden corrected by genes’ length and number.
#' @param freq_matrix The mutations matrix,generated by `get_mut_matrix`.
#' @param genesmbol The genes' length matrix,generated by `get_gene_length`.
#' @param path_mut_cutoff A threshold value(zero percent as the default value).Pathways with a given mutation frequency equal or greater than the threshold value are retained for the following analysis.
#' @param gene_path User input pathways geneset list.
#' @return Return the Pathway-based Tumor Mutational Burden matrix.
#' @export
#' @examples
#' #get the path of the mutation annotation file and samples' survival data
#' \donttest{maf<-system.file("extdata","data_mutations_extended.txt",package = "pathwayTMB")
#' sur_path<-system.file("extdata","sur.csv",package = "pathwayTMB")
#' sur<-read.csv(sur_path,header=TRUE,row.names = 1)
#' #perform the function 'get_mut_matrix'
#' mut_matrix<-get_mut_matrix(maffile=maf,mut_fre = 0.01,is.TCGA=FALSE,sur=sur)}
#' #perform the function `get_PTMB`
#' PTMB_matrix<-get_PTMB(freq_matrix=mut_matrix,genesmbol=genesmbol,gene_path=gene_path)
get_PTMB<-function(freq_matrix,genesmbol,path_mut_cutoff=0,gene_path){
  inter_genesymbol<-intersect(rownames(freq_matrix),names(genesmbol))
  inter_freq<-freq_matrix[inter_genesymbol,]
  inter_length<-genesmbol[inter_genesymbol]
  leng<-c()
  for(i in 1:length(inter_length)){
    a<-inter_length[[i]]
    leng<-c(leng,a)
  }
  inter_length<-as.data.frame(cbind(symbol=as.character(names(inter_length)),length=as.numeric(as.character((leng)))))
  rownames(inter_length)<-inter_length[,1]
  inter_length[,2]<-as.numeric(as.character(inter_length[,2]))
  path_TMB<-data.frame()
  for(i in 1:length(gene_path)){
    gene<-intersect(gene_path[[i]],rownames(inter_freq))
    if(length(gene)>1){
      b<-c()
      for(j in 1:length(gene)){
        a<-inter_freq[gene[j],]/(as.numeric(inter_length[gene[j],2])/1000000)
        b<-rbind(b,a)
      }
      c<-colSums(b)/length(gene_path[[i]])
    }
    else if(length(gene)==1){
      b<-inter_freq[gene,]/(as.numeric(inter_length[gene,2])/1000000)
      c<-b/length(gene_path[[i]])}
    else{c=rep(0,length(inter_freq[1,]));names(c)<-colnames(inter_freq)}
    path_TMB<-rbind(path_TMB,c)
  }
  colnames(path_TMB)<-colnames(inter_freq)
  rownames(path_TMB)<-names(gene_path)
  a<-apply(path_TMB,1,function(x){length(which(x!=0))})
  if(length(which(a>=(length(path_TMB[1,])*path_mut_cutoff)))>0){
    path_TMB_cutoff<-path_TMB[which(a>=(length(path_TMB[1,])*path_mut_cutoff)),]
  }
  else{path_TMB_cutoff<-path_TMB}
  return(path_TMB_cutoff)
}


#' @title Filter cancer-specific dysfunction pathways.
#' @description The function `get_final_signature` , using to filter cancer-specific dysfunction pathways (a potential marker for cancer prognostic and immunotherapy), is the main function of our analysis.
#' @param PTMB The pathway tumor mutation burden matrix,generated by`get_PTMB`.
#' @param sur A nx2 data frame of samples' survival data,the first line is samples' survival event and the second line is samples' overall survival.
#' @importFrom stats wilcox.test
#' @importFrom stats median
#' @importFrom caret rfe
#' @importFrom randomForest randomForest
#' @importFrom caret predictors
#' @importFrom caret rfeControl
#' @importFrom caret rfFuncs
#' @importFrom survival Surv
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#' @importFrom survival survfit
#' @importFrom survminer surv_pvalue
#' @return Return the final PTMB signature,could be a potential marker for prognostic and immunotherapy prediction.
#' @export
#' @examples
#' #get the path of the mutation annotation file and samples' survival data
#' maf<-system.file("extdata","data_mutations_extended.txt",package = "pathwayTMB")
#' \donttest{sur_path<-system.file("extdata","sur.csv",package = "pathwayTMB")
#' sur<-read.csv(sur_path,header=TRUE,row.names = 1)
#' #perform the function 'get_mut_matrix'
#' mut_matrix<-get_mut_matrix(maffile=maf,mut_fre = 0.01,is.TCGA=FALSE,sur=sur)
#' #perform the function `get_PTMB`
#' PTMB_matrix<-get_PTMB(freq_matrix=mut_matrix,genesmbol=genesmbol,gene_path=gene_path)
#' set.seed(1)
#' final_character<-get_final_signature(PTMB=PTMB_matrix,sur=sur)}
get_final_signature<-function(PTMB,sur){
  inter<-intersect(colnames(PTMB),rownames(sur))
  if(length(inter)==0){
    stop("please input the same sample data.")
  }
  path_TMB_inter<-PTMB[,inter]
  sur<-sur[inter,]
  colnames(sur)<-c("event","survival")
  pvalue<-apply(path_TMB_inter,1,function(x){wilcox.test(as.numeric(x[which(sur[,1]==0)]),as.numeric(x[-which(sur[,1]==0)]))$p.value})
  FC<-apply(path_TMB_inter[,which(sur[,1]==0)],1,mean)/apply(path_TMB_inter[,-which(sur[,1]==0)],1,mean)
  DE_path<-names(pvalue)[intersect(which(pvalue<0.01),which(FC>2|FC<1/2))]
  #DE_path<-names(pvalue)[which(pvalue<0.05)]
  DE_path_sur<-as.data.frame(t(rbind(path_TMB_inter[DE_path,],event=sur[,1])))
  rfProfile <- rfe(DE_path_sur[,1:(dim(DE_path_sur)[2]-1)],as.factor(DE_path_sur[,(dim(DE_path_sur)[2])]), sizes = c(1:(dim(DE_path_sur)[2]-1)),rfeControl = rfeControl(functions = rfFuncs,method = "cv"))
  final_character<-predictors(rfProfile)
  #lasso回归
  data3<-as.data.frame(DE_path_sur[,c(final_character,"event")])
  x = as.matrix(data3[,-dim(data3)[2]])
  sur[,1]<-as.double(sur[,1])
  sur[,2]<-as.double(sur[,2])
  surs<-Surv(time = sur[,2],event = sur[,1])
  lasso_fit<-cv.glmnet(x,surs,family="cox",type.measure = "deviance")
  coeff<-coef(lasso_fit,s=lasso_fit$lambda.min)
  final_character<-rownames(coeff)[which(as.numeric(coeff)!=0)]
  #
  final_sur<-as.data.frame(cbind(DE_path_sur[,final_character],event=sur[,1],survival=sur[,2]))
  rownames(final_sur)<-rownames(DE_path_sur)
  DEpathname<-colnames(final_sur)[1:length(final_character)]
  colnames(final_sur)<-c(paste0("a",1:length(final_character)),"event","survival")
  name2name<-cbind(pathname=DEpathname,na=paste0("a",1:length(final_character)))
  rownames(name2name)<-name2name[,2]
  data1<-final_sur
  pvalue_logRank<-c()
  for(i in 1:length(final_character)){
    data<-data1[,c(i,dim(data1)[2]-1,dim(data1)[2])]
    a<-which(data[,1]<median(data[,1]))
    data[which(data[,1]<median(data[,1])),1]<-"low_risk"
    data[-a,1]<-"high_risk"
    colnames(data)[1]<-"class"
    fit<-survfit(Surv(survival, event) ~class,data = as.data.frame(data))
    b<-surv_pvalue(fit,data = data)$pval
    pvalue_logRank<-c(pvalue_logRank,b)
  }
  names(pvalue_logRank)<-colnames(data1)[1:length(final_character)]
  final_character<-names(pvalue_logRank)[which(pvalue_logRank<0.05)]
  final_character<-name2name[final_character,1]
  return(final_character)
}

