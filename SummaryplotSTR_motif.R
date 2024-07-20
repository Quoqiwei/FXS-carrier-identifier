library(data.table)
library(ggplot2)
library(argparser)
options(warn=-1)
wkdir=getwd()
p <- arg_parser("Statistically analyze and plot STR results")
p <- add_argument(p, "--name", short = '-n',help='sample name',default = 'Sample')
p <- add_argument(p, "--input", short = '-i',help='str summary file',default = 'STR_summary_hcov.tsv')
p <- add_argument(p, "--input2", short = '-i2',help='str summary file',default = 'STR_summary_hcov.tsv')
p <- add_argument(p, "--output", short = '-o',help='str summary file',default = wkdir)
p <- add_argument(p, "--strnum", short = '-S',help='str summary file',default = 10 )
p <- add_argument(p, "--lowper", short = '-l',help='threshold of AGG interruption',default = 0 )
p <- add_argument(p, "--thrvalue", short = '-tv',help='threads value',default = 53 )
p <- add_argument(p, "--motif", short = '-m',help='motif',default = 'AGG' )
p <- add_argument(p, "--motif_raw", short = '-mr',help='motif_raw',default = 'CGG' )
p <- add_argument(p, "--Mut", short = '-t',help='mut',default = 'Mut' )
argv <- parse_args(p)


num = args
afile <- fread(argv$input,col.names = c('r1','s1','e1','r2','s2','e2','strand','strlen','strnum','qvalue','tarlen','cov','sam','STR_num'))
Alleledb <- fread(argv$input2,col.names = c('ID','STR_len','STR_num','qvalue','tarlen','cov','sam'))
afile$STR_num <- round(afile$STR_num)
Alleledb$STR_num <- round(Alleledb$STR_num)
afile <- na.omit(afile)
afile <- afile[afile$STR_num > 0,]

getmode <- function(ndb,tmpx) {
      tse = ceiling(tmpx/10)
      topldb <-  ndb[(ndb$STR_num>tmpx-tse)&(ndb$STR_num<tmpx+tse),]
      v = topldb$STR_num
      uniqv <- unique(v)
      maxx <- uniqv[which.max(tabulate(match(v, uniqv)))]
      nrow(topldb)
    }
testfun <- function(x,nnu,rawfile,maxy){
  list1 = 0
  list2 = 0
  topl <- which.max(x$y)
  toplx <- x$x[topl]
  maxnum <- maxy
  for(i in c(2:(length(x$y)-1))){
    tmpx = ceiling(x$x[i])
    tmpy = getmode(rawfile,tmpx)
    if(x$y[i]>x$y[i-1]){
      if(x$y[i]>x$y[i+1]){
	if(argv$lowper==0){
	  lowper = 0
	}else{
	  if(argv$lowper > 1){
		  lowper = argv$lowper
	  }else{
	  	lowper = maxnum*argv$lowper
	  }
	}
        if(tmpy >= lowper){
          list1 <- c(list1,x$y[i])
          list2 <- c(list2,round(x$x[i],digits = 2))
        }
        
      }
    }
  }
  newfile <- data.frame(cbind(list1,list2))
  newfile <- newfile[newfile$list1!=0,]
  newfilet <- newfile[newfile$list1 %in% sort(newfile$list1,decreasing = T)[1:nnu],]
  newfilet <- newfilet[order(newfilet$list1,decreasing = T),]
  return(newfilet)
}
SumNdb <- read.table('STR_statis.tsv',header=T)
SumNdb <- SumNdb[SumNdb$Sample_name==argv$name,]
reads_number <- list()

for(j in SumNdb[,sprintf('%s_repeat',argv$motif_raw)]){
  tmpalist <- Alleledb[Alleledb$STR_num==as.numeric(j),]$ID
  tmpafile <- afile[afile$r1 %in% tmpalist, ]
  Allnum = length(tmpalist)
  print(sprintf('%s_repeats: %s, reads number: %s',argv$motif_raw,j,Allnum))
  if(nrow(tmpafile)>0){
    defile <- density.default(tmpafile$STR_num,bw=1,kernel = 'gaussian',n = 1024)
    defile$x[which.max(defile$y)]
    Strnum <- testfun(defile,argv$strnum,tmpafile,Allnum)
    label1 <- round(Strnum$list2)
    Strnum1 <- testfun(defile,(argv$strnum)*2,tmpafile,Allnum)
    label2 <- round(Strnum1$list2)
    getbig <- function(tdata,num) {
      return(tdata)
    }
    Rsitef = 0
    newlabel = c()
    drawlim = c()

   
      for(i in Strnum[order(Strnum$list2),]$list2){
  Rsite_old = 0
  Rsite = round(i)

  while (Rsite != Rsite_old) {
    Se = ceiling(Rsite/10)
    tmpfile = tmpafile[tmpafile$STR_num>Rsite-Se,]
    tmpfile = tmpfile[tmpfile$STR_num<Rsite+Se,]
    tdb = data.frame(table(tmpfile$STR_num))
    if(nrow(tdb)>=1){
	    Rsite_old = Rsite
	    Rsite = as.numeric(as.character(tdb[which.max(tdb$Freq),'Var1']))
    }else{
	    Rsite_old = Rsite
    }
  }
  
    tmpnum = nrow(tmpfile)
    tmpnum2 = nrow(tmpafile[tmpafile$STR_num == Rsite,])
    tmpnum3 = nrow(tmpafile[tmpafile$STR_num == Rsitef,])
    reads_number<- c(tmpnum2,tmpnum3)
    max_y = max(reads_number)
    max_y = round(max_y/0.8)
  
    if(Rsite != Rsitef){
	  if(Rsite <= j){
    print(sprintf('The position of AGG interruption within the CGG repeat is %s, with %s reads. The number of reads within ±%s of %s AGG interruption are %s.',Rsite,tmpnum2,Se,Rsite,tmpnum))
    Rsitef = Rsite
    newlabel <- c(newlabel,Rsite)
    drawlim <- c(drawlim,c((Rsite-Se):(Rsite+Se)))
  }}
}
    if(argv$Mut != 'Mut'){
    if(length(newlabel)==0){
      newlabel = '-'
    }
    SumNdb[which(SumNdb$CGG_repeat==j),sprintf('%s_number',argv$motif)] = length(newlabel)
    SumNdb[which(SumNdb$CGG_repeat==j),sprintf('%s_interruption',argv$motif)] = paste(newlabel,collapse = ',')
}else{
    SumNdb[which(SumNdb$CGG_repeat==j),sprintf('%s_number',argv$motif)] = '-'
    SumNdb[which(SumNdb$CGG_repeat==j),sprintf('%s_interruption',argv$motif)] = '-'

  }
   }
     
  if(exists('drawlim')){
    drawfile = tmpafile[tmpafile$STR_num %in% drawlim]
  if(argv$lowper == 0){
    p1 <- ggplot(data = tmpafile,aes(x=STR_num))+
    geom_bar(position = 'identity',col='red')+
    xlab('AGG position')+
    ylab('Reads number')+
    theme_classic()+
    theme(axis.line = element_line(colour = 'black'),
          axis.title = element_text(size = 15),
          axis.text = element_text(face = "bold"))+
    theme(axis.text.x = element_text(angle = 90))+
    scale_y_continuous(limits=c(0,max_y),expand = c(0.01,1.2))
  }else{ 
  p1 <- ggplot(data = drawfile,aes(x=STR_num))+
  geom_bar(position = 'identity',col='red')+
  xlab('AGG position')+
  ylab('Reads number')+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black'),
        axis.title = element_text(size = 15),
        axis.text = element_text(face = "bold"))+
  theme(axis.text.x = element_text(angle = 90))+
  scale_y_continuous(limits =c(0,max_y),expand = c(0.01,1.2))}
  if(newlabel != '-'){
    p1 <- p1+scale_x_continuous(breaks =newlabel,labels = newlabel,limits = c(0,240))
  }else{
    p1 <- p1+scale_x_continuous(labels=NULL)
  }
suppressMessages(ggsave(sprintf('%s/%s_%s_All_copy_%s_plot_distribute.pdf',argv$output,argv$name,argv$motif,j),plot=p1,width = 15))

}
  
}

if(file.exists('STR_statis_new.tsv')){
	write.table(SumNdb,'STR_statis_new.tsv',quote=F,row.names=F,sep='\t',append=T,col.names=F)
}else{
	write.table(SumNdb,'STR_statis_new.tsv',quote=F,row.names=F,sep='\t')
}


  


