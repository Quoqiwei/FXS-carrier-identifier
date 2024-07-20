library(data.table)
library(ggplot2)
library(argparser)
wkdir=getwd()
p <- arg_parser("Statistically analyze and plot STR results")
p <- add_argument(p, "--name", short = '-n',help='sample name',default = 'Sample')
p <- add_argument(p, "--input", short = '-i',help='str summary file',default = 'STR_summary_hcov.tsv')
p <- add_argument(p, "--output", short = '-o',help='str summary file',default = wkdir)
p <- add_argument(p, "--strnum", short = '-S',help='maximum str number',default = 10 )
p <- add_argument(p, "--lowper", short = '-l',help='threshold of premutation and full mutation allele',default = 0 )
p <- add_argument(p, "--lowper53", short = '-l53',help= 'threshold of wildtype allele',default = 10000)
p <- add_argument(p, "--thrvalue", short = '-tv',help='threads value',default = 53 )
p <- add_argument(p, "--motif", short = '-m',help='motif',default = 'AGG' )
p <- add_argument(p, "--SE", short = '-se',help='SE',default = '/home/gqw/STR_soft/STR_test/script/SE.tsv')
argv <- parse_args(p)

num = args
afile <- fread(argv$input,col.names = c('ID','STR_len','STR_num','qvalue','tarlen','cov','sam'))
afile$STR_num <- round(afile$STR_num)
afile <- na.omit(afile)
afile <- afile[afile$STR_num > 0,]
defile <- density.default(afile$STR_num,bw=0.5,kernel = 'gaussian',n = 1024)
sefile = argv$SE
defile$x[which.max(defile$y)]
getSE <- function(sefile,site){
  tmpdb <- fread(sefile)
  tmpdb$lab <- apply(tmpdb,1,function(x){
    ttr <- x['V1'][[1]][1]
    if(grepl('-',ttr)){
      trange <- c(strsplit(ttr,'-')[[1]][1]:strsplit(ttr,'-')[[1]][2])
      if(floor(site) %in% trange){
        return(x['V2'][[1]][1])
      }else{
        return('notin')
      }
    }else{
      trange = as.numeric(gsub('\\+','',ttr))
      if(site > trange){
        return(x['V2'][[1]][1])
      }else{
        return('notin')
      }
    }
  })
  tmpdb <- tmpdb[tmpdb$lab!='notin']
  tts <- tmpdb$lab[[1]][1]
  tts <- as.numeric(tts)
  if(tts >= 1){
    ts = tts
  }else{
    ts = ceiling(site*tts)
  }
  return(ts)
}
getmode <- function(ndb,tmpx) {
  tse <- getSE(sefile,tmpx)
  topldb <-  ndb[(ndb$STR_num>tmpx-tse)&(ndb$STR_num<tmpx+tse),]
  v = topldb$STR_num
  uniqv <- unique(v)
  maxx <- uniqv[which.max(tabulate(match(v, uniqv)))]
  tse = round(maxx/10)
  topldb <-  ndb[(ndb$STR_num>maxx-tse)&(ndb$STR_num<maxx+tse),]
  nrow(topldb)
}
getmodex <- function(ndb,tmpx) {
  tmpx = ceiling(tmpx)
  tse = ceiling(tmpx/10)
  topldb <-  ndb[(ndb$STR_num>tmpx-tse)&(ndb$STR_num<tmpx+tse),]
  v = topldb$STR_num
  uniqv <- unique(v)
  maxx <- uniqv[which.max(tabulate(match(v, uniqv)))]
  nrow(topldb)
}

testfun <- function(x,nnu){
  list1 = 0
  list2 = 0
  topl <- which.max(x$y)
  toplx <- x$x[topl]
  maxnum <- getmode(afile,toplx)
  for(i in c(2:(length(x$y)-1))){
    tmpx = ceiling(x$x[i])
    tmpy = getmode(afile,tmpx)
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
        if(tmpy >lowper){
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
Strnum <- testfun(defile,argv$strnum)
label1 <- round(Strnum$list2)
Strnum1 <- testfun(defile,(argv$strnum)*2)
label2 <- round(Strnum1$list2)
getbig <- function(tdata,num) {
  return(tdata)
}
Rsitef = 0
newlabel = c()
SumNdb = data.frame(matrix(ncol=6,nrow = 0))
names(SumNdb) <- c('Sample_name','Allele','CGG_repeat',sprintf('%s_number',argv$motif),sprintf('%s_interruption',argv$motif),'Alarm')
SumNdb[1,'Sample_name'] <- argv$name
SumNdb[is.na(SumNdb)] <- '-'
reads_number <- list()
for(i in Strnum[order(Strnum$list2),]$list2){
  Rsite_old = 0
  Rsite = round(i)
  while (Rsite != Rsite_old) {
    Se <- getSE(sefile,Rsite)
    tmpfile = afile[afile$STR_num>Rsite-Se,]
    tmpfile = tmpfile[tmpfile$STR_num<Rsite+Se,]
    tdb = data.frame(table(tmpfile$STR_num))
    Rsite_old = Rsite
    if(nrow(tdb)>0){
    	Rsite = as.numeric(as.character(tdb[which.max(tdb$Freq),'Var1']))
    }else{
    	Rsite = Rsite
    }
    
  }
  
  tmpnum = nrow(tmpfile)
  tmpnum2 = nrow(afile[afile$STR_num == Rsite,])
  reads_number[[(length(reads_number)+1)]] <- tmpnum2
  if(Rsite != Rsitef){
    if(argv$lowper != 0){
      if(Rsite < 53 && tmpnum > argv$lowper53){
        print(sprintf('%s reads are %s CGG repeats. The number of reads within ±%s of %s CGG repeats are %s.',tmpnum2,Rsite,Se,Rsite,tmpnum))
        Rsitef = Rsite
        newlabel <- c(newlabel,Rsite)
      }else{
        if(Rsite>= 53){
           print(sprintf('%s reads are %s CGG repeats. The number of reads within ±%s of %s CGG repeats are %s.',tmpnum2,Rsite,Se,Rsite,tmpnum))
        Rsitef = Rsite
        newlabel <- c(newlabel,Rsite)
        }
      }
    }else{
      print(sprintf('%s reads are %s CGG repeats. The number of reads within ±%s of %s CGG repeats are %s.',tmpnum2,Rsite,Se,Rsite,tmpnum))
      Rsitef = Rsite
      newlabel <- c(newlabel,Rsite)
    }
  }
}
if(is.null(newlabel)){
	print('Threshold is too high to output the result.')
}else{

SumNdb[c(1:length(newlabel)),'CGG_repeat'] <- newlabel
SumNdb$Allele <- c(1:length(newlabel))
SumNdb$Alarm <- apply(SumNdb,1,function(x){
  if(as.numeric(x['CGG_repeat'])>=argv$thrvalue){
    return('carrier')
  }else{
    return('-')
  }
  })}
SumNdb$Sample_name <- argv$name
SumNdb[is.na(SumNdb)] <- '-'
if(file.exists('STR_statis.tsv')){
 write.table(SumNdb,'STR_statis.tsv',quote=F,row.names=F,sep='\t',append=T,col.names=F)
}else{
 write.table(SumNdb,'STR_statis.tsv',quote=F,row.names=F,sep='\t')
}
#----Plot------
max_y <- max(unlist(lapply(reads_number,function(x){return(as.numeric(x))})))
max_y = round(max_y/0.8)
min_y = round(tmpnum2/0.8)
p1 <- ggplot(data = afile,aes(x=STR_num))+
  geom_bar(position = 'identity',col='red')+
  geom_line(stat = 'count',col='#0000CD')+
  xlab('CGG repeat number')+
  ylab('Reads number')+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black'),
        axis.title = element_text(size = 15),
        axis.text = element_text(face = "bold"))+
  scale_x_continuous(limits = c(0,240),breaks =newlabel,labels = newlabel)+
  theme(axis.text.x = element_text(angle = 90))+
  scale_y_continuous(limits = c(0,max_y), expand = c(0.01,1.2))

if (round(Strnum$list2[1])!= Rsite && (newlabel[which.max(newlabel)]-newlabel[which.min(newlabel)])>21){
  bfile <- afile[afile$STR_num > round(Rsite-21),]
  bfile <- bfile[bfile$STR_num < round(Rsite+21),]
  p2 <- ggplot(data = bfile,aes(x=STR_num))+
    geom_bar(position = 'identity',col='red')+
    geom_line(stat = 'count',col='#0000CD')+
    xlab('CGG repeat number')+
    ylab('Reads number')+
    theme(axis.line = element_line(colour = 'black'),
          axis.title = element_text(size = 15),
          axis.text = element_text(face = "bold"))+
    scale_x_continuous(limits = c(Rsite - 21,Rsite + 21),breaks =newlabel,labels =newlabel)+
    theme(axis.text.x = element_text(angle = 90))+
    scale_y_continuous(limits = c(0, min_y), expand = c(0.01,1.2))
}else{
   p2 <- ggplot(data = afile,aes(x=STR_num))+
         geom_bar(position = 'identity',col='red')+
         geom_line(stat = 'count',col='#0000CD')+
         xlab('CGG repeat number')+
         ylab('Reads number')+
         theme(axis.line = element_line(colour = 'black'),
               axis.title = element_text(size = 15),
               axis.text = element_text(face = "bold"))+
         scale_x_continuous(limits = c(0,240),breaks =newlabel,labels =newlabel)+
         theme(axis.text.x = element_text(angle = 90))+
         scale_y_continuous(limits = c(0, max_y), expand = c(0.01,1.2))
}

cfile <- afile[afile$STR_num %in% c(0:250)]
p3 <- ggplot(data = cfile,aes(x=STR_num))+
   geom_bar(position = 'identity',col='red')+
   geom_line(stat = 'count',col='#0000CD')+
   xlab('CGG repeat number')+
   ylab('Reads number')+
   theme_classic()+
   theme(axis.line = element_line(colour = 'black'),
	       axis.title = element_text(size = 15),
	       axis.text = element_text(face = "bold"))+
   scale_x_continuous(limits = c(0,240),breaks = newlabel,labels = newlabel)+
   theme(axis.text.x = element_text(angle = 90))+
   scale_y_continuous(limits = c(0,max_y),expand = c(0.01,1.2))
ggsave(sprintf('%s/%s_All_plot_distribute.pdf',argv$output,argv$name),plot=p1,width = 15)
ggsave(sprintf('%s/%s_BigSTR_plot_distribute.pdf',argv$output,argv$name),plot=p2,width = 15)
ggsave(sprintf('%s/%s_minSTR_plot_distribute.pdf',argv$output,argv$name),plot=p3,width = 15)
  


