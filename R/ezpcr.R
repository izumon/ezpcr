
 #' ezread
 #'
 #' This function reads data exported by QuantStudio.
 #' QuantStudioからエクスポートされたデータを読み込む
 #'
 #' @param dir folder/directory.
 #' @param skip skip rows.
 #' @param SMAPLENAME column of sample names.
 #' @param TARGETNAME column of target names.
 #' @param Ct column of Ct values.
 #' @return dataframe.
#' @export
ezRead<-function(dir="./",skip=40,SAMPLENAME='Sample Name',TARGETNAME='Target Name',Ct="CT")
{
    library(dplyr)
    library(readxl)
    files<-list.files(dir,pattern="*.csv|*.txt|*.xlsx|*.xls")
    cat(sprintf("files \033[32m %s \033[0m have been found",paste(files,collapse=",")))
    tables<-lapply(files,function(f)
    {
        A0<-NULL
        if(endsWith(f,"txt")||endsWith(f,".csv"))
        {
            A0<-read.csv(paste(dir,f,sep="/"),header=T,sep=",",skip=skip)
        }
        if(endsWith(f,"xlsx")){
            A0<-read_xlsx(paste(dir,f,sep="/"),skip=skip)
        }
        if(endsWith(f,"xls")){
            A0<-read_xls(paste(dir,f,sep="/"),skip=skip)
        }
        if(is.null(A0)){
            print(sprintf("%s was not found.",paste(dir,f,sep="/")))
        }
    #A0<-subset(A0,subset=!is.na(Omit))
    A0<-as.data.frame(A0)
    A0<-data.frame(Samples=as.vector(A0[,SAMPLENAME]),Targets=as.vector(A0[,TARGETNAME]),Ct=as.numeric(as.vector(A0[,Ct])))
    A0<-A0[!is.na(A0$Samples),]
    return(as.data.frame(A0))
    })
    atable<-bind_rows(tables)
    return(atable)
}#end function ezread


 #' ezcalc
 #'
 #' This function calculates dCt ddCt RQ.
 #' Ct値から比較定量値Relative Quantity(RQ)を求める
 #' @param df dataframe obtained from ezread.
 #' @param internalControl a sample name of biological-control.
 #' @param internalControl a gene of internal-control.
 #' @param CtMax if Ct Value is undetermined(NA), this value sets to CtMax.
 #' @param CtTh if Ct Value is more than CtTh, this value sets to CtMax.
 #' @return dataframe.
#' @export
ezCalc<-function(df,biologicalControl,internalControl="GAPDH",CtMax=40,CtTh=40)
{
library(dplyr)
    if(length( subset(df,subset=Targets==biologicalControl))==0)
    {
        print(sprintf(" %s is not containd in this data!"))
    }

#Ctを数値化
    df<-df %>% dplyr::mutate(Ct=as.numeric(Ct)) %>% dplyr::mutate(Ct=if_else(is.na(Ct) | Ct>CtTh,CtMax,Ct))

#CtMean
    x3<-df%>% group_by(Samples,Targets) %>% dplyr::mutate(CtMean=mean(Ct,na.rm=TRUE)) %>% ungroup

#calc internal control mean
    icontrol <- tapply(x3$CtMean,list(x3$Samples,x3$Targets),mean)
    icontrol<- icontrol[,internalControl]

#dCt
    x3$dCt<-x3$Ct
    for (x in names(icontrol))
    {
           x3<-x3%>%mutate(dCt=if_else(Samples==x,dCt-icontrol[x],dCt))
           #filt<-which(x3$Samples==x,);x3[filt,"dCt"]<-x3[filt,"dCt"]-icontrol[x] #で最初に変更する行を抽出するか
           #x3[x3$Samples==x,] %<>% mutate(dCt=dCt-icontrol[x]) #magrittrを使ったパイプの方がスマートでは
    }

#dCtMean
    x3<-x3%>%
        group_by(Samples,Targets)%>%
        mutate(dCtMean=mean(dCt))%>% 
        ungroup()

#ddCt
    bcons<-paste(unique(grep(biologicalControl,x3$Samples,value=TRUE)),collapse=", ")
    cat(sprintf("Sample name \033[31m %s\033[m was used as biological control \r\n", bcons))
    x3[grepl(biologicalControl,x3$Samples),"Samples"]<-biologicalControl
    bcontrol<-tapply(x3$dCt,list(x3$Samples,x3$Targets),mean)
    bcontrol <- bcontrol[biologicalControl,]
    x3$ddCt<-x3$dCtMean
    for(x in names(bcontrol))
    {
        x3<-x3 %>% transform(ddCt=if_else(Targets==x,dCt-bcontrol[x],ddCt))
    }

#ddCtMean
    x3<-x3%>% group_by(Samples,Targets) %>% mutate(ddCtMean=mean(ddCt)) 
    x3<-x3%>% mutate(RQ=2^(-ddCt))
    x3<-x3%>% mutate(RQMEAN=mean(RQ)) %>% ungroup() #%>% as.data.frame()

    x3<-x3%>% mutate(rdCtMean=-dCtMean)
    return(as.data.frame(x3))

}#end function ezcalc



 #' ezgraph
 #'
 #' This function creates bar plots.
 #' ezpcrで出力されたデータを棒グラフにする.
 #' @param data dataframe that contains Samples,Targets,RQ,Signif
 #' @param samplenames sample names that is contained in dataframe$Samples.
 #' @param dot TRUE/FALSE, indicationg whether the plot includes a dotplot .
 #' @param linewidth line thickness of axis lines and bar plot lines.
 #' @param textSize text size of tick labels.
 #' @param labelSize text size of sample labels.
 #' @param titleSize text size of gene name.
 #' @param legendPosition position of legend. "none"=no legends "right","top","bottom","left"=legend positions
 #' @param genes selected genes that are contained in dataframe$Targets.
 #' @param titleSize significance ex) list(c("sampleA","sampleB")).
 #' @param color Bar colors. the number of bar colors has to be greater than number of samples.
 #' @param dotsize Size of dots. if dot option is TRUE, this option is applied to dot plot.
 #' @param newline line feed characters. the characters separate sample name and add new line character. 
 #' @param y_extension scale that extends max value of y axis. 
 #' @return list(plot); function names(returned value) returns gene names.
 #' @examples p <- ezGraph(dataframe,samplenames=c("S1","S2","S3"),dot=FALSE,genes=c("gene1","gene2"),color=c("red","blue","white"))
 #' plot(p)
 #' @export
ezGraph<-function(data,samplenames=NULL,dot=FALSE,linewidth=2,textSize=22,labelSize=26,titleSize=32,legendPosition="none",genes=NULL,signifs=list(),color=c(),dotsize=3,newline=" ",y_extension=1.05)
{
    library(ggplot2)
    library(dplyr)
    library(ggpubr)

    if(is.null(genes)){
        genes<-unique(data$Targets)
    }
    #samplenamesに含まれるサンプルのデータのみを抽出
    if(is.null(samplenames))
    {
        samplenames<-unique(data$Samples)
    }
        datas<-lapply(samplenames,function(x) {
            d<-data[grepl(x,data$Samples),]
            cat(sprintf("\033[31m%s\033[0m is renamed to \033[32m%s\033[0m\r\n",paste(unique(d$Samples),collapse=", "),x))
            d[grepl(x,d$Samples),"Samples"]<-x
            return(d)
        })#end lapply samplenames

        data0<-bind_rows(datas)
        data0$Samples<-gsub(newline,"\r\n",data0$Samples)
        samplenames <- gsub(newline,"\r\n",samplenames)
        data0$Samples<-factor(data0$Samples,levels=samplenames)

    PS<-lapply(genes,function(g)
        {
            data1<-subset(data0,subset=Targets==g)
	p <- ggplot(data = data1 , aes(x = Samples, y = RQ, fill= Samples)) + #aesで使用するパラメータを指定

    #plot ===
	stat_summary(geom="bar", fun=mean, color="black", linewidth=linewidth,width=0.7 )+ #自動的に平均化した棒グラフを作ってくれる stat_summary
    #group化する際には barとerrorbarにpositoin = position_dodge(0.5)などを付けること
	stat_summary(geom="errorbar",fun.data=mean_sdl,fun.args = list(mult = 1), width=0.5 ,size=linewidth)+ #エラーバーも自動で計算してくれる 標準誤差(mean_se)を使用 標準偏差は(mean_sdl,fun.args = list(mult=1)) あるいはggpubrのmean_sd
	scale_y_continuous(limit=c(0,max(data1$RQ)*y_extension) , expand=c(0,0))+ #x軸の最小値を０に固定 NAをmax(data)*xにすると、最大値を拡張できる
	ggtitle(g)+ 
	theme_classic()+ #シンプルなデザインに変更
	theme(
        plot.title = element_text(hjust=0.5), #theme : 軸の太さなどの細かい点を指定
	axis.title.x = element_blank(),
    axis.title.y = element_text( size = labelSize, vjust = 2),
	axis.text.x = element_text(size=textSize, color="black" ),
    axis.text.y = element_text(size=textSize, color="black" ),

	title=element_text(size=titleSize),

	axis.line = element_line(linewidth =linewidth),
	axis.ticks = element_line(linewidth=linewidth),
    axis.ticks.x = element_blank(),
	axis.ticks.length.y = unit(2,"mm"),
    legend.position = legendPosition #without legend
    )

    #optional
    if(dot==TRUE)
    {
      p <-  last_plot()+geom_jitter(color = "black", fill = "red", size = dotsize,shape = 21 , width = 0.3)  #position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.2)# groupingのときにつける width,fillは除外すること
    }
 
    #引数で有意差を追加
    if(length(signifs)!=0){
    #add text ===
  p<-last_plot()+geom_signif(
    comparisons = list(signifs),  # 比較するペア
    test = "wilcox.test",                  # t検定を使用
    map_signif_level = TRUE           # p<0.05,* / p<0.01,** / p<0.001,*** に自動変換
  )}

    if(NROW(color)>=NROW(unique(data1$Samples))){
        #settings ===
        p<-last_plot()+scale_fill_manual(values=color) #色指定
    }else if(NROW(color)<NROW(unique(data1$Samples)) && NROW(color)!=0){
        print("color list is less than number of samples")
        print(sprintf("color=%s, samples=%s",NROW(color),NROW(unique(data1$Samples))))
    }

    #signifがすでに入力されているとき有意差を表示
    if("signif" %in% colnames(data))
    {
        p<-last_plot()+stat_summary(fun="max", geom="text",aes(label=signif),vjust=0.5,size=textSize*1.5,color="black",fontface="bold")
    }

    return(p)

    })#end lapply genes
    names(PS)<-genes
    return(PS)
} #end function ezgraph


#' extract sample names from the dataframe generated by ezpcr
#' @export
lS<-function(df)
{
    return(unique(df$Samples))
}

#' extract sample names from the dataframe generated by ezpcr
#' @export
lSamples<-function(df)
{
    return(unique(df$Samples))
}

#' extract gene targets from the dataframe generated by ezpcr
#' @export
lT<-function(df)
{
    return(unique(df$Targets))
}

#' extract gene targets from the dataframe generated by ezpcr
#' @export
lTargets<-function(df)
{
    return(unique(df$Targets))
}


#' extract data from the dataframe generated by ezpcr
#' @export
ext<-function(data,Sample="",Target="")
{
    return(data[grepl(Sample,data$Samples)&grepl(Target,data$Targets),])
}

#' write plots as png
#' @export
ezPng<-function(plots,plotname="plot",width=5,height=5,dpi=350)
{
    for(w in names(P))
    {
        ggsave(sprintf("./%s-%s.png",plotname,w),plot=P[[w]],device=png,width=width,height=height,dpi=dpi,bg="white")
    }
}


#' write plots as svg
#' @export
ezSvg<-function(plots,plotname="plot",width=5,height=5,dpi=350)
{
    if(!("svglite" %in% installed.packages()))
    {
        install.packages("svglite")
    }
    library("svglite")

    for(w in names(P))
    {
        ggsave(sprintf("./%s-%s.svg",plotname,w),plot=P[[w]],device=svglite,width=width,height=height,dpi=dpi,bg="white")
    }
}


#' checking outlier depends on Turkey's IQR(Interquartile Range) rule. 
#' 四分位範囲による外れ値検定
#' @param data frame which contains Samples,Targets and RQ.
#' @export
check_outlier <- function( data ){
    library(dplyr)
    data <- data %>% group_by(Samples,Targets) %>% mutate(outlier = (function(rq){
            Q3 <- quantile(rq,0.75)
            Q1 <- quantile(rq,0.25)
            IQR <- Q3-Q1
            lower <- Q1 -1.5*IQR
            upper <- Q3 + 1.5*IQR
            return(if_else((rq >= lower & rq<=upper),FALSE,TRUE))
              })(RQ)) 
}

#' removing outlier depends on Turkey's IQR(Interquartile Range) rule. 
#' 四分位範囲による外れ値検定
#' @param data frame which contains Samples,Targets and RQ.
#' @export
remove_outlier <- function(data){
    if(!("outlier" %in% colnames(data)))
    {
        data<-check_outlier(data)
    }
    return(subset(data,subset=outlier==FALSE))
}

#' t.test 
#' t検定. 比較対象が複数あるときはBH法による多重検定
#' @param data data frame which contains Samples,Targets and dCt. @param control control name which is contaied in column Samples.
#' @param samples sample names which is contaied in column Samples.
#' @export
ttest <- function(data, control,samples){

    con <- data[grepl(control,data$Samples),]
    cat(sprintf("\033[33m %s \033[0m is used as control\r\n",paste(unique(con,collapse=","))))

    if(NROW(samples)==1){ pair<-TRUE}else{pair<-FALSE}

    #begin lapply ===========================
    results<-lapply(data$Targets,function(tg){
        tdata <- subset(data,subset=Targets==tg)
        #get control data
        con <- data[grepl(control,data$Samples),]
        for(w in samples)
        {
            tdata<-mutate( p_val=if_else(grepl(w,Samples),t.test(con$dCt,dCt)$p.value,1))
            tdata<-mutate( p_adjust=if_else(grepl(w,Samples),pairwize.t.test(con$dCt,dCt,p.adjust.method="BH")$p.value,1))
        }
    })%>% bind_rows() #end lapply============

    # 列signif に * ** を追加
    if(pair==TRUE)
    {
        results<-mutate( signif = if_else(p_val<0.05,"*",""),
            signif=if_else(p_val<0.01,"**",signif))
    }else
    {
        results<-mutate( signif = if_else(p_adjust<0.05,"*",""),
            signif=if_else(p_adjust<0.01,"**",signif))
    }

    return(result)
}


#' wilcox.test 
#' wilcox検定 ウィルコクソンの順位和検定.比較対象が複数あるときはBH法による多重検定
#' @param data data frame which contains Samples,Targets and dCt.
#' @param control control name which is contaied in column Samples.
#' @param samples sample names which is contaied in column Samples.
#' @export
wilcoxtest <- function(data, control,samples){

    con <- data[grepl(control,data$Samples),]
    cat(sprintf("\033[33m %s \033[0m is used as control\r\n",paste(unique(con,collapse=","))))

    if(NROW(samples)==1){ pair<-TRUE}else{pair<-FALSE}

    #begin lapply ===========================
    results<-lapply(data$Targets,function(tg){
        tdata <- subset(data,subset=Targets==tg)
        #get control data
        con <- data[grepl(control,data$Samples),]
        for(w in samples)
        {
            tdata<-mutate( p_val=if_else(grepl(w,Samples),wilcox.test(con$dCt,dCt)$p.value,1))
            tdata<-mutate( p_adjust=if_else(grepl(w,Samples),pairwize.wilcox.test(con$dCt,dCt,p.adjust.method="BH")$p.value,1))
        }
    })%>% bind_rows() #end lapply============

    # 列signif に * ** を追加
    if(pair==TRUE)
    {
        results<-mutate( signif = if_else(p_val<0.05,"*",""),
            signif=if_else(p_val<0.01,"**",signif))
    }else
    {
        results<-mutate( signif = if_else(p_adjust<0.05,"*",""),
            signif=if_else(p_adjust<0.01,"**",signif))
    }

    return(result)
}
