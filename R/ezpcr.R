# update予定
# biologicalControlが複数サンプルある場合、それらすべての平均をコントロールとするが、サンプル名は変更しない仕様にする

#' ezUpdate
#' ezpcrをアップデートする
#' @export
ezUpdate<-function()
{
    if(!"pak" %in% installed.packages())
    {
        install.packages("pak")
    }
    library(devtools)
    library(pak)
    tryCatch(
    pak::pak("izumon/ezpcr"),
    error = function(e) e,
    finally = devtools::install_github("izumon/ezpcr")
    )
    return()
}

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
    #cat(sprintf("files \033[32m %s \033[0m have been found",paste(files,collapse=",")))
    counter <<- 1
    tables<-lapply(files,function(f)
    {
        print(sprintf("file: %s is loading", f))
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
    A0<-A0[!is.na(A0[,SAMPLENAME]),]
    A0<-data.frame(Samples=as.vector(A0[,SAMPLENAME]),Targets=as.vector(A0[,TARGETNAME]),Ct=as.numeric(as.vector(A0[,Ct])),ID=counter)
    A0<-A0[!is.na(A0$Samples),]

    counter <<- counter +1
    return(as.data.frame(A0))
    })#end lapply
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
 #' @param CtMean and dCt calculation separated by Sample,Target,ID.
 #' @param technical TRUE...technical replicate is remained / FALSE...technical control is removed
 #' @param RQcontrolIsOne TRUE... RQ of the biological control is set to 1.
 #' @return dataframe.
#' @export
ezCalc<-function(df,biologicalControl,internalControl="GAPDH",CtMax=40,CtTh=40,ID=TRUE,technical=TRUE,RQcontrolIsOne=TRUE)
{
library(dplyr)
    if(NROW(df[grepl(internalControl,df$Targets),]) ==0 )
    {
        print(sprintf("Internal control %s is not containd in this data!"))
        return(NULL)
    }

    if(NROW(df[grepl(biologicalControl,df$Samples),]) ==0 )
    {
        print(sprintf("Biological control %s is not containd in this data!"))
        return(NULL)
    }

    if(!("ID" %in% colnames(df))){df$ID<-0}

#Ctを数値化
    df<-df %>% dplyr::mutate(Ct=as.numeric(Ct)) %>% dplyr::mutate(Ct=if_else(is.na(Ct) | Ct>CtTh,CtMax,Ct))

#internal controlを保存
    df$internalControl<-internalControl

#IDを付加
    df$Samples2<-paste(df$Samples,df$ID,sep="---")

#CtMean
    x3<-df%>% group_by(Samples2,Targets,ID) %>% dplyr::mutate(CtMean=mean(Ct ,na.rm=TRUE)) %>% ungroup

#calc internal control mean
    icontrol <- tapply(x3$CtMean,list(x3$Samples2,x3$Targets),mean)
    if(internalControl %in% rownames(icontrol))
    {
        icontrol<-t(icontrol)
    }

    icon <- icontrol[,grepl(internalControl,colnames(icontrol))]
        if(NROW(unlist(strsplit(internalControl,"\\|")))>1)
        {
            icontrol<-  icon %>% apply(1,mean)
        }else{
            icontrol<- icon %>% mean
        }

#dCt
    x3$dCt<-x3$Ct
    for (x in names(icontrol))
    {
           x3<-x3%>%mutate(dCt=if_else(Samples2 == x, dCt-icontrol[x], dCt))
           #filt<-which(x3$Samples==x,);x3[filt,"dCt"]<-x3[filt,"dCt"]-icontrol[x] #で最初に変更する行を抽出するか
           #x3[x3$Samples==x,] %<>% mutate(dCt=dCt-icontrol[x]) #magrittrを使ったパイプの方がスマートでは
    }

#IDを除去
    x3$Samples2<-NULL

#dCtMean
    x3<-x3%>%
        group_by(Samples,Targets)%>%
        mutate(dCtMean=mean(dCt))%>% 
        ungroup()

#biological controlを保存
    x3$biologicalControl <- x3$Samples

#ddCt
    bcons<-paste(unique(grep(biologicalControl,x3$Samples,value=TRUE)),collapse=", ")
    cat(sprintf("Sample name \033[31m %s\033[m was used as biological control \r\n", bcons))
    x3[grepl(biologicalControl,x3$Samples),"Samples"]<-biologicalControl
    bcontrol<-tapply(x3$dCt, list(x3$Samples,x3$Targets),mean)
    bcontrol <- bcontrol[biologicalControl,]
    x3$ddCt<-x3$dCt
    for(x in names(bcontrol))
    {
        x3<-x3 %>% transform(ddCt=if_else(Targets==x, dCt-bcontrol[x],ddCt))
    }

#ddCtMean
    x3<-x3%>% group_by(Samples,Targets) %>% mutate(ddCtMean=mean(ddCt)) 
    x3<-x3%>% mutate(RQ=2^(-ddCt))
    if(RQcontrolIsOne==TRUE){
        x3[x3$Samples==biologicalControl,"RQ"]=1
    }
    x3<-x3%>% mutate(RQMEAN=mean(RQ)) %>% ungroup() #%>% as.data.frame()

    #control を 1にする

    #Sample nameを元に戻し、biologicalControl名を記録
    x3$Samples<-x3$biologicalControl
    x3$biologicalControl <- biologicalControl

    if(technical==FALSE)
    {
        x3<-ezShrink(x3)
        print("technical replicates were grouped together.")
    }


    #x3<-x3%>% mutate(rdCtMean=-dCtMean)
    return(as.data.frame(x3))

}#end function ezcalc


 #' ezgraph
 #'
 #' This function creates bar plots.
 #' ezpcrで出力されたデータを棒グラフにする.
 #' @param data dataframe that contains Samples,Targets,RQ,Signif
 #' @param samplenames sample names that is contained in dataframe$Samples.
 #' @param targets gene names that is contained in dataframe$Targets.
 #' @param dot TRUE/FALSE, indicationg whether the plot includes a dotplot .
 #' @param LOG10 TRUE/FALSE, create the plot on a log scale .
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
 #' @param technical TRUE technical replicate is remained / FALSE technical control is removed
 #' @param significant_column text column that is put above significance bar.
 #' @param txtNotSignif text if this data is not significant, txtNotSignif is put above significance bar.
 #' @param signifSize text size of text of significant.
 #' @return list(plot); function names(returned value) returns gene names.
 #' @examples p <- ezGraph(dataframe,samplenames=c("S1","S2","S3"),dot=FALSE,genes=c("gene1","gene2"),color=c("red","blue","white"))
 #' plot(p)
 #' @export
ezGraph<-function(data,samplenames=NULL,targets=c(),dot=FALSE,linewidth=2,textSize=22,labelSize=26,titleSize=32,legendPosition="none",genes=NULL,signiflist=list(),color=c(),dotsize=3,newline=" ",y_extension=1.05,controlIsOne=TRUE,technical=TRUE,significant_column="signif",textNotSignif="n.s.",signifSize=10,LOG10=FALSE)
{
    library(ggplot2)
    library(dplyr)
    library(ggpubr)
    library(scales)

    if(controlIsOne==TRUE)
    {
        bc <- unique(data$biologicalControl)
        data[data$Samples==bc,"RQ"]<-1
    }
    if(technical==FALSE)
    {
        data<-ezShrink(data)
    }

custom_scale <- function(x) {
      index_zero <- which(x == 0)
      under_10k <- which(x<10**4)
      label2 <- x
      label <- scientific_format()(x)
      label <- gsub( "e", " %*% 10^",label)
      label <- gsub("\\^\\+", "\\^",label)
      label[index_zero] <- "0"
      label[under_10k] <- label2[under_10k]
      parse(text=label)
}
    #ログスケールか否か
    if(LOG10==TRUE)
    {
        tran<-"log10"
    }else{
        tran<-"identity"
    }

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
            d[grepl(x,d$Samples),"Samples"] <- x
            return(d)
        })#end lapply samplenames

        data0<-bind_rows(datas)
        data0$Samples<-gsub(newline,"\r\n", data0$Samples)
        samplenames <- gsub(newline,"\r\n", samplenames)
        data0$Samples<-factor(data0$Samples, levels=samplenames)

        #グラフ化する遺伝子を選択
        if(NROW(targets)!=0){
            genes <- intersect(genes,targets)
        }
    PS<-lapply(genes,function(g)
        {
            data1<-subset(data0,subset=Targets==g)
	p <- ggplot(data = data1 , aes(x = Samples, y = RQ, fill= Samples)) + #aesで使用するパラメータを指定

    #plot ===
	stat_summary(geom="bar", fun=mean, color="black", linewidth=linewidth, width=0.7 )+ #自動的に平均化した棒グラフを作ってくれる stat_summary
    #group化する際には barとerrorbarにpositoin = position_dodge(0.5)などを付けること
	stat_summary(geom="errorbar",position="dodge",
               fun = mean, 
               fun.min = function(x) pmax(mean(x) - sd(x), 0), 
               fun.max = function(x) mean(x) + sd(x),
        width=0.5 ,size=linewidth)+ #エラーバーも自動で計算してくれる 標準誤差(mean_se)を使用 標準偏差は(mean_sdl,fun.args = list(mult=1)) あるいはggpubrのmean_sd
#	scale_y_continuous(limit=c(0,max(data1$RQ)*y_extension) , expand=c(0,0))+ #x軸の最小値を０に固定 NAをmax(data)*xにすると、最大値を拡張できる

	ggtitle(g)+ 
	theme_classic()+ #シンプルなデザインに変更
    coord_cartesian(ylim = c(0, NA))+ #エラーバーが切れるのを防ぐ
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

    if(LOG10==TRUE)
    {
    p<-last_plot()+ylab("Log10 Relative Quantity") + scale_y_continuous(limit=c(NA,max(log10(data$RQ))*y_extension) , trans=tran)

    }else
    {
        p<-last_plot()+ylab("Relative Quantity") + scale_y_continuous(limit=c(0,max(data1$RQ))*y_extension , expand=c(0,0), breaks = pretty_breaks(n=3), label=custom_scale, n.breaks=3, trans=tran)
    }

    #optional
    if(dot==TRUE)
    {
      p <-  last_plot()+geom_jitter(color = "black", fill = "red", size = dotsize,shape = 21 , width = 0.3)  #position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.2)# groupingのときにつける width,fillは除外すること
    }
 
    #引数で有意差を追加
    PSIGNIF<-NULL #有意差プロット
    if(length(signiflist)!=0){
        PSIGNIF <- ezSignif(data1,pairs=signiflist,textsize=labelSize,significant_column=significant_column,txtNotSignif=textNotSignif,signifSize=signifSize,LOG10=LOG10)

    #add text ===
#  p<-last_plot()+geom_signif(
#    aes(dCt),
#    comparisons = list(signifs),  # 比較するペア
#    test = "wilcox.test",                  # t検定を使用
#    map_signif_level = TRUE           # p<0.05,* / p<0.01,** / p<0.001,*** に自動変換
#  )
    }

    if(NROW(color)>=NROW(unique(data1$Samples))){
        #settings ===
        p<-last_plot()+scale_fill_manual(values=color) #色指定
    }else if(NROW(color)<NROW(unique(data1$Samples)) && NROW(color)!=0){
        print("color list is less than number of samples")
        print(sprintf("color=%s, samples=%s",NROW(color),NROW(unique(data1$Samples))))
    }

    if(!is.null(PSIGNIF))
    {
        p<- last_plot()+PSIGNIF
    }

    #signifがすでに入力されているとき有意差を表示
    if(("signif" %in% colnames(data)) && length(signiflist)==0 && is.null(PSIGNIF))
    {
        p<-last_plot()+stat_summary(fun="max", geom="text",aes(label=signif),vjust=0.5,size=textSize*1.1,color="black",fontface="bold")
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
#' @param plots list of bar plot out put by ezGraph
#' @param plotname text of filename
#' @param dir save directory. if plotname contains directory, dir is ignored.
#' @param width horizonal size of plot
#' @param height vertical size of plot
#' @param dpi dpi
#' @export
ezPng<-function(plots,plotname="plot",dir="./",width=5,height=5,dpi=350)
{
    for(w in names(P))
    {
        if(!grepl("/",plotname)){
            ggsave(sprintf("%s/%s-%s.png",dir,plotname,w),plot=P[[w]],device=png,width=width,height=height,dpi=dpi,bg="white")
        }else{
            ggsave(sprintf("%s-%s.png",plotname,w),plot=P[[w]],device=png,width=width,height=height,dpi=dpi,bg="white")
        }
    }
}


#' write plots as svg
#' @param plots list of bar plot out put by ezGraph
#' @param plotname text of filename
#' @param dir save directory. if plotname contains directory, dir is ignored.
#' @param width horizonal size of plot
#' @param height vertical size of plot
#' @param dpi dpi
#' @export
ezSvg<-function(plots,plotname="plot",dir="./",width=5,height=5,dpi=350)
{
    if(!("svglite" %in% installed.packages()))
    {
        install.packages("svglite")
    }
    library("svglite")

    for(w in names(P))
    {
        if(!grepl("/",plotname)){
            ggsave(sprintf("%s/%s-%s.svg",dir,plotname,w),plot=P[[w]],device=svglite,width=width,height=height,dpi=dpi,bg="white")
        }else{
            ggsave(sprintf("%s-%s.svg",plotname,w),plot=P[[w]],device=svglite,width=width,height=height,dpi=dpi,bg="white")
        }
    }
}


#' checking outlier depends on Turkey's IQR(Interquartile Range) rule. 
#' 四分位範囲による外れ値検定
#' @param data frame which contains Samples,Targets and RQ.
#' @param samples sample names.
#' @param using RQ to evaluate outlier.
#' @export
check_outlier <- function( data,samples=c(),RQ=FALSE ){
    library(dplyr)
    data$Samples2<-data$Samples
    if(NROW(samples)!=0)
    {
        for(w in samples)
        {
            print(sprintf("%s =  %s ",w, paste(unique(data[grepl(w,data$Samples),"Samples2"]),collapse=", ")))
            data[grepl(w,data$Samples),"Samples2"]<-w
        }
    }
    outlier<-function(rq){
    }
    if(RQ==TRUE){
    data <- data %>% group_by(Samples2,Targets) %>% mutate(outlier = (function(rq){
            Q3 <- quantile(rq,0.75)
            Q1 <- quantile(rq,0.25)
            IQR <- Q3-Q1
            lower <- Q1 -1.5*IQR 
            upper <- Q3 + 1.5*IQR
            return(if_else((rq >= lower & rq<=upper),FALSE,TRUE))
              })(RQ)) 
    }else{
    data <- data %>% group_by(Samples2,Targets) %>% mutate(outlier = (function(rq){
            Q3 <- quantile(rq,0.75)
            Q1 <- quantile(rq,0.25)
            IQR <- Q3-Q1
            lower <- Q1 -1.5*IQR 
            upper <- Q3 + 1.5*IQR
            return(if_else((rq >= lower & rq<=upper),FALSE,TRUE))
              })(dCt)) 
    }
    data$Samples2<-NULL

    return(as.data.frame(data))
}

#' removing outlier depends on Turkey's IQR(Interquartile Range) rule. 
#' 四分位範囲による外れ値検定
#' @param data frame which contains Samples,Targets and RQ.
#' @export
remove_outlier <- function(data,samples=c(),RQ=FALSE){
    if(!("outlier" %in% colnames(data)))
    {
        data<-check_outlier(data,samples,RQ)
    }
    print(as.data.frame(subset(data,subset=outlier==FALSE)))
    print("removed")
    return(as.data.frame(subset(data,subset=outlier==FALSE)))
}

#' t.test 
#' t検定. 比較対象が複数あるときはBH法による多重検定
#' @param data data frame which contains Samples,Targets and dCt. 
#' @param control control name which is contaied in column Samples.
#' @param samples sample names which is contaied in column Samples.
#' @export
ttest <- function(data, control,samples=c()){
    library(stats)

    if(!("p_value" %in% colnames(data)))
    {
            data$p_val<-1
            data$p_adjust<-1
    }
    con <- data[grepl(control,data$Samples),]
    cat(sprintf("\033[33m %s \033[0m is used as control\r\n",paste(unique(con$Samples,collapse=","))))

    if(NROW(samples)==1){ pair<-TRUE}else{pair<-FALSE}

    if(NROW(samples)==0)
    {
        samples <- setdiff(unique(data$Samples),unique(con$Samples))
    }

    #begin lapply ===========================
        results<-lapply(unique(data$Targets),function(tg){
        print(tg)
        tdata <- subset(data,subset=Targets==tg)
        #get control data
        con <- tdata[grepl(control,tdata$Samples),]
        tdata$compare<-tdata$Samples
        for(w in samples)
        {
            tdata[grepl(w,tdata$Samples),"compare"]<-w
            compdata<-tdata[tdata$compare==w,"dCt"]
            p_value <- t.test(con$dCt,compdata)$p.value
            tdata<- mutate(tdata, p_val=if_else(compare == w, p_value, p_val))
        }
        tdata$p_adjust <- p.adjust(tdata$p_val,method="BH")
        return(tdata)
    })%>% bind_rows() #end lapply============

    # 列signif に * ** を追加
    if(pair==TRUE || NROW(samples)>1)
    {
        results<-mutate(results, signif = if_else(p_val<0.05,"*",""),
            signif=if_else(p_val<0.01,"**",signif))
    }else
    {
        results<-mutate(results, signif = if_else(p_adjust<0.05,"*",""),
            signif=if_else(p_adjust<0.01,"**",signif))
    }

    results$compare<-NULL
    return(results)
}


#' wilcox.test 
#' wilcox検定 ウィルコクソンの順位和検定.比較対象が複数あるときはBH法による多重検定
#' @param data data frame which contains Samples,Targets and dCt.
#' @param control control name which is contaied in column Samples.
#' @param samples sample names which is contaied in column Samples.
#' @export
wilcoxtest <- function(data, control,samples){
    library(stats)

    if(!("p_value" %in% colnames(data)))
    {
            data$p_val<-1
            data$p_adjust<-1
    }

    con <- data[grepl(control,data$Samples),]
    cat(sprintf("\033[33m %s \033[0m is used as control\r\n",paste(unique(con$Samples,collapse=","))))

    if(NROW(samples)==1){ pair<-TRUE}else{pair<-FALSE}

    if(NROW(samples)==0)
    {
        samples <- setdiff(unique(data$Samples),unique(con$Samples))
    }

    #begin lapply ===========================
        results<-lapply(unique(data$Targets),function(tg){
        print(tg)
        tdata <- subset(data,subset=Targets==tg)
        #get control data
        con <- tdata[grepl(control,tdata$Samples),]

        tdata$compare<-tdata$Samples
        for(w in samples)
        {
            tdata[grepl(w,tdata$Samples),"compare"]<-w

            compdata<-tdata[tdata$compare==w,"dCt"]
            p_value <- t.test(con$dCt,compdata)$p.value
            tdata<- mutate(tdata, p_val=if_else(compare == w, p_value, p_val))
        }
        tdata$p_adjust <- p.adjust(tdata$p_val,method="BH")
        return(tdata)
    })%>% bind_rows() #end lapply============

    # 列signif に * ** を追加
    if(pair==TRUE || NROW(samples)>1 ) #サンプル数が複数だと強制的にpairwizeになる
    {
        results<-mutate(results, signif = if_else(p_val<0.05,"*",""),
            signif=if_else(p_val<0.01,"**",signif))
    }else
    {
        results<-mutate(results, signif = if_else(p_adjust<0.05,"*",""),
            signif=if_else(p_adjust<0.01,"**",signif))
    }


    results$compare<-NULL
    return(results)
}

#' join 
#' join two datas that are get from ezpcr::ezCalc.
#' ezpcrから得られたdataframeのdata1とdata2を結合する。
#' @export
ezJoin<-function(data1,data2)
{
    return(rbind(data1,data2))
}

#' heatmap 
#' heatmapを作成する。
#' @param data data frame which contains Samples,Targets and dCt.
#' @param samples sample names which is contaied in column Samples.
#' @param titlename title name of plot
#' @param legend add heat map legend
#' @param flip turn the data matrix.
#' @export
ezHeatmap<-function(data,samplenames=c(),titlename="",legend=FALSE,flip=FALSE)
{
    if(!("ComplexHeatmap" %in% installed.packages()))
    {
        install.packages(c("heatmap3","viridisLite","reshape2","RColorBrewer","gpar"))
    }
    if(!("genefilter" %in% installed.packages()))
    {
        BiocManager::install("genefilter")
    }
    if(NROW(samplenames) >0)
    {
        for(w in samplenames)
        {
            data[grepl(w,data$Samples),"Samples"]<-w
        }
    }
    library(heatmap3)
    library(genefilter)
    library(RColorBrewer)
    library(viridisLite)
    library(ComplexHeatmap)
    library(grid)

    mat<-tapply(data$dCt,list(data$Targets,data$Samples),mean)
    mat.z<-genescale(-mat,1,method="Z")
    if(flip==TRUE){mat.z<-t(mat.z)}
	P<-ComplexHeatmap::Heatmap(
	mat.z, #matrix 形式の data
	gap = unit(20, "mm"), #グループ間の境目 
	show_heatmap_legend = legend, #legendの非表示
	clustering_method_rows = "ward.D2",  #行のクラスタリングメソッド
	clustering_method_columns ="ward.D2", #列のクラスタリングメソッド
	column_names_gp = grid::gpar(fontsize = 18), #列ラベルの間隔
    column_names_rot = 90, #角度
	row_names_gp = grid::gpar(fontsize = 18), #行ラベルの間隔
	row_dend_width = unit(35, "mm"), #行デンドログラムの長さ
    row_dend_gp = gpar(lwd = 1, col = "black"), #デンドログラムのライン幅
    row_names_rot = 0,
	column_dend_height = unit(35,"mm"), #列デンドログラムの長さ
    column_dend_gp = gpar(lwd = 1, col = "black"), #デンドログラムのライン幅
	column_title = titlename, #列のタイトル名
	column_title_gp=grid::gpar(fontsize=18), #列のタイトルのフォントサイズ

   row_names_max_width = unit(20, "cm"),#行のマージン
    name = "Z-score" #legend title
	)
    return(P)
}

#' removing technical control 
#' 遺伝子名/サンプル名が重複するデータを1つにまとめる
#' @param data data frame which contains Samples,Targets and RQ.
#' @export
ezShrink<-function(data)
{
    data<-data %>% mutate(Ct=CtMean,
                            dCt=dCtMean,
                            ddCt=ddCtMean,
                            RQ=RQMEAN)
    data<-data[!duplicated(data[,c("Samples","Targets")]),]
    print("The data of technical replicates were removed.")
    return(data)
}


#' create significant plot
#' 有意差を表示する
#' @param data data frame which contains Samples,Targets and RQ.
#' @param pairs exp) argument list(c("A","B"),c("A","C")) creates line among A-B, A-C
#' @param textsize text size of significant annotations.
#' @param size significant bar size,
#' @param tip_length tip of significant bar size,
#' @param significant_column you can set optional data column as annotation, default="signif"
#' @return ggsignif::geom_signif
#' @export
ezSignif<-function(data,pairs=list(),textsize=3.88,size=0.5,tip_length=0,significant_column="signif",txtNotSignif="n.s.",signifSize=25,LOG10=FALSE)
{
    if(!(significant_column %in% colnames(data))){
        print("significant data is not present in the dataframe...")
        return(NULL)
    }
    Samples <- levels(data$Samples)
    #sampleが存在するかどうかのチェック
    #flags <- lapply(pairs,function(x){
    #                return( (x[1] %in% Samples) & (x[2] %in% Samples))
    #                        }) %>% unlist()
    #if(!all(flags))
    #{
    #    print(sprintf("signif annotation error::argument signif list is invalid samples=%s",Samples))
    #    return(NULL)
    #}

    xmins <- lapply(pairs,function(x){
                        return(x[1])
                            }) %>% unlist()
    xmax <- lapply(pairs,function(x){
                        return(x[2])
                            }) %>% unlist()

    #yの下限、上限を取得
    if(log10==TRUE)
    {
        RQ=log10(data$RQ)
        print("y axis is set to log10 sale") 
    }else{
        RQ=data$RQ
    }
    increment <- max(RQ)*(0.02*(textsize/5)) #y軸の拡張値
    ymin <- max(RQ) + increment #yの最大値
    ymax <- ymin + increment * NROW(xmax) + increment

    #yの位置を決定
    yextends<-seq(ymin,ymax,by=increment)
    if(NROW(yextends)>NROW(xmax)){yextends<-yextends[1:NROW(xmax)]}
    if(NROW(yextends)<NROW(xmax)){yextends<-append(yextends,yextends[NROW(yextends)]+increment)}

    #表示データを取得
    shrink <- data[!duplicated(data$Samples),]
    sig <- lapply(xmax,function(x){
                    return(unique(shrink[shrink$Samples==x,significant_column]))
                            }) %>% unlist()
    if(NROW(sig)!=NROW(xmax))
    {
        print(sprintf("signif value is invalid. sample=%s:significant=%s",NROW(xmax),NROW(sig)))
        return(NULL)
    }
    Samples<-levels(shrink)

    #データフレームを作成
    sigdata <- data.frame(Samples="",y=yextends,xmin=factor(xmins,levels=levels(shrink$Samples)),xmax=factor(xmax,levels=levels(shrink$Samples)),signifc=sig)
    sigdata[sigdata$signifc=="" | is.na(sigdata$signifc),"signifc"] <- txtNotSignif

    returner<-geom_signif(
        textsize=signifSize,
        xmin=sigdata$xmin,
        xmax=sigdata$xmax,
        y_position=sigdata$y,
        annotations=sigdata$signifc,
        size=size,
        tip_length=tip_length,
        hjust=0.5,
        vjust=0.10
    )
    return(returner)
}
