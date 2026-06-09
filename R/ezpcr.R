
library(tidyverse)
library(dplyr)
#グラフ用パッケージ
library(ggplot2)

#dplyrを使った方法
#現状一番シンプル
library(tidyverse)
library(readxl)



 #' ezread
 #'
 #' This function reads data exported by QuantStudio.
 #'
 #' @param dir folder/directory.
 #' @param skip skip rows.
 #' @param samplename column of sample names.
 #' @param targetname column of target names.
 #' @param Ct column of Ct values.
 #' @return dataframe.
 #' @export
 #' @examples

#' @export
ezread<-function(dir="./",skip=40,samplename='Sample Name',targetname='Target Name',CT="CT")
{
    files<-list.files(dir,pattern="*.csv|*.txt|*.xlsx|*.xls")
    tables<-lapply(files,function(f)
    {
        if(endsWith(".txt",f)||endsWith(".csv",f))
        {
            A0<-read.csv(f,header=T,sep=",",skip=skip)
        }
        if(endsWith(".xlsx",f)){
            A0<-read_xlsx(f,skip=40)
        }
        if(endsWith(".xls",f)){
            A0<-read_xls(f,skip=40)
        }

    A0<-subset(A0,subset=!is.na(Omit))
    A0<-data.frame(Samples=A0[,samplename],Targets=A0[,targetname],Ct=as.numeric(A0[,"CT"]))
    return(A0)
    })
    atable<-bind_rows(tables)
    return(atable)
}#end function ezread


 #' ezcalc
 #'
 #' This function calculates dCt ddCt RQ.
 #'
 #' @param df dataframe obtained from ezread.
 #' @param internalControl a sample name of biological-control.
 #' @param internalControl a gene of internal-control.
 #' @param CtMax if Ct Value is undetermined(NA), this value sets to CtMax.
 #' @return dataframe.
 #' @examples

#' @export
ezcalc<-function(df,biologicalControl,internalControl="GAPDH",CtMax=40)
{
    if(length( subset(df,subset=Targets==biologicalControl))==0)
    {
        print(sprintf(" %s is not containd in this data!"))
    }

#Ctを数値化
    df<-df %>% dplyr::mutate(Ct=as.numeric(Ct)) %>% dplyr::mutate(Ct=if_else(is.na(Ct) | Ct>38,CtMax,Ct))

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
    bcons<-paste(grep(biologicalControl,x3$Samples,value=TRUE),collapse=",")
    print(sprintf("%s was used as biological control",bcons))
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
    return(x3)

}#end function ezcalc



 #' ezgraph
 #'
 #' This function creates plots.
 #'
 #' @param df dataframe obtained from ezread.
 #' @param internalControl a sample name of biological-control.
 #' @param internalControl a gene of internal-control.
 #' @param CtMax if Ct Value is undetermined(NA), this value sets to CtMax.
 #' @return dataframe.
 #' @examples

 #' @export
ezGraph<-function(data,samplenames=NA,dot=FALSE,linewidth=2,textSize=22,titleSize=32,legendPosition="none",genes=NA)
{
    if(is.na(genes)){
        genes<-unique(data$Targets)
    }
    #samplenamesに含まれるサンプルのデータのみを抽出
    if(!is.na(samplenames))
    {
        samplenames<-unique(data$Samples)
    }
        datas<-lapply(samplenames,function(x)
        {
            d<-data[grepl(x,data$Samples),]
        })
        data0<-bind_rows(datas)
        data0$Samples<-levels(data0$Samples,levels=samplenames)

    PS<-lapply(genes,function(g)
        {
            data1<-subset(data0,subset=Targets==g)
	p <- ggplot(data = data1 , aes(x = Samples, y = RQ, fill= Samples)) + #aesで使用するパラメータを指定

    #plot ===
	stat_summary(geom="bar", fun=mean, color="black",size=linewidth,width = 0.8 )+ #自動的に平均化した棒グラフを作ってくれる stat_summary
    #group化する際には barとerrorbarにpositoin = position_dodge(0.5)などを付けること
	stat_summary(geom="errorbar",fun.data=mean_sdl,fun.args = list(mult = 1), width=0.5 ,size=linewidth)+ #エラーバーも自動で計算してくれる 標準誤差(mean_se)を使用 標準偏差は(mean_sdl,fun.args = list(mult=1)) あるいはggpubrのmean_sd
    geom_jitter(color = "black", fill = "red", size = 2.5,shape = 21 , width = 0.2) + #position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.2)# groupingのときにつける width,fillは除外すること
 
    #add text ===
    stat_summary( fun = "max", geom = "text",
       aes( label = TEXT_ ), vjust = -0.5, 
       size = 5,  # テキストのサイズを設定
       color = "black",  # テキストの色を青に設定
       fontface = "bold"  # テキストを太字に設定
     )+    
  geom_signif(
    comparisons = list(c("SampleA", "SampleB")),  # 比較するペア
    test = "wilcox.test",                  # t検定を使用
    map_signif_level = TRUE           # p<0.05,* / p<0.01,** / p<0.001,*** に自動変換
  )+
    #settings ===
	scale_y_continuous(limit=c(0,NA) , expand=c(0,0))+ #x軸の最小値を０に固定 NAをmax(data)*xにすると、最大値を拡張できる
	scale_fill_manual(values=c("black","black"))+ #色指定
	ggtitle(GENENAME)+ 

	theme_classic()+ #シンプルなデザインに変更

	theme(
        plot.title = element_text(hjust=0.5), #theme : 軸の太さなどの細かい点を指定
	axis.title.x = element_blank(),
    axis.title.y = element_text( size = textSize, vjust = 2),
	axis.text.x = element_text(size=textSize, color="black" ),
    axis.text.y = element_text(size=textSize, color="black" ),

	title=element_text(size=titleSize),

	axis.line = element_line(linewidth =linewidth),
	axis.ticks = element_line(size=linewidth),
    axis.ticks.x = element_blank(),
	axis.ticks.length.y = unit(2,"mm"),
    legend.position = legendPosition #without legend
    )
    })#end lapply genes
    names(PS)<-genes
    return(PS)
} #end function ezgraph


ezPng<-function(plots,plotname="plot",width=5,height=5)
{

    if(!("svglite" %in% installed.packages()))
    {
        install.packages("svglite")
    }
    library("svglite")

    for(w in names(P))
    {
        ggsave(sprintf("./%s-%s.png",plotname,w),plot=P[[w]],device=png,width=width,height=height,dpi=350,bg="white")
    }
}


ezSvg<-function(plots,plotname="plot",width=5,height=5)
{
    if(!("svglite" %in% installed.packages()))
    {
        install.packages("svglite")
    }
    library("svglite")

    for(w in names(P))
    {
        ggsave(sprintf("./%s-%s.svg",plotname,w),plot=P[[w]],device=svglite,width=width,height=height,dpi=350,bg="white")
    }

}

