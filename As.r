##### 分散分析Ａｓ(１要因参加者間)

# パッケージ psych が必要

## js-STARの入力: As
levA = 3
n    = c(5,5,5)

data = c(
 24, 13, 11, 0, 7,
 16, 82, 80, 59, 79,
 53, 90, 98, 76, 63
 )

taju="BH"  # BH, holm, b

## js-STARの入力おわり

## スタック
A= factor(rep(1:levA, n) )
N= sum(n)
s= gl(N, 1)
As= data.frame(A, s, data)

## 基本統計量
Title="平均と標準偏差(不偏分散の平方根)"
library( psych )
tx0=c(); txD=c()
for( i in 1:levA )
{
 kx= describe(subset(As,A==i)$data )
 tx0= c(tx0, kx$n,kx$mean,kx$sd,kx$min,kx$max )
 txD= rbind(txD, kx[1,c(5,10:13)] )
}
tx0= matrix(round(tx0,4), nc=5,by=1 )
rownames(tx0)= paste("A.",1:levA, sep="" )
colnames(tx0)= c("   n","    Mean","     SD","   Min","   Max" )
rownames(txD)= rownames(tx0)
colnames(txD)= c("Median","Range","歪度","尖度","SE" )
hk=tx0[ ,2]; sd=tx0[ ,3]

## 作図
mitsdo= c(0:(levA-1))*10*log(1:levA) # バーの斜線密度
xjik= c(0,   levA+1) # ｘ軸の拡幅：levA+1
bai = 2              # ｙ軸の拡幅倍率
shoh= min(min(hk),0) # 最小平均か,ゼロ
dais= max(sd)        # 最大ＳＤ
yjik=c(shoh+dais*bai*sign(shoh),max(hk)+dais*bai) # ｙ軸の範囲

for( L in 1:2 )
{
if( L==1 ) Sd=sd else { Sd=txD[,5]; Title="平均とSE(標準誤差)" }
windows()
antx= barplot(hk, xli=xjik, yli=yjik,
  las=1, tck=0.02,   # ｙ軸の目盛
  col=1, den=mitsdo, # バー斜線密度
  names.arg=paste("A",1:levA,sep=""),
  xlab="要因Ａ", # ｘ軸ラベル名
  ylab="平　均", # ｙ軸ラベル名
  main=Title, cex.main=1.5,cex.lab=1.2,cex.names=1.2 )
arrows(antx,hk, antx,hk+Sd*sign(hk), lwd=L,ang=90,len=0.15 )
if(L==2) arrows(antx,hk,antx,hk-Sd*sign(hk),lwd=L,ang=90,len=0.15 )
}

## 分散分析
NAM="要因Ａ"
kx= anova(lm(data~A,As) )
SS=kx$S; DF=kx$D; MS=SS/DF; Fch=kx$F; Pch=kx$P
MS=SS/DF; eta=SS[1]/sum(SS)
tx1= c()
tx1= matrix(round(c(SS,DF,MS,Fch,Pch,eta,NA), 4), nc=6 )
colnames(tx1)= c("SS","df","MS","Ｆ","ｐ","η2" )
rownames(tx1)= c(NAM,"  ｓ")
d1=DF[1];d2=DF[2]; Fch=Fch[1];Pch=Pch[1]

tx3= c() # 多重比較
if(Pch<0.10 & levA>2)
{
 tx3=round(pairwise.t.test(As$data,As$A,p.ad=taju)$p.va, 4 )
 dimnames(tx3)= list(paste("A.",2:levA,sep=""), paste("A.",1:(levA-1),sep="") )
}

# 分散の均一性
tx8= c()
kx= bartlett.test(data~A, data=As )
Bst=kx$stat; Bpv=kx$p.va
tx8= matrix(round(c(Bst,DF[1],Bpv), 4), nc=3 )
rownames(tx8)="Bartlett Test"
colnames(tx8)= c("    χ2","  df","     ｐ" )

## パワ
ESf= sqrt(Fch*d1/d2 )
Pow= 1-pf(qf(0.95,d1,d2),d1,d2, ncp=ESf^2*N )

nn= NA
if( Pow<0.80 )
{
for(i in (N-1):10000 )
{
 D2= i-levA
 beta= pf(qf(0.95,d1,D2),d1,D2, ESf^2*i )
 if( beta<0.20 ){ nn=c(i, nn ); break }
}
} else {
for(i in (N+1):2 )
{
 D2= i-levA
 beta= pf(qf(0.95,d1,D2),d1,D2, ESf^2*i )
 if( beta>0.20 ) break
 nn= c(i, nn )
}
} # if(Pow)
tx4= c()
tx4= matrix(round(c(ESf,Pow,N,nn[1]), 4), nc=4 )
rownames(tx4)= NAM
colnames(tx4)= c("効果量ｆ"," 検出力","今回のＮ","次回のＮ" )

# 修正Ｆ検定
kx= oneway.test(data~A, data=As )
fch= kx$stat
tx5= matrix(round(c(fch,kx$para,kx$p.val), 4), nc=4 )
Df2= kx$para[2]
tx5= rbind(tx5, matrix(round(c(Fch,DF[1],Df2,pf(Fch,d1,Df2,low=F)), 4), nc=4) )
rownames(tx5)= c("再計算のＦ比検定","元のＦ比検定")
colnames(tx5)= c("Ｆ"," df1"," df2","     ｐ" )


# ■結果の書き方
PV= function(x) floor(x*1000)/1000
txt=c()
if( sum(sd==0)>0 ) txt="※分散=0 の群があります。全データが同じ値です。以下の分析は信頼性がありません。\n\n"
Tari=0 # 多重あり
fuki= paste(" (F(",d1,",",d2,")=",round(Fch,3),", p=",PV(Pch),", η2=",round(eta,3),", 1-β=",round(Pow,3),")。", sep="" )
if( Pch<0.05 ) tmp="であった" else
if( Pch<0.10 ) tmp="傾向であった" else tmp="でなかった"
txt= paste(txt, "　各群の○○得点について基本統計量をTable(tx0) or Fig.■に示す。\n　分散分析の結果 (Table(tx1)参照)，",NAM,"は有意",tmp,fuki, sep="" )

goku=""; fuki=""; tmp=""
if( Pch<0.10 )
{
 if( Pow<0.70 ){ goku="ただし"; tmp="不十分であり，信頼性が低い。" } else
 if( Pow<0.76 ) tmp="やや低いが0.70以上あり不十分ではない。" else
 if( Pow<0.80 ) tmp="ほぼ十分といえる。" else tmp="十分である。"
 txt= paste(txt, goku,"検出力 (1-β) は",tmp, sep="" )
}

# 分散均一性
ika= 0
if( Pch<0.10 )
{
 fuki= paste(" (χ2(",d1,")=",round(Bst,3),", p=",PV(Bpv),")。", sep="" )
 if( Bpv>0.05 ) tmp="ないことを確認した" else tmp="あった"
 txt= paste(txt, "\n　参加者間の分散の均一性についてBartlett検定を行った結果 (Table(tx8)参照)，有意で",tmp,fuki, sep="" )
 if( Bpv<0.05 )
{
 tmp= tx5[1,1:4]
 pch= PV(tmp[4])
 fuki= paste(" (F(",d1,",",round(tmp[3],3),")=",round(tmp[1],3),", p=",pch,")。", sep="" )
 goku=""; ika=0
 if( pch<0.10 ){ tmp="であることを確認した"; if(pch>0.05) goku="傾向" } else
               { tmp="性を得ることができなかった"; ika=1 }
 txt= paste(txt, "Welchの方法による修正Ｆ検定 (Ｒのoneway.test関数) を行った結果，要因Ａは有意",goku,tmp,fuki, sep="" )
}
} # if(Pch)
if( ika==1 ) txt=paste(txt, "以下，参考までに分析を進める。\n", sep="" )

# 多重比較
if( Pch<0.10 )
{
if( levA==2 )
{
 if( ika==1 ) goku="参考知見" else goku="結果"
 if( hk[1]>hk[2] ) tmp="大きい" else tmp="小さい"
 if( Pch>0.05 ) fuki="傾向がある" else fuki=""
 txt= paste(txt, "\n　したがって，A1の平均",round(hk[1],3),"がA2の平均",round(hk[2],3),"よりも有意に",tmp,fuki,"ことが見いだされた。", sep="" )
} else {
 Tari= 1
 okjun= c()
 okjun= order(rownames(tx0)[ rank(tx0[,2]) ], decreasing=1 )
 fuki=c(); tmp=c( rep(",",levA-1), "" ) # 区切り
 for(i in 1:levA) fuki=paste(fuki,"A",okjun[i],tmp[i], sep="" )
 txt= paste(txt, "\n　プールドSDを用いたｔ検定による多重比較 (α=0.05, 両側検定) を行った結果 (Table(tx3)参照)，", sep="" )
 kazu= sum(tx3<0.10, na.rm=1 )
 minP= min(tx3, na.rm=1)
 if( kazu<1 ) txt= paste(txt, "どの群の平均の間にも有意差は見いだされなかった（adjusted ps>",PV(minP),")。", sep="" ) else
{ #####
 txt= paste(txt, "平均の大きい順 (",fuki,") に記述すると，", sep="" )
 owa= 0
 for(i in 1:(levA-1))
   {
    if( owa==0 )
   {
    if( i==2 ) mata="また，" else mata=""
    hk1=round(hk[okjun[i]], 3); n1=n[okjun[i]]
    txt= paste(txt, mata,"A",okjun[i],"の平均",hk1,"は", sep="" )
    yuiP=c(); sumi=0
 for(j in (i+1):levA)
   {
    if( sumi==0 )
   {
    gyo=okjun[j]; rez=okjun[i]
    if( gyo>rez ) adP=tx3[gyo-1,rez] else adP=tx3[rez-1,gyo]
    if( adP<0.10 )
{
keko="い傾向があ"; if( adP<0.05 ) keko="か"
hk2=round(hk[okjun[j]], 3); n2=n[okjun[j]]
tch= round( abs(hk1-hk2)/sqrt( (1/n1+1/n2)*MS[2] ), 3)
if( i+1<j )
  {
   kt=1+(i-1)*3; bs=j-i-1; aida=substr(fuki,kt+3*1,1+kt+3*bs )
   kazu=length(yuiP)
   if(kazu>1) goku="ps>" else goku="p="
   txt= paste(txt, aida,"の平均と有意差がなく (adjusted ",goku,PV(yuiP[kazu]),")，次の", sep="" )
  }
txt= paste(txt,"A",okjun[j],"の平均",hk2,"よりも有意に大き",keko,"った (t(",d2,")=",tch,", adjusted p=",PV(adP),")。", sep="" )
sumi= 1
} else {
 yuiP= c(yuiP, adP )
 kazu=length(yuiP); if(kazu>1) goku="ps>" else goku="p="
 if( j==levA ){ txt=paste(txt, "以降の平均と有意差がなかった (adjusted ",goku,PV(yuiP[kazu]),")。", sep="" ); owa=1 }
}
  }  }  }  } # for(i)(j)
} #####
} # if(levA)
} # if(Pch)

if( taju=="BH" ) TAJ="Benjamini & Hochberg (1995) "
if( taju=="holm" ) TAJ="Holm"
if( taju=="b" ) TAJ="Bonferroni"
if( Tari>0 ) txt= paste(txt, "\n　以上のp値の調整には",TAJ,"の方法を用いた。", sep="" )
txt= paste(txt,"\n",sep="")
if( Tari>0 & taju=="BH" ){ txt=paste(txt, "\n［引用文献］\n",
"Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 58, 289-300.\n", sep="" ) }
H=matrix(0,nr=1)
rownames(H)=""
colnames(H)=paste(" (多数回検定時のｐ値調整は",taju,"法による) ", sep="" )

txt= paste(txt, "\n\n⇒データ分布に正規性を仮定できない場合は下の記述を参照してください。\n", sep="" )
kx= kruskal.test(data~A, As )
pc= kx$p.va
fuki= paste(" (χ2(",kx$para,")=",round(kx$stat,3),", p=",PV(pc),")。", sep="" )
if( pc<0.05 ) tmp="であった" else 
if( pc<0.10 ) tmp="傾向であった" else tmp="でなかった"
txt= paste(txt, "\n　Kruskal-Wallisの順位和検定 (rank sum test) を行った結果，各群の順位和の差は有意",tmp,fuki, sep="" )
if( pc<0.10 & levA>2 ) txt=paste(txt, "Wilcoxonの順位和検定を用いた多重比較の結果は，Table(オプション『Wilcox多重比較(調整後ｐ値)』)の通りである (データにタイがある場合ｐ値は近似値)。p値の調整には",TAJ,"の方法を用いた。", sep="" )
txt= paste(txt,"\n",sep="" )
wlx=function()
{
 options(warn=-1 )
 round(pairwise.wilcox.test(As$data,As$A,p.ad=taju)$p.va, 4 )
}
TesT=function(群=1, 相手=2 )
{
 tch= abs(hk[群]-hk[相手])/sqrt((1/n[群]+1/n[相手])*MS[2])
 nPV= pt(tch,d2,low=0)*2
 cat(paste("t(",d2,")=",round(tch,4),", non-adjusted p=",round(nPV,4),"(両側)\n", sep="") )
}

if( levA>2 )
{
kx9= cor.test(as.numeric(As$A), As$data )
if( kx9$p.va<0.10 )
txt= paste(txt, "\n⇒水準番号（群）とデータに有意な相関が見られます (r=",round(kx9$est,3),", t(",kx9$para,")=",round(kx9$stat,3),", p=",PV(kx9$p.va),")。水準が年齢段階や所得段階のような連続変量ならば，相関係数の計算または回帰分析を実行してみることも有望です。\n", sep ="" )
}

options(digits=5)
options(scipen=5)
##############################
#
#   分散分析 Ａｓ-design
#   １要因  参加者間計画
#
H#############################
tx0 # 基本統計量（SD=不偏分散の平方根）
# 歪度,尖度,SEなどは『オプション』参照

tx1 # 分散分析Ａｓ
    # 効果量η2 はイータ２乗

tx4 # 効果量ｆ,検出力(1-β),追試用Ｎ
    # 効果量ｆの評価：大=0.40,中=0.25,小=0.10
    # 「次回のＮ」はα=0.05,検出力=0.80を想定

tx3 # ｔ検定による多重比較の調整後ｐ値

tx8 # 分散の均一性検定

tx5 # Welchの方法を用いた修正Ｆ検定

cat(txt) # 結果の書き方

# ■オプション：［↑］⇒行頭の♯を消す⇒［Enter］
# windows();boxplot(data~A,d=As) # 箱ひげ図
# txD   # Median,歪度,尖度,SEなど
# wlx() # Wilcox多重比較(調整後ｐ値)
# TesT(群=1, 相手=2 ) # 任意の２群のｔ検定
# write(txt,file.choose(),ap=T) # 結果のファイル保存

# _/_/_/ Powered by js-STAR _/_/_/
