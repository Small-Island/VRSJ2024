##### 分散分析ｓＡ(１要因参加者内)

# パッケージ psych が必要

## js-STARの入力: sA
n    = c(6,6,6)

levA = 3

data = c(
 24, 16, 53,
 13, 82, 90,
 11, 80, 98,
 0, 59, 76,
 74, 66, 80,
 7, 79, 63
)

taju="holm"  # BH, holm, b

## js-STARの入力おわり

## スタック
s = gl(mean(n), levA)
A = gl(levA, 1, n[1]*levA)
sA= data.frame(s, A, data)
N = mean(n)

## 基本統計量
Title="平均と標準偏差(不偏分散の平方根)"
library( psych )
tx0=c(); txD=c()
for( i in 1:levA )
{
 kx= describe(subset(sA,A==i)$data )
 tx0= c(tx0, kx$n,kx$mean,kx$sd,kx$min,kx$max )
 txD= rbind(txD, kx[1,c(5,10:13)] )
}
tx0= matrix(round(tx0,4), nc=5,by=1 )
rownames(tx0)= paste("A.",1:levA, sep="" )
colnames(tx0)= c("   n","    Mean","     SD","   Min","   Max" )
rownames(txD)= rownames(tx0)
colnames(txD)= 
colnames(txD)= c("Median","Range","歪度","尖度","SE" )
hk=tx0[ ,2]; sd=tx0[ ,3]

## 作図
suij= 1:levA # 水準番号
bai=2 # ｙ軸の拡幅倍率

for( L in 1:2 )
{
if( L==1 ) Sd=sd else { Sd=txD[,5]; Title="平均とSE(標準誤差)" }
windows() # ライン図
xjik= c(1-.3, levA+.3) # ｘ軸の拡幅
yjik= c(min(hk)-max(sd)*bai,max(hk)+max(sd)*bai)
plot(suij,hk, bty="l",xli=xjik,yli=yjik,
  tck=0.03,  # 目盛り突起
  las=1,     # 目盛り正立
  xlab="水準番号", ylab="平　均",
  xaxt="n",
  main=Title, cex.main=1.5,cex.lab=1.2 )

meisho=paste("A",1:levA,sep="")
axis(side=1, at=1:levA, tck=0.03,labels=meisho )
arrows(suij,hk, suij,hk+Sd, lwd=L,ang=90,len=0.15 )
arrows(suij,hk, suij,hk-Sd, lwd=L,ang=90,len=0.15 )
lines (suij,hk, lwd=2,col=1)
points(suij,hk, pch=21,cex=3,bg=1 )

windows() # バー図
mitsdo= c(0:(levA-1))*10*log(1:levA) # バーの斜線密度
xjik= c(0,   levA+1) # ｘ軸の拡幅：levA+1
shoh= min(min(hk),0) # 最小平均か,ゼロ
dais= max(sd)        # 最大ＳＤ
yjik=c(shoh+dais*bai*sign(shoh),max(hk)+dais*bai) # ｙ軸の範囲

antx= barplot(hk, xli=xjik, yli=yjik,
 las=1, tck=0.02,     # ｙ軸の目盛
 col=1, den=mitsdo,   # バー斜線密度
 names.arg=paste("A",1:levA,sep=""),
 xlab="要因Ａ", # ｘ軸ラベル名
 ylab="平　均", # ｙ軸ラベル名
 main=Title, cex.main=1.5, cex.lab=1.2 )
arrows(antx,hk, antx,hk+Sd*sign(hk),lwd=L,ang=90,len=0.15 )
if(L==2) arrows(antx,hk,antx,hk-Sd*sign(hk),lwd=L,ang=90,len=0.15 )
}

# ANOVA
NAM="要因Ａ"
kx= anova(lm(data~s+A,sA) )
SS=kx$S; DF=kx$D; MS=kx$M; Fch=kx$F; Pch=kx$P
eta= SS[2]/sum(SS[2:3])
tx1= c()
tx1= matrix(round(c(SS,DF,MS,Fch,Pch,NA,eta,NA),4), nc=6 )
rownames(tx1)= c("  ｓ",NAM,"ｓ×Ａ" )
colnames(tx1)= c("SS","df","MS","Ｆ","ｐ","   ηp2" )
Fch=Fch[2]; Pch=Pch[2]

tx3= c() # 多重比較
if( Pch<0.10 & levA>2 )
{
 kx= pairwise.t.test(sA$data,sA$A,p.ad=taju, pair=1 )
 tx3= round(kx$p.v, 4)
 dimnames(tx3)= list(paste("A.",2:levA,sep=""),paste("A.",1:(levA-1),sep="") )
}

# hkR
dmat= matrix(data,nc=levA, by=1 )
rmat= cor(dmat)
keis= c()
for(x in 1:(levA-1)){ for(y in (x+1):levA) keis=c(keis,rmat[x,y] ) }
hkR= tanh( sum(nrow(dmat)*atanh(keis)) /( nrow(dmat)*length(keis) ) ) # Fisher.Z
if( sum(rmat<0)==0 | (sum(rmat<0)>0 & sum(rmat>0)==levA) ) mixR=0 else mixR=1
tmp= paste("A",1:levA, sep="" )
dimnames(rmat)= list(tmp,tmp)

# 球面性
lmod= lm(dmat~1)
kx= mauchly.test(lmod, X=~1 )
MsW=kx$stat; MsP=kx$p.va
kx= anova(lmod,X=~1,test="Spherical" )
eps= as.numeric(substr(attributes(kx)$head[5:6],29,34) )
tx9= c()
tx9= matrix(round(c(MsW,MsP,eps), 4), nc=4 )
rownames(tx9)="要因Ａ"
colnames(tx9)= c("Mauchly's W","   ｐ値","  G-G_ε"," H-F_ε" )

crP= c()
eps[ eps>1 ]=1
GGp= pf(Fch,DF[2]*eps[1],DF[3]*eps[1],low=0 )
HFp= pf(Fch,DF[2]*eps[2],DF[3]*eps[2],low=0 )
crP= matrix(round(c(DF[2],Fch,Pch,GGp,HFp),4),nc=5 )
rownames(crP)= NAM
colnames(crP)= c("df","Ｆ","ｐ"," G-G_ｐ"," H-F_ｐ")

# パワ
ESf= sqrt(eta/(1-eta) )
NCP= ESf^2*N*levA
PWz=c(); PWr=c()
d1=DF[2]; d2=DF[3]
crF= qf(0.95,d1,d2 )
PWz= 1-pf(crF,d1,d2,NCP )
PWr= 1-pf(crF,d1,d2,NCP*( 1/(1-abs(hkR)) ))
tx2= c()
tx2= matrix(round(c(ESf,PWz,PWr,hkR),4), nc=4 )
rownames(tx2)= NAM
colnames(tx2)= c("効果量ｆ"," 検出力0"," 検出力r"," 水準間相関")
if( mixR>0 ) Pow=PWz else Pow=PWr


# ■結果の書き方
PV= function(x) floor(x*1000)/1000
txt=c()
if( sum(sd==0)>0 ) txt="※分散=0 の水準があります。全データが同じ値です。以下の分析は信頼性がありません。\n\n"
Tari=0 # 多重あり
fuki= paste(" (F(",d1,",",d2,")=",round(Fch,3),", p=",PV(Pch),", ηp2=",round(eta,3),", 1-β=",round(Pow,3),")。", sep="" )
if( Pch<0.05 ) tmp="であった" else
if( Pch<0.10 ) tmp="傾向であった" else tmp="でなかった"
txt= paste(txt, "　各水準の○○得点について基本統計量をTable(tx0) or Fig.■に示す。\n　分散分析の結果 (Table(tx1)参照)，",NAM,"は有意",tmp,fuki, sep="" )

goku=""
if( Pch<0.10 )
{
 if( Pow<0.70 ){ goku="ただし"; tmp="不十分であり，信頼性が低い。" } else
 if( Pow<0.76 ) tmp="やや低いが0.70以上あり不十分ではない。" else
 if( Pow<0.80 ) tmp="ほぼ十分といえる。" else tmp="十分である。"
 txt= paste(txt, goku,"検出力 (1-β) は",tmp, sep="" )
 if( mixR>0 ) tmp="水準間の相関係数に正負が混在していたため平均相関を0と仮定し" else tmp="Fisherの重み付きZ変換値による平均相関を用いて"
 txt= paste(txt, "なお検出力の値は",tmp,"算出した。", sep="" )

if( levA>2 )
{
shus= 0
if( is.nan(MsP)==1 ){ tmp="計算不能であった"; MsP=NA; shus=shus+1 } else
  { if( MsP>=0.05 ) tmp="有意でなかったことを確認した" else { tmp="有意であった"; shus=shus+1 } }
fuki= paste(" (Mauchly's W=",round(MsW,3),", p=",PV(MsP),")。", sep="" )
txt= paste(txt, "\n　参加者内誤差についてMauchlyの球面性検定を行った結果 (Table(tx9)参照)，",tmp,fuki, sep="" )
if( (MsP<0.05 | is.na(MsP)>0) & shus>0 )
  {
   crp=crP[4]; fuki=""; ika=""
   if( crp<0.10 ){ tmp="であることを確認した"; if(crp>0.05) fuki="傾向" } else { tmp="性を得られなかった"; ika="以下，参考までに分析を進める。" }
   txt= paste(txt, "このためGreenhouse-Geisserの自由度調整係数による修正検定を行った結果 (Table(crP)参照)，",NAM,"は有意",fuki,tmp," (G-G corrected p=",PV(crp),")。",ika, sep="" )
  }
} # if(levA>)

if( levA==2 )
{
 if( hk[1]>hk[2] ) tmp="大きい" else tmp="小さい"
 if( Pch>0.05 ) fuki="傾向がある" else fuki=""
 txt= paste(txt, "\n　したがって，A1の平均",round(hk[1],3),"がA2の平均",round(hk[2],3),"よりも有意に",tmp,fuki,"ことが見いだされた。", sep="" )
} else {
 Tari= 1
 okjun= c()
 okjun= order(rownames(tx0)[ rank(tx0[,2]) ], decreasing=1 )
 fuki=c(); tmp=c( rep(",",levA-1), "" ) # 区切り
 for(i in 1:levA) fuki=paste(fuki,"A",okjun[i],tmp[i], sep="" )
 txt= paste(txt, "\n　対応のあるｔ検定を用いた多重比較 (α=0.05, 両側検定) を行った結果 (Table(tx3)参照)，", sep="" )
 kazu= sum(tx3<0.10, na.rm=1 )
 minP= min(tx3, na.rm=1)
 if( kazu<1 ) txt= paste(txt, "どの水準の平均の間にも有意差は見いだされなかった（adjusted ps>",PV(minP),")。", sep="" ) else
 { #####
txt= paste(txt, "平均の大きい順 (",fuki,") に記述すると，", sep="" )
owa= 0
for(i in 1:(levA-1))
  {
   if( owa==0 )
  {
   if( i==2 ) mata="また，" else mata=""
   hk1= round(hk[okjun[i]], 3 )
   txt= paste(txt, mata,"水準A",okjun[i],"の平均",hk1,"は", sep="" )
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
 if( i+1<j )
   {
    kt=1+(i-1)*3; bs=j-i-1; aida=substr(fuki,kt+3*1,1+kt+3*bs )
    kazu=length(yuiP); if(kazu>1) goku="ps>" else goku="p="
    txt= paste(txt, aida,"の平均と有意差がなく (adjusted ",goku,PV(yuiP[kazu]),")，次の", sep="" )
   }
 hk2= round(hk[okjun[j]], 3 )
 tch= t.test(data~A, subset(sA, A==okjun[i] | A==okjun[j]), pair=1 )$stat
 txt= paste(txt,"A",okjun[j],"の平均",hk2,"よりも有意に大き",keko,"った (t(",N-1,")=",abs(round(tch,3)),", adjusted p=",PV(adP),")。", sep="" )
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
if( Tari>0 ) txt= paste(txt, "
　以上のp値の調整には",TAJ,"の方法を用いた。", sep="" )
txt= paste(txt,"\n",sep="")

if( Tari>0 & taju=="BH" ){ txt=paste(txt, "\n［引用文献］\n",
"Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 58, 289-300.\n", sep="" ) }
H=matrix(0,nr=1)
rownames(H)=""
colnames(H)=paste(" (多数回検定時のｐ値調整は",taju,"法による) ", sep="" )

txt= paste(txt, "\n\n⇒正規分布が仮定できない場合，または水準間の球面性が仮定できない場合は下の記述を参照あるいは追加してください。\n", sep="" )
kx= friedman.test(data~A|s, sA )
pc= kx$p.va
fuki= paste(" (χ2(",kx$para,")=",round(kx$stat,3),", p=",PV(pc),")。", sep="" )
if( pc<0.05 ) tmp="であった" else 
if( pc<0.10 ) tmp="傾向であった" else tmp="でなかった"
txt= paste(txt, "\n　Friedmanの順位和検定 (rank sum test) を行った結果，各水準の順位和の差は有意",tmp,fuki, sep="" )
if( pc<0.10 & levA>2 ) txt=paste(txt, "Wilcoxonの符号化順位検定を用いた多重比較の結果は，Table(オプションのWilcox多重比較)の通りである (データにタイがある場合ｐ値は近似値)。p値の調整には",TAJ,"の方法を用いた。", sep="" )
txt= paste(txt,"\n",sep="" )
fgoT=function()
{
 options(warn=-1 )
 round(pairwise.wilcox.test(sA$data,sA$A,p.ad=taju,pair=1)$p.va, 4 )
}
TesT= function(水準=1,相手=2)
{
 tx= t.test(data~A,subset(sA,A==水準 | A==相手), pair=1 )
 ttx= paste("t(",N-1,")=",abs(round(tx$stat,4)),", non-adjusted p=",round(tx$p.va,4),"(両側)\n", sep="" )
 cat(ttx)
}

options(digits=5)
options(scipen=5)
##############################
#
#   分散分析 ｓＡ-design
#   １要因  参加者内計画
#
H#############################
tx0 # 基本統計量（SD=不偏分散の平方根）
# 歪度,尖度,SEなどは『オプション』参照

tx1 # 分散分析ｓＡ
    # 効果量ηp2 は偏イータ２乗

tx2 # 効果量ｆと検出力(1-β)
# 検出力0：水準間相関=0 (正負混在の場合)
# 検出力r：水準間相関=r (標本値から算出)

tx3 # ｔ検定による多重比較の調整後ｐ値

tx9 # 球面性検定（df=1は不要）と自由度調整係数ε
crP # 球面性検定が有意のときの修正ｐ値
    # ■参加者数Ｎ≦10なら H-F_ｐ参照可

cat(txt) # 結果の書き方

# ■オプション：［↑］⇒行頭の♯を消す⇒［Enter］
# windows();boxplot(data~A,sA) # 箱ひげ図
# txD # Median,歪度,尖度,SEなど
# round(rmat,3) # 相関行列
# fgoT()        # Wilcox多重比較(調整後ｐ値)
# TesT(水準=1, 相手=2 ) # 任意の２水準のｔ検定
# write(txt,file.choose(),ap=T) # 結果のファイル保存

# _/_/_/ Powered by js-STAR _/_/_/
