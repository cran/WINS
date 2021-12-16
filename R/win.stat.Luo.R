win.stat.Luo<-function(y1,y2,d,g){
  n<-length(y1)
  r11<-r10<-r21<-r20<-r31<-r30<-r41<-r40<-r51<-r50<-r61<-r60<-rep(0,n)
  for (i in 1:n){
    r11[i]<-sum(y2>=y2[i]&g==1);r10[i]<-sum(y2>=y2[i]&g==0)
    r21[i]<-sum(y2>=y2[i]&y1>=y1[i]&g==1);r20[i]<-sum(y2>=y2[i]&y1>=y1[i]&g==0)
    r31[i]<-sum(y2<=y2[i]&d<=2&g==1);r30[i]<-sum(y2<=y2[i]&d<=2&g==0)
    r41[i]<-sum(y1<=y1[i]&y2<=y2[i]&d==4&g==1)
    r40[i]<-sum(y1<=y1[i]&y2<=y2[i]&d==4&g==0)
    r51[i]<-sum(y1>=y1[i]&y2<y2[i]&d>2&g==1)
    r50[i]<-sum(y1>=y1[i]&y2<y2[i]&d>2&g==0)
    r61[i]<-sum(y1<=y1[i]&y2>y2[i]&abs(d-2.5)>1&g==1)
    r60[i]<-sum(y1<=y1[i]&y2>y2[i]&abs(d-2.5)>1&g==0)
  }
  LNA<-sum(r10[d<=2&g==1]);LNB<-sum(r11[d<=2&g==0])
  LNC<-sum(r20[d==4&g==1])+sum(r50[abs(d-2.5)>1&g==1])
  LND<-sum(r21[d==4&g==0])+sum(r51[abs(d-2.5)>1&g==0])
  wr<-(LNB+LND)/(LNA+LNC)
  wd<-(LNB+LND)-(LNA+LNC)
  la<-lb<-lc<-ld<-w<-v<-rep(0,n)
  for (i in 1:n){
    la[i]<-(1-g[i])*r31[i];lb[i]<-g[i]*r30[i]
    if (d[i]<=2) {la[i]<-la[i]+g[i]*r10[i];lb[i]<-lb[i]+(1-g[i])*r11[i]}
    lc[i]<-(1-g[i])*r41[i];ld[i]<-g[i]*r40[i]
    if (d[i]==4) {lc[i]<-lc[i]+g[i]*r20[i];ld[i]<-ld[i]+(1-g[i])*r21[i]}
    if (abs(d[i]-2.5)>1) {lc[i]<-lc[i]+g[i]*r50[i];ld[i]<-ld[i]+(1-g[i])*r51[i]}
    if (d[i]>2) {lc[i]<-lc[i]+(1-g[i])*r61[i];ld[i]<-ld[i]+g[i]*r60[i]}
  }
  la<-la/n;lb<-lb/n;lc<-lc/n;ld<-ld/n
  w<-(lb+ld)-(la+lc);barw<-mean(w);w<-w-barw
  v<-((lb+ld)-wr*(la+lc))/(LNA+LNC)*(n^2)
  vr<-sum(v^2)/n
  vd<-sum(w^2)/n
  list(wratio=wr,vwratio=vr,wdiff=wd,vwdiff=vd)
}
