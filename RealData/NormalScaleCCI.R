

xx<-scan(what = "character",sep="\n")
1  &  -.118 & .364  \\  
2  &  -.079 & .453  \\  
3  &  .000 & .588  \\  
4  &  .085 & .679  \\  
5  &  .037 & .295  \\  
6  &  -.027 & .614  \\  
7  &  -.029 & .683  \\  
8  &  -.078 & .709  \\  
9  &  -.057 & .102  \\  
10  &  5.437 & .696  \\  
11  &  6.580 & .584  \\  
12  &  .213 & 2.143  \\  
13  &  .463 & 3.473  \\  
14  &  .002 & .101  \\  
15  &  -.064 & .340  \\  
16  &  .035 & .200  \\  
17  &  .010 & .231  \\  
18  &  -.118 & .787  \\  
19  &  -.113 & .933  \\  
20  &  -.137 & .274  \\  
21  &  -.045 & .380  \\  
22  &  .059 & .098  \\ 

df<-data.frame(item=as.numeric(sapply(xx,function(x) strsplit(x,"&")[[1]][1])),
           A1 = as.numeric(sapply(xx,function(x) strsplit(x,"&")[[1]][2])),
           A2 = as.numeric(sapply(xx,function(x) substr(strsplit(x,"&")[[1]][3],1,5))))
df[,2]<-df[,2]/sqrt(1+df[,2]^2)         
df[,3]<-df[,3]/sqrt(1+df[,3]^2)         
for (i in 1:nrow(df)) cat(df[i,1],"&",apaformat(df[i,2],digits = 3),"&",apaformat(df[i,3],digits = 3),"\\\\ \n")
