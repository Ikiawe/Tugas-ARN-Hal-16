# Tugas-ARN-Hal-16
Program R Analisis Regresi NonParametrik Halaman 16


#_____________________________AWAL ESTIMASI KERNEL________________________

kernel1=function(u)
{
  exp(-0.5*(u^2))/sqrt(6.28)
}
kernelgaus=function(x,y)
{
  data=cbind(x,y)
  x=data[,1]
  y=data[,2]
  hl=as.numeric(readline("Bandwidth Optimal="))
  n=length(y)
  u=matrix(0,ncol = n,nrow = n)
  for (i in 1:n) {
    for (j in 1:n) {
      u[i,j]=(x[i]-x[j])/hl
    }
  }
  w=matrix(0,ncol=n,nrow = n)
  for(i in 1:n){
    for(j in 1:n){
      w[i,j]=kernel1(u[i,j])
    }
  }
  H=matrix(0,ncol = n,nrow = n)
  for(i in 1:n){
    for(j in 1:n){
      H[i,j]=w[i,j]/sum(w[i,])
    }
  }
  mhlambda=H%*%y
  cat("\n")
  cat("---------------------------------- \n")
  cat("         y*        y*topi error \n")
  cat("---------------------------------- \n")
  for(i in 1:n)
  {
    cat("\n",y[i],"",round(mhlambda[i],3),"",round(y[i]-mhlambda[i],3),sep = "\t")
  }
  win.graph()
  plot(x,y,type = "p",xlim = c(min(x),max(x)),ylim = c(min(y),max(y)),ylab = "Rings",xlab = "Length")
  par(new=T)
  plot(x,mhlambda,type = "l",xlim = c(min(x),max(x)),ylim = c(min(mhlambda),max(mhlambda)),ylab = "Rings",xlab = "Length")
}

#_____________________________AKHIR ESTIMASI KERNEL________________________



#_____________________________AWAL KERNEL GCV________________________

kernelGaussian=function(u)
{
 exp(-0.5*(u^2))/sqrt(2*pi)
}
kernelGCV=function(x,y)
{
	bb=as.numeric(readline("batas bawah bandwidth = "))
	ba=as.numeric(readline("batas atas bandwidth = "))
	inc=as.numeric(readline("increment = "))
	vh=seq(bb,ba,inc)
	nvh=length(vh)
	n=length(y)
	GCV=rep(0,nvh)
	cat("\n=========================")
	cat("\n lamda            GCV")
	cat("\n=========================")
	for(k in 1:nvh)
	{
		u=matrix(0,ncol=n,nrow=n)
		for(i in 1:n)
		{
		 for(j in 1:n)
		 { u[i,j]=(x[i]-x[j])/vh[k]
		 }
		}
		jumlah=matrix(0,ncol=1,nrow=n)
		for(i in 1:n)
		{
		 for(j in 1:n)
		 { jumlah[i,1]=jumlah[i,1]+kernelGaussian(u[i,j])
		 }
		}
		H=matrix(0,ncol=n,nrow=n)
		for(i in 1:n)
		{
		 for(j in 1:n)
		 { H[i,j]=kernelGaussian(u[i,j])/jumlah[i,1]
		 }
		}
		mhlamda=H%*%y
		atas=t(y-mhlamda)%*%(y-mhlamda)/n
		bawah=(1-(sum(diag(H)))/n)^2
		GCV[k]=atas/bawah
	}
s=matrix(c(vh,GCV),nvh,2)
print(s)
GCVminimum=min(GCV)
h.optimal=s[s[,2]==GCVminimum,1]
cat("BERDASARKAN TABEL DIATAS DIKETAHUI:\n")
cat("bandwidth optimal = ",h.optimal,"\n","Nilai GCV minimum = ",GCVminimum)
cat("\n")
plot(vh,GCV,type="l",lwd=2)
}

#_____________________________AKHIR KERNEL GCV________________________
