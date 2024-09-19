
ggprojection_M <- function(MA,No1=1,No2=2,OX=1,OY=1)
{
  m <- MA[["scores"]][,c(No1,No2)]
  m <- scale(m,scale=F)
  m <- as.data.frame(m)
  m[,2] <- OX*m[,2]
  m[,1] <- OY*m[,1]
  N <- nrow(m)
  gg <- ggplot(data = m,aes(x=m[,1],y=m[,2],label=1:N)) +
    geom_point(pch=21,size=9) +
    geom_text(size=4) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    labs(title='',x=bquote(U[1]),y=bquote(U[2])) +
    theme_bw() +
    theme(text=element_text(size=14))
  gg
}

ggchoro <- function(MA,MS,No1,No2,mapa,OX=1,OY=1)
{
  zw_S <- eigen(MS)
  m <- MA%*%zw_S$vectors[,c(No1,No2)]
  m <- scale(m,scale=F)
  m[,2] <- OX*m[,2]
  m[,1] <- OY*m[,1]
  n <- nrow(MA)
  ct <- vector("numeric",n)
  for (i in 1:n)
  {
    if (m[i,1]>=0 & m[i,2]>=0)
    {
      ct[i] <- 1
    } else if (m[i,1]>=0 & m[i,2]<0)
    {
      ct[i] <- 2
    } else if (m[i,1]<0 & m[i,2]>=0)
    {
      ct[i] <- 3
    }
    else if (m[i,1]<0 & m[i,2]<0)
    {
      ct[i] <- 4
    }
  }
  ggplot(data = mapa) +
    geom_sf(fill=ct) +
    theme_tufte() +
    theme_void()
}


