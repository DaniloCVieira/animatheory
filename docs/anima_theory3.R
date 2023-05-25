library(colorRamps)
library(magick)
library(raster)
library(animation)
#' @export
plotpool<-function(pool_table,immigs=NA)
{
  pool.positions=pool_table[,1:2]
  pool=pool_table[,3]
  cols.pool=pool_table[,4]


  plot(pool.positions, axes=F, ann=F,pch=pool, col=cols.pool, cex=runif(8, max=2, min=1.5), xlim=c(0,6), ylim=c(0,6))
  box(lty=2,col='gray50')
  par(xpd=T)



}

#' @export
plotrichness<-function(species_table,pool_table,time=0, n_sp=length(table(na.omit(species_table[,3]))), title2='Community',cycles=100,main)
{
  par(mar=c(5.1, 4.1, 9, 2.1))
  time<-attr(species_table,"S")
  #layout(matrix(c(1,2,3),ncol=3), widths = c(2,2,2))
  { plot(seq(0,cycles, length.out=cycles),seq(0,max(n_sp), length.out=cycles), type="n", xlab="Time steps", ylab="Number of species", bty="l", ylim=c(0,nrow(species_table)), las=1)
    title("Species Richness",font.main=1)
    if(length(time)==1){

      points(0,n_sp,pch=16,col="red",cex=.7)
      text(0,y=n_sp,labels=n_sp,col="red",adj=0)
    } else {
      text(length(time),y=n_sp,labels=n_sp,col="red")
      lines(1:length(time),time,pch=16,col="red")
    }



  }
}

#' @export
plotrank<-function(species_table,pool_table,time=0, n_sp=length(table(na.omit(species_table[,3]))), title2='Community',cycles=100,main)
{
  {

    species_table[c(2,8,15),3]<-NA
    time<-attr(species_table,"S")
    spes<-na.omit(species_table[,3])
    colspes<- species_table[-which(is.na(species_table[,3])),4]
    colpool<-pool_table[,4]
    names(colpool)<-pool_table[,3]
    sprank<-table(spes)[order(table(spes), decreasing = T)]
    sprank<-sprank/sum(sprank)*100
    colrank<-colpool[names(sprank)]




  par(xpd=F)
  par(mar=c(5.1, 4.1, 9, 2.1))
  plot(sprank,type='b', axes=F, xlab="species abundance rank", ylab="number of species", ylim=c(0,100))
  title("Rank abundance",font.main=1)
  axis(1, labels=NA, at=1:length(names(sprank)))
  par(xpd=T)
  text(1:length(names(sprank)),rep(-8,length(sprank)),names(sprank),col=colrank)
  axis(2, las=1, at=c(0,20,40,60,80,100), labels=c(0,20,40,60,80,100))
  par(xpd=T)

  mtext(main, side = 3, line = -3, outer = TRUE,cex=1.5,font=2)
  par(xpd=F)
  }
}

#' @export
plotplay<-function(species_table,pool_table,time=0, n_sp=length(table(na.omit(species_table[,3]))), title2='Community',cycles=100,main, xdim,ydim)

{

  par(mfrow=c(1,3),mar=c(5.1, 4.1, 9, 2.1), cex=1.25)

  time<-attr(species_table,"S")

  {
    #plotcomm

    par(xpd=F, mar=c(2, 1, 6, 0), cex=1.3)

    positions=species_table[,1:2]
    species=species_table[,3]
    cols=species_table[,4]


    lims=c(min(species_table[,1])-1, max(species_table[,1])+1)
    plot( positions[,1], positions[,2], xlim=lims, ylim=lims, pch=c(species), cex=1.5, col=cols, ann=F, axes=F, type="n")
    par(xpd=T)
    points(rep(-.5,nrow(pool_table)),seq(min(species_table[,2])-2,max(species_table[,2])+4,length.out=(nrow(pool_table))), pch=pool_table[,3],col=pool_table[,4], cex=.9)
    par(xpd=F)
    if(ncol(species_table)!=4)
    {
      m<-matrix(species_table[,5],max(species_table[,1]), byrow = F)
      m<-apply(m,2,rev)
      r <- raster(xmn = min(species_table[,1])-.5, xmx = max(species_table[,1])+.5, ymn = min(species_table[,2])-.5, ymx = max(species_table[,1])+.5, nrows = xdim, ncols = xdim)
      r[] <-rev(unlist(m))

      plot(r,  col=adjustcolor(matlab.like(length(table(species_table[,5]))),0.4),  add=T,legend=F)


      par(xpd=T,cex=1)
      points(rep(-.5,nrow(pool_table)),seq(min(species_table[,2])-2,max(species_table[,2])+4,length.out=(nrow(pool_table))), pch=pool_table[,3],col=pool_table[,4], cex=.9)
      text(rep(-.1,nrow(pool_table)),seq(min(species_table[,2])-2,max(species_table[,2])+4,length.out=(nrow(pool_table))), label=round(pool_table[,5],1),col=pool_table[,4],cex=.9)



      text(species_table[,1],species_table[,2]+.3, label=paste("EV",species_table[,"EV"]), cex=.7)
      text(species_table[,1],species_table[,2]-.15, label=paste0("SO",round(species_table[,"SO"],1)), cex=.7)
      text(species_table[,1],species_table[,2]-.3, label=round(species_table[,"probdie"],2), cex=.7)
      par(xpd=F, cex=1.3)
    }




    points( positions[,1], positions[,2], pch=c(species), cex=1.5, col=cols)


  getlims<-function(lims){
    lims<-cbind( do.call(rbind,lapply(strsplit( paste(lims[1]-.5,seq(lims[1]+.5,lims[2]+.5, by=1))," "), as.numeric)),
                 do.call(rbind,lapply(strsplit( paste(lims[2]+.5,seq(lims[1]+.5,lims[2]+.5, by=1))," "), as.numeric)))
    return(lims)
  }

  xlines<-getlims(range(species_table[,1]))
  ylines<-getlims(range(species_table[,2]))

  segments(xlines[,1],xlines[,2],xlines[,3],xlines[,4], lty=2, col="darkgreen")
  segments(ylines[,2],ylines[,1],ylines[,4],ylines[,3], lty=2, col="darkgreen")
  segments(0.5,0.5,0.5,max(species_table[,2])+.5, lty=2, col="darkgreen")
  segments(0.5,0.5,max(species_table[,1])+.5,0.5, lty=2, col="darkgreen")
  }

  title(title2,font.main=1, line=1)

}

#' @export
initiate<-function(n_sp=25,xdim=5,ydim=5, n_pool=35,emoj_cat="Animals & Nature", cycles, build_envi=T, npatchs=3,sel=0.1, main,seed0)
{



  emoj<-data.frame(emo::jis)
  em<-emoj$emoji[which(emoj$group%in%emoj_cat)]
  spcol<-matlab.like(length(em))
  emojdf<-data.frame(em,spcol)
  set.seed(seed0)
  pool.rows<-order(sample(1:nrow(emojdf), n_pool))
  pool.sp<-emojdf[pool.rows,1]
  cols.pool<-matlab.like(length(pool.sp))

  pool_table<-data.frame(x=runif(length(pool.sp), min=1, max=xdim),y= runif(length(pool.sp), min=1, max=ydim))
  pool_table[,3]<-pool.sp
  pool_table[,4]<-cols.pool

  seed0<-seed0+1
  set.seed(seed0)
  sps.rows<-sample(1:nrow(pool_table), xdim*ydim, replace = T)
  sps<-pool_table[sps.rows,3]
  spcolors<-pool_table[sps.rows,4]

  species_table<-expand.grid(data.frame(x=1:xdim, y=1:ydim))
  species_table[,3]<-sps
  species_table[,4]<-spcolors


  S<-length(table(na.omit(pool_table[,3])))
  attr(species_table,"S")<-S
  attr(species_table,"seed")<-seed0
  if(isTRUE(build_envi)){
    species_table<-makeenvi(npatchs,species_table,pool_table)

    pool_table[,4]<- matlab.like(nrow(pool_table))
    #set.seed(1)
    #sps.rows<-sample(1:nrow(pool_table), xdim*ydim, replace = T)
    species_table[,3]<-pool_table[sps.rows,3]
    species_table[,4]<-pool_table[sps.rows,4]
    fitness<-seq(1,max(species_table[,5]), length.out=nrow(pool_table))
    pool_table[,"SO"]<-fitness
    names(fitness)<-pool_table[,3]
    spcom<-species_table[,3]
    names(spcom)<-species_table[,3]
    fit<-round(as.vector(fitness[names(spcom)]),1)
    species_table[,'SO']<-fit
    probdie<-get_selection(species_table,pool_table, sel=sel)
    species_table[,'probdie']<-probdie



    # text(species_table[,1],species_table[,2]+.3, label=paste(species_table[,5]), cex=.7)

  }
  par(mfrow=c(1,3),mar=c(5.1, 4.1, 9, 2.1), cex=1.25)
  plotplay(species_table,pool_table,title2="Start community", cycles=cycles, main=main, xdim=xdim,ydim=ydim)

  plotrichness(species_table,pool_table,title2="Start community", cycles=cycles, main=main)
  plotrank(species_table,pool_table,title2="Start community", cycles=cycles, main=main)
  attr(species_table,"seed")<-seed0


  return(list(sp=species_table,pool=pool_table))
}


#' @export
die<-function(species_table,pool_table,prob=0.1, cycles,main, xdim=xdim, ydim=ydim, n_pool=n_pool)
{

  seed0<-attr(species_table,"seed")+1

  if(prob>=1)
  {
    set.seed(seed0)
    sampled.rows=sample(1:nrow(species_table),prob)
    deaths.rows<-sampled.rows
  } else {
    repeat {

      set.seed(seed0)

      sampled.rows=sample(c(T,F), length(na.omit(species_table[,3])),prob=c(prob,1-prob), replace = T)
      seed0<-seed0+1

      if(sum(sampled.rows)!=0){ break}
    }
    deaths.rows<- which(sampled.rows==T)
  }





  #deaths.rows<-sample(1:nrow(species_table),n_ext)
  species_table[deaths.rows,3:4]<-NA
  n_sp<-length(table(na.omit(species_table[,3])))
  attr(species_table,"S")[length(attr(species_table,"S"))]<-n_sp
  par(mfrow=c(1,3),mar=c(5.1, 4.1, 9, 2.1), cex=1.25)
  plotplay(species_table,pool_table,title2="Death", n_sp = n_sp, cycles=cycles, main=main, xdim=xdim,ydim=ydim)
  points(species_table[ which(is.na(species_table[,3])),1:2], pch=  16, cex=1.5, col="black")
  plotrichness(species_table,pool_table,title2="Start community", cycles=cycles, main=main)
  plotrank(species_table,pool_table,title2="Start community", cycles=cycles, main=main)
  attr(species_table,"seed")<-seed0
  return(species_table)

}

#' @export
selection<-function(species_table,pool_table,selprob, cycles,main, xdim, ydim, n_pool)
{

  repeat {

    sampled.rows<-which(selprob>=0.5)
    #sampled.rows<-  sample(rep(c(T,F), each=length(selprob)),length(na.omit(species_table[,3])),prob=c(selprob, 1-selprob))
    if(sum(sampled.rows)!=0){ break}
  }


  deaths.rows<- sampled.rows
  #deaths.rows<-sample(1:nrow(species_table),n_ext)
  species_table[deaths.rows,3:4]<-NA
  n_sp<-length(table(na.omit(species_table[,3])))
  attr(species_table,"S")[length(attr(species_table,"S"))]<-n_sp
  par(mfrow=c(1,3),mar=c(5.1, 4.1, 9, 2.1), cex=1.25)
  plotplay(species_table,pool_table,title2="Selection", n_sp = n_sp, cycles=cycles,main=main, xdim=xdim,ydim=ydim)
  points(species_table[ which(is.na(species_table[,3])),1:2], pch=  16, cex=1.5, col="black")

  plotrichness(species_table,pool_table,title2="Start community", cycles=cycles, main=main)
  plotrank(species_table,pool_table,title2="Start community", cycles=cycles, main=main)


  return(species_table)

}


#' @export
birth<-function(species_table,pool_table,targ, cycles, displim,main, xdim=xdim, ydim=ydim, n_pool=n_pool)
{
  distgrid<-as.matrix(dist(species_table[,1:2]))
  neigh<-which( distgrid[ targ,]>0& distgrid[ targ,]<=1.5)
  nei<-as.numeric(row.names(na.omit(species_table[neigh,])))
  seed0<-attr(species_table,"seed")+1
  set.seed(seed0)
  newborns<-sample(nei, 1)

  if(isFALSE(displim))
  {
    set.seed(seed0)
    newborns<-sample(as.numeric(rownames(na.omit(species_table))), 1)
  }

  species_table[targ,3:4]<-species_table[newborns,3:4]
  n_sp<-length(table(na.omit(species_table[,3])))
  attr(species_table,"S")[length(attr(species_table,"S"))]<-n_sp
  par(mfrow=c(1,3),mar=c(5.1, 4.1, 9, 2.1), cex=1.25)
  #grid
  {
    plotplay(species_table,pool_table,title2="Birth", n_sp = n_sp, cycles=cycles,main=main, xdim=xdim,ydim=ydim)
    points(species_table[ which(is.na(species_table[,3])),1:2], pch=  16, cex=1.5, col="black")
    origin<-unlist(species_table[newborns,1:2])
    target<-unlist(species_table[targ,1:2])
    arrows(origin[1],origin[2],target[1],target[2], lty=2, length=.1, code=2)}

  plotrichness(species_table,pool_table,title2="Start community", cycles=cycles, main=main)
  plotrank(species_table,pool_table,title2="Start community", cycles=cycles, main=main)



  attr(species_table,"seed")<-seed0
  return(species_table)


}

#' @export
immig<-function(species_table,pool_table,targ, cycles,main, xdim=xdim, ydim=ydim, n_pool=n_pool)
{

  pool=pool_table[,3]
  cols.pool=pool_table[,4]
  seed0<-attr(species_table,"seed")+1
  set.seed(seed0)
  immig<-sample(1:nrow(pool_table),1)
  species_table[targ,3:4]<-pool_table[immig,3:4]

  im<-pool_table[immig,3]
  n_sp<-length(table(na.omit(species_table[,3])))
  attr(species_table,"S")[length(attr(species_table,"S"))]<-n_sp

  par(mfrow=c(1,3),mar=c(5.1, 4.1, 9, 2.1), cex=1.25)

  #grid
  {plotplay(species_table,pool_table,title2="Immig",n_sp=n_sp, cycles=cycles,main=main, xdim=xdim,ydim=ydim)
    points(species_table[ which(is.na(species_table[,3])),1:2], pch=  16, cex=1.5, col="black")
  colimmig<-  pool_table[immig,4]

  y<-c(mean(range(pool_table[,2])))

  y<-seq(min(species_table[,2])-2,max(species_table[,2])+4,length.out=(nrow(pool_table)))
  y<-y[immig]
  target<-unlist(species_table[targ,1:2])
  par(xpd=T)
  origin<-c(-1.5,  y)
  arrows(origin[1],origin[2],target[1],target[2], lwd=1.5,lty=1, length = 0.15, col="gray50")
  par(xpd=T)
  points(-1.5,  y, pch=im, col=colimmig, cex=2)
  par(xpd=F)}

  plotrichness(species_table,pool_table,title2="Start community", cycles=cycles, main=main)
  plotrank(species_table,pool_table,title2="Start community", cycles=cycles, main=main)

  attr(species_table,"seed")<-seed0



  return(species_table)
}

#' @export
disp<-function(species_table, pool_table,targ, cycles,main, xdim=xdim, ydim=ydim, n_pool=n_pool)
{
  distgrid<-as.matrix(dist(species_table[,1:2]))
  neigh<-which( distgrid[ targ,]>0& distgrid[ targ,]<=1.5)
  nei<-as.numeric(row.names(na.omit(species_table[neigh,1:2])))
  seed0<-attr(species_table,"seed")+1
  set.seed(seed0)
  spmove<-sample(nei,1)
  species_table[targ,3:4]<-  species_table[spmove,3:4]
  species_table[spmove,3:4]<-NA

  n_sp<-length(table(na.omit(species_table[,3])))
  attr(species_table,"S")[length(attr(species_table,"S"))]<-n_sp
  par(mfrow=c(1,3),mar=c(5.1, 4.1, 9, 2.1), cex=1.25)
  plotplay(species_table,pool_table,title2="Dispersal", n_sp = n_sp, cycles=cycles,main=main, xdim=xdim,ydim=ydim)
  points(species_table[ which(is.na(species_table[,3])),1:2], pch=  16, cex=1.5, col="black")
  origin<-unlist(species_table[spmove,1:2])
  target<-unlist(species_table[targ,1:2])
  arrows(origin[1],origin[2],target[1],target[2], lty=2, length=.1, code=2)

   plotrichness(species_table,pool_table,title2="Start community", cycles=cycles, main=main)
  plotrank(species_table,pool_table,title2="Start community", cycles=cycles, main=main)

  attr(species_table,"seed")<-seed0

  return(species_table)


}

#' @export
simulate<-function(species_table,pool_table, cycles, actions=c("death","birth","immig","disp"),prob=NULL,env_selection=F,random_death=T,sel=100,  bd=.505, path=NULL,displim=T, main="anima_theory", record=T, xdim, ydim, n_pool)
{
  seed0<-attr(species_table,"seed")+1
  if(is.null(prob)){
    die.prob<-0.1
    prob.actions<-NULL
  } else {   die.prob<-prob[1]
  prob.actions<-prob[-1]
  }


  time<-attr(species_table,"S")
  for (t in 1:cycles)
  {

    if(isTRUE(env_selection)){

      {


        fitness=pool_table[,"SO"]
        names(fitness)<-pool_table[,3]
        spcom<-species_table[,3]
        names(spcom)<-species_table[,3]
        fit<-(as.vector(fitness[names(spcom)]))
        species_table[,'SO']<-fit
        probdie<-get_selection(species_table,pool_table, sel=sel)*bd
        species_table[,'probdie']<-probdie
        if(is.null(path)==F){   png(paste0(path,recordpic(path)) , 900, 500,pointsize = 13)}
        par(mfrow=c(1,3),mar=c(5.1, 4.1, 9, 2.1), cex=1.25)
        species_table<-selection(species_table,pool_table, selprob=probdie, cycles=cycles,main=main, xdim=xdim, ydim=ydim,n_pool=n_pool)

        if(is.null(path)==F){   graphics.off()}
        if(is.null(path)){if(isTRUE(record)){ani.pause()}}


      }

    }
    if(random_death==T){
      if(is.null(path)==F){   png(paste0(path,recordpic(path)) , 900, 500,pointsize = 13)}
      species_table<-die(species_table,pool_table, prob=die.prob, cycles=cycles,main=main, xdim=xdim, ydim=ydim, n_pool=n_pool)
      if(is.null(path)==F){   graphics.off()}
 if(is.null(path)){if(isTRUE(record)){ani.pause()}}
    }

    empty_spaces=which(is.na(species_table[,3]))




    repeat {
      for( e in 1:length(empty_spaces))
      {
        targ<-empty_spaces[e]

        set.seed(seed0)
        ac<-sample(actions[-1],1, prob=prob.actions)
        seed0<-seed0+1

        if(ac=="birth"){
          distgrid<-as.matrix(dist(species_table[,1:2]))
          neigh<-which( distgrid[ targ,]>0& distgrid[ targ,]<=1.5)


          if(sum( is.na(species_table[neigh,3]))==length(neigh))
          {
            if(is.null(path)==F){   png(paste0(path,recordpic(path)) , 900, 500,pointsize = 13)}
            species_table<-immig(species_table,pool_table,targ, cycles=cycles,main=main, xdim=xdim, ydim=ydim, n_pool=n_pool)
            if(is.null(path)==F){   graphics.off()}
       if(is.null(path)){if(isTRUE(record)){ani.pause()}}
          } else {
            if(is.null(path)==F){   png(paste0(path,recordpic(path)) , 900, 500,pointsize = 13)}
            species_table<-birth(species_table,pool_table,targ, cycles=cycles,displim=displim,main=main, xdim=xdim, ydim=ydim, n_pool=n_pool)
            if(is.null(path)==F){   graphics.off()}
       if(is.null(path)){if(isTRUE(record)){ani.pause()}}
          }



        }
        if(ac=="immig"){
          if(is.null(path)==F){   png(paste0(path,recordpic(path)) , 900, 500,pointsize = 13)}
          species_table<-immig(species_table,pool_table,targ, cycles=cycles,main=main, xdim=xdim, ydim=ydim, n_pool=n_pool)
          if(is.null(path)==F){   graphics.off()}
     if(is.null(path)){if(isTRUE(record)){ani.pause()}}

        }
        if(ac=="disp"){
          temp<-species_table
          if(is.null(path)==F){   png(paste0(path,recordpic(path)) , 900, 500,pointsize = 13)}
          species_table<-disp(species_table,pool_table,targ, cycles=cycles,main=main, xdim=xdim, ydim=ydim, n_pool=n_pool)
          if(is.null(path)==F){   graphics.off()}
     if(is.null(path)){if(isTRUE(record)){ani.pause()}}


        }


      }
      empty_spaces<-which(is.na(species_table[,3])==T)
      if(length(empty_spaces)==0){
        break }
      Sys.sleep(0.05)
    }


    time[length(time)+1]<- length(table(na.omit(species_table[,3])))
    attr(species_table,"S")<-time

  }
  attr(species_table,"seed")<-seed0
  par(mfrow=c(1,3),mar=c(5.1, 4.1, 9, 2.1), cex=1.25)
  plotplay(species_table,pool_table,title2="Final community", cycles=cycles, main=main, xdim=xdim,ydim=ydim)

  plotrichness(species_table,pool_table,title2="Final community", cycles=cycles, main=main)
  plotrank(species_table,pool_table,title2="Final community", cycles=cycles, main=main)


  return(species_table)
}

#' @export
recordpic<-function(path=getwd()){

  if(length(list.files(path))==0){
    filename<-"1.png"
  } else {
    pics<-list.files(path,".png")
    filename<-paste0( max(as.numeric(gsub(".png","",pics)))+1,".png")


  }
  return(filename)
}



#' @export
makeenvi<-function(npatchs=5,species_table,pool_table)
{

  seed0<-attr(species_table,"seed")+1
  set.seed(seed0)
  starts<-sample(1:nrow(species_table),npatchs)
  seed0<-seed0+1
  species_temp=species_table
  species_temp[as.numeric(rownames(species_table))[-starts],3:4]<-NA

  sps.init<-levels(as.factor(species_temp[,3]))
  distgrid<-as.matrix(dist(species_temp[,1:2]))
  repeat{

    for(i in 1:npatchs)
    {
      tobirth<-which(species_temp[,3]== sps.init[i])
      neigh<-which( distgrid[ tobirth,]>0& distgrid[ tobirth,]<=1.5)
      if(length(tobirth)>1)
      {
        neigh<-which( distgrid[ tobirth,]>0& distgrid[ tobirth,]<=1.5,arr.ind = T)[,2]
        neigh<- neigh[which(neigh%in%which(is.na(species_temp[,3])))]
      }

      if(length(neigh%in%which(is.na(species_temp[,3])))>0)

      {
        set.seed(seed0)
      targ<-sample(neigh,1)
      seed0<-seed0+1
      species_temp[targ,3:4]<-unique(species_temp[tobirth,3:4])
      }

    }
    empty_spaces=sum(is.na(species_temp[,3]))
    if(empty_spaces==0){
      break
    }
  }
  ev<-as.numeric(as.factor(species_temp[,3]))
  species_table[,"EV"]<-ev

  quantile(1:nrow(pool_table))






  cols=    matlab.like(nrow(pool_table))[round(seq(1,nrow(pool_table), length.out= nrow(pool_table)/npatchs))][as.factor(species_table[,5])]

  species_table[,"colEV"]<-cols


  attr(species_table,"seed")<-seed0
  return(species_table)
}

#' @export
get_selection<-function(species_table,pool_table, sel=1.3)
{

  fitness=pool_table[,"SO"]
  names(fitness)<-pool_table[,3]
  spcom<-species_table[,3]
  spcom[which(is.na(spcom))]<-0
  names(spcom)<-species_table[,3]
  fit<-as.vector(fitness[names(spcom)])
  fit[which(is.na(fit))]<-0
  fd<-(abs(fit-species_table[,5]))

  probdie=(fd+sel)/(max(fd)+sel)



  return(probdie)
}


#' @export
maketheory<-function(filename,n_sp=25,cycles=100, build_envi = F,npatchs=8,  actions=c("death","birth","immig","disp"),prob=c(0.01,0.6,0.4,0.1), env_selection = F,random_death=T, sel=.3,displim = T,main="anima_theory",xdim=5,ydim=5, n_pool=35, seed0){

  path=getwd()
  path<-paste0(path,"/temp","/")
  file.remove(list.files(path, full.names = F))

  png(paste0(path,recordpic(path)) ,900,500,pointsize = 13)
  init<-initiate(n_sp=n_sp,xdim=xdim,ydim=ydim, n_pool=n_pool,emoj_cat="Animals & Nature", cycles=cycles, build_envi=build_envi, npatchs=npatchs,sel=sel, main=main, seed0 = seed0)
  graphics.off()
  species_table<-init$sp
  pool_table<-init$pool

  attr(species_table,"seed")<-1

  species_table<-simulate(species_table,pool_table, cycles=cycles, actions=c(actions),prob=c(prob), env_selection =env_selection,random_death=random_death, sel=sel,path=NULL,displim = displim, xdim=xdim, ydim=dim, n_pool=n_pool)




  imgs <- list.files(path, full.names = F)
  ## list file names and read in
  imgs <- paste0(path,sort(as.numeric(gsub(".png","",imgs))),".png")
  img_list <- lapply(imgs, image_read)
  ## join the images together
  img_joined <- image_join(img_list)
  ## animate at 2 frames per second
  img_animated <- image_animate(img_joined, fps = 10, optimize = T)
  image_write(image = img_animated,
              path = paste(filename,".gif"))

  file.remove(imgs)
  graphics.off()}




#' @export
simulatheory<-function(n_sp=25,cycles=100, build_envi = T,npatchs=8,  actions=c("death","birth","immig","disp"),prob=c(0.01,0.6,0.4,0.1), env_selection = T,random_death=F, sel=.3,displim = T,path=NULL,main="anima_theory",xdim=5,ydim=5, n_pool=35, seed0){

  init<-initiate(n_sp=n_sp,xdim=xdim,ydim=ydim, n_pool=n_pool,emoj_cat="Animals & Nature", cycles=cycles, build_envi=build_envi, npatchs=npatchs,sel=sel, main=main, seed0 = seed0)
  species_table<-init$sp
  pool_table<-init$pool
  species_table<-simulate(species_table,pool_table, cycles=cycles, actions=c(actions),prob=c(prob), env_selection =env_selection,random_death=random_death, sel=sel,path=path,displim = displim, main=main, xdim=xdim, ydim=dim, n_pool=n_pool)
  }

#' @export
animatheory<-function(title,n_sp=25,cycles=100, build_envi = F,npatchs=8, actions=c("death","birth","immig","disp"),prob=c(0.01,0.6,0,0), env_selection = F,random_death=T, sel=.3,displim = F, path=NULL,main="anima_theory", xdim, ydim, n_pool, seed0) {
  des="\n\n\n\n"
  saveHTML({


    ani.options(verbose=F)
    simulatheory(n_sp=n_sp,cycles=cycles, build_envi = build_envi,npatchs=npatchs, actions=c(actions),prob=c(prob), env_selection = env_selection,random_death=random_death, sel=sel,displim = displim, path=path, main=main, xdim=xdim,ydim=ydim, n_pool=n_pool, seed0=seed0)

  }, img.name = title, imgdir = 'title', htmlfile = paste0('title',".html"),
  autobrowse = T,interval = 0.2,ani.height = 500, ani.width = 1000, loop=1, autoplay = FALSE,verbose=F,description = des, loop=F)


}



#' @export
animatheory<-function(title,n_sp=25,cycles=100, build_envi = F,npatchs=8, actions=c("death","birth","immig","disp"),prob=c(0.01,0.6,0,0), env_selection = F,random_death=T, sel=.3,displim = F, path=NULL,main="anima_theory", xdim, ydim, n_pool, seed0) {
  des="\n\n\n\n"
  saveHTML({


    ani.options(verbose=F)
    simulatheory(n_sp=n_sp,cycles=cycles, build_envi = build_envi,npatchs=npatchs, actions=c(actions),prob=c(prob), env_selection = env_selection,random_death=random_death, sel=sel,displim = displim, path=path, main=main, xdim=xdim,ydim=ydim, n_pool=n_pool, seed0=seed0)

  }, img.name = title, imgdir = title, htmlfile = paste0(title,".html"),
  autobrowse = T,interval = 0.2,ani.height = 500, ani.width = 1000, loop=1, autoplay = FALSE,verbose=F,description = des, loop=F)


}




