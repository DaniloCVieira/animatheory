
rm(list=ls())
#roxygen2::roxygenize()

library(animatheory)
#neu_imm0disp0

main = c("Neutral model:  \n Random death and birth")
title<-'neu_imm0disp0'
animatheory(title,n_sp=50,cycles=200, build_envi = F,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(0.1,0.6,0,0), env_selection = F,random_death=T, sel=.3,displim = F, path=NULL, main=main, xdim=7,ydim=7,n_pool=100)



#neu_imm0dispT
title<-'neu_imm0dispT'
main = c("Neutral model:  \n Random death, birth and dispersal limitation")
set.seed(1)
animatheory(title,n_sp=50,cycles=200, build_envi = F,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(0.1,0.6,0,0), env_selection = F,random_death=T, sel=.3,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=100)

#neu_imm0dispTmove
title<-'neu_imm0dispTmove'
main = c("Neutral model:  \n Random death, birth and dispersal limitation (move)")
set.seed(1)
animatheory(title,n_sp=50,cycles=200, build_envi = F,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(0.1,0.6,0,0.3), env_selection = F,random_death=T, sel=.3,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=100)



#neu_immTdispT
title<-'neu_immTdispT'
main = c("Neutral model:  \n Random death, birth, dispersal limitation and immigration")
set.seed(1)
animatheory(title,n_sp=50,cycles=200, build_envi = F,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(0.1,0.6,0.4,0), env_selection = F,random_death=T, sel=.3,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=100)

#neu_immTdispTmove
title<-'neu_immTdispTmove'
main = c("Neutral model:  \n Random death, birth, dispersal limitation and immigration (move)")
set.seed(1)
animatheory(title,n_sp=50,cycles=200, build_envi = F,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(0.1,0.6,0.4,0.3), env_selection = F,random_death=T, sel=.3,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=100)




#sel_imm0disp0
main = c("Environmental filtering:  \n Random death and birth")
title<-'sel_imm0disp0'
set.seed(1)
animatheory(title,n_sp=50,cycles=200, build_envi = T,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(0.1,0.6,0,0), env_selection = T,random_death=F, sel=.5,displim = F, path=NULL, main=main, xdim=7,ydim=7,n_pool=100)


#sel_imm0dispT
title<-'sel_imm0dispT'
main = c("Environmental filtering:  \n Random death, birth and dispersal limitation")
set.seed(1)
animatheory(title,n_sp=50,cycles=200, build_envi = T,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(0.1,0.6,0,0), env_selection = T,random_death=F, sel=.5,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=100)


#sel_imm0dispTmove
title<-'sel_imm0dispTmove'
main = c("Environmental filtering:  \n Random death, birth and dispersal limitation (move)")
set.seed(1)
animatheory(title,n_sp=50,cycles=200, build_envi = T,npatchs=6, actions=c("death","birth","immig","disp"),prob=c(0.1,0.6,0,0.3), env_selection = T,random_death=F, sel=.5,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=100)


#source("anima_theory.R")
#sel_immTdispT
title<-'sel_immTdispT'
main = c("Environmental filtering: \n Random death, birth, dispersal limitation and immigration")
set.seed(1)
animatheory(title,n_sp=50,cycles=200, build_envi = T,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(0.1,0.6,0.4,0), env_selection = T,random_death=F, sel=.5,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=100)

#sel_immTdispTmove
title<-'sel_immTdispTmove'
main = c("Environmental filtering: \n Random death, birth, dispersal limitation and immigration (move)")
set.seed(1)
animatheory(title,n_sp=50,cycles=200, build_envi = T,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(0.1,0.6,0.4,0.3), env_selection = T,random_death=F, sel=.5,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=100)



