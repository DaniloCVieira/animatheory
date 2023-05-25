source('anima_theory4.R')
#rm(list=ls())
#roxygen2::roxygenize()

#neu_imm0disp0

main = c("Neutral model:  Random death and birth")

title<-'neu_imm0disp0'
animatheory(title,n_sp=10,cycles=60, build_envi = F,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(10,0.6,0,0), env_selection = F,random_death=T, sel=.3,displim = F, path=NULL, main=main, xdim=7,ydim=7,n_pool=25, seed0=5)



#neu_imm0dispT
title<-'neu_imm0dispT'
main = c("Neutral model:  Random death, birth and dispersal limitation")
set.seed(1)
animatheory(title,n_sp=10,cycles=60, build_envi = F,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(10,0.6,0,0), env_selection = F,random_death=T, sel=.3,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=25, seed0=5)


#neu_immTdispT
title<-'neu_immTdispT'
main = c("Neutral model:  Random death, birth, dispersal limitation and immigration")
set.seed(1)
animatheory(title,n_sp=10,cycles=60, build_envi = F,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(10,0.6,0.05,0), env_selection = F,random_death=T, sel=.3,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=25, seed0=5)

##############################################################################



#sel_imm0disp0
main = c("Environmental filtering and birth")
title<-'sel_imm0disp0'
set.seed(1)
animatheory(title,n_sp=10,cycles=100, build_envi = T,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(10,0.6,0,0), env_selection = T,random_death=F, sel=20,displim = F, path=NULL, main=main, xdim=7,ydim=7,n_pool=25, seed0=5, selforce = 10)


#sel_imm0dispT
title<-'sel_imm0dispT'
main = c("Environmental filtering, birth and dispersal limitation")
set.seed(1)
animatheory(title,n_sp=10,cycles=100, build_envi = T,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(10,0.6,0,0), env_selection = T,random_death=F, sel=20,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=25, seed0=5, selforce=10)


#sel_imm0dispTmove
title<-'sel_imm0dispTmove'
main = c("Environmental filtering, birth and dispersal limitation (move)")
set.seed(1)
animatheory(title,n_sp=10,cycles=100, build_envi = T,npatchs=6, actions=c("death","birth","immig","disp"),prob=c(10,0.6,0,0.3), env_selection = T,random_death=F, sel=20,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=25, seed0=5, selforce=10)


#source("anima_theory.R")
#sel_immTdispT
title<-'sel_immTdispT'
main = c("Environmental filtering, birth, dispersal limitation and immigration")
set.seed(1)
animatheory(title,n_sp=10,cycles=100, build_envi = T,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(10,0.6,0.4,0), env_selection = T,random_death=F, sel=20,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=25, seed0=5, selforce=10)

#sel_immTdispTmove
title<-'sel_immTdispTmove'
main = c("Environmental filtering, birth, dispersal limitation and immigration (move)")
set.seed(1)
animatheory(title,n_sp=10,cycles=100, build_envi = T,npatchs=5, actions=c("death","birth","immig","disp"),prob=c(10,0.6,0.4,0.3), env_selection = T,random_death=F, sel=20,displim = T, path=NULL, main=main, xdim=7,ydim=7,n_pool=25, seed0=5, selforce=10)



