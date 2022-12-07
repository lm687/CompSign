# require(TMB)
# TMB::compile("R/mm_multinomial/fullRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
# dyn.load(dynlib("R/mm_multinomial/fullRE_ME_dirichletmultinomial"))

# .onLoad <- function(lib) {
#   cat("Loading compiled code...\n")
#   dyn.load(dynlib("/singleRE_dirichlet_multinomial"))
#   # library.dynam("diagRE_ME_dirichletmultinomial", "Compsign")
# }
# 
# .onLoad()

dummyfun <- function(){cat('hehe\n')}
