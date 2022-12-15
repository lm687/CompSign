all: diagRE_dirichletmultinomial_single_lambda.so diagRE_ME_dirichletmultinomial.so diagRE_ME_multinomial.so FE_dirichletmultinomial_single_lambda.so FE_dirichletmultinomial.so fullRE_dirichletmultinomial_single_lambda.so fullRE_ME_dirichletmultinomial_onefixedlambda.so fullRE_ME_dirichletmultinomial.so fullRE_ME_multinomial.so singleRE_dirichlet_multinomial.so


diagRE_dirichletmultinomial_single_lambda.so: diagRE_dirichletmultinomial_single_lambda.cpp
	Rscript --vanilla -e "TMB::compile('diagRE_dirichletmultinomial_single_lambda.cpp')"

diagRE_ME_dirichletmultinomial.so: diagRE_ME_dirichletmultinomial.cpp
	Rscript --vanilla -e "TMB::compile('diagRE_ME_dirichletmultinomial.cpp')"

diagRE_ME_multinomial.so: diagRE_ME_multinomial.cpp
	Rscript --vanilla -e "TMB::compile('diagRE_ME_multinomial.cpp')"

FE_dirichletmultinomial_single_lambda.so: FE_dirichletmultinomial_single_lambda.cpp
	Rscript --vanilla -e "TMB::compile('FE_dirichletmultinomial_single_lambda.cpp')"

FE_dirichletmultinomial.so: FE_dirichletmultinomial.cpp
	Rscript --vanilla -e "TMB::compile('FE_dirichletmultinomial.cpp')"

fullRE_dirichletmultinomial_single_lambda.so: fullRE_dirichletmultinomial_single_lambda.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_dirichletmultinomial_single_lambda.cpp')"

fullRE_ME_dirichletmultinomial_onefixedlambda.so: fullRE_ME_dirichletmultinomial_onefixedlambda.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_dirichletmultinomial_onefixedlambda.cpp')"

fullRE_ME_dirichletmultinomial.so: fullRE_ME_dirichletmultinomial.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_dirichletmultinomial.cpp')"

fullRE_ME_multinomial.so: fullRE_ME_multinomial.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_multinomial.cpp')"

singleRE_dirichlet_multinomial.so: singleRE_dirichlet_multinomial.cpp
	Rscript --vanilla -e "TMB::compile('singleRE_dirichlet_multinomial.cpp')"

clean:
	rm -rf *o

