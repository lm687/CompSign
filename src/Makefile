all: fullRE_ME_singlelambda_dirichletmultinomial.so diagRE_dirichletmultinomial_single_lambda.so diagRE_ME_dirichletmultinomial.so diagRE_ME_dirichletmultinomialreducedmodel.so diagRE_ME_multinomial.so FE_dirichletmultinomial_single_lambda_copy.so FE_dirichletmultinomial_single_lambda.so FE_dirichletmultinomial.so FE_multinomial.so fullRE_dirichletmultinomial_single_lambda_REv2.so fullRE_dirichletmultinomial_single_lambda.so fullRE_dirichletmultinomial_single_lambda2.so fullRE_ME_dirichletmultinomial_onefixedlambda.so fullRE_ME_dirichletmultinomial_onefixedlambda2.so fullRE_ME_dirichletmultinomial_onefixedlambda3.so fullRE_ME_dirichletmultinomial_sparsecov.so fullRE_ME_dirichletmultinomial_sparsecov2.so fullRE_ME_dirichletmultinomial.so fullRE_ME_dirichletmultinomialnoscaling.so fullRE_ME_halfdirichletmultinomial.so fullRE_ME_multinomial.so ME_dirichletmultinomial_single_lambda.so singleRE_dirichlet_multinomial.so

fullRE_ME_singlelambda_dirichletmultinomial.so: fullRE_ME_singlelambda_dirichletmultinomial.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_singlelambda_dirichletmultinomial.cpp')"

diagRE_dirichletmultinomial_single_lambda.so: diagRE_dirichletmultinomial_single_lambda.cpp
	Rscript --vanilla -e "TMB::compile('diagRE_dirichletmultinomial_single_lambda.cpp')"

diagRE_ME_dirichletmultinomial.so: diagRE_ME_dirichletmultinomial.cpp
	Rscript --vanilla -e "TMB::compile('diagRE_ME_dirichletmultinomial.cpp')"

diagRE_ME_dirichletmultinomialreducedmodel.so: diagRE_ME_dirichletmultinomialreducedmodel.cpp
	Rscript --vanilla -e "TMB::compile('diagRE_ME_dirichletmultinomialreducedmodel.cpp')"

diagRE_ME_multinomial.so: diagRE_ME_multinomial.cpp
	Rscript --vanilla -e "TMB::compile('diagRE_ME_multinomial.cpp')"

FE_dirichletmultinomial_single_lambda_copy.so: FE_dirichletmultinomial_single_lambda_copy.cpp
	Rscript --vanilla -e "TMB::compile('FE_dirichletmultinomial_single_lambda_copy.cpp')"

FE_dirichletmultinomial_single_lambda.so: FE_dirichletmultinomial_single_lambda.cpp
	Rscript --vanilla -e "TMB::compile('FE_dirichletmultinomial_single_lambda.cpp')"

FE_dirichletmultinomial.so: FE_dirichletmultinomial.cpp
	Rscript --vanilla -e "TMB::compile('FE_dirichletmultinomial.cpp')"

FE_multinomial.so: FE_multinomial.cpp
	Rscript --vanilla -e "TMB::compile('FE_multinomial.cpp')"

fullRE_dirichletmultinomial_single_lambda_REv2.so: fullRE_dirichletmultinomial_single_lambda_REv2.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_dirichletmultinomial_single_lambda_REv2.cpp')"

fullRE_dirichletmultinomial_single_lambda.so: fullRE_dirichletmultinomial_single_lambda.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_dirichletmultinomial_single_lambda.cpp')"

fullRE_dirichletmultinomial_single_lambda2.so: fullRE_dirichletmultinomial_single_lambda2.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_dirichletmultinomial_single_lambda2.cpp')"

fullRE_ME_dirichletmultinomial_onefixedlambda.so: fullRE_ME_dirichletmultinomial_onefixedlambda.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_dirichletmultinomial_onefixedlambda.cpp')"

fullRE_ME_dirichletmultinomial_onefixedlambda2.so: fullRE_ME_dirichletmultinomial_onefixedlambda2.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_dirichletmultinomial_onefixedlambda2.cpp')"

fullRE_ME_dirichletmultinomial_onefixedlambda3.so: fullRE_ME_dirichletmultinomial_onefixedlambda3.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_dirichletmultinomial_onefixedlambda3.cpp')"

fullRE_ME_dirichletmultinomial_sparsecov.so: fullRE_ME_dirichletmultinomial_sparsecov.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_dirichletmultinomial_sparsecov.cpp')"

fullRE_ME_dirichletmultinomial_sparsecov2.so: fullRE_ME_dirichletmultinomial_sparsecov2.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_dirichletmultinomial_sparsecov2.cpp')"

fullRE_ME_dirichletmultinomial.so: fullRE_ME_dirichletmultinomial.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_dirichletmultinomial.cpp')"

fullRE_ME_dirichletmultinomialnoscaling.so: fullRE_ME_dirichletmultinomialnoscaling.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_dirichletmultinomialnoscaling.cpp')"

fullRE_ME_halfdirichletmultinomial.so: fullRE_ME_halfdirichletmultinomial.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_halfdirichletmultinomial.cpp')"

fullRE_ME_multinomial.so: fullRE_ME_multinomial.cpp
	Rscript --vanilla -e "TMB::compile('fullRE_ME_multinomial.cpp')"

ME_dirichletmultinomial_single_lambda.so: ME_dirichletmultinomial_single_lambda.cpp
	Rscript --vanilla -e "TMB::compile('ME_dirichletmultinomial_single_lambda.cpp')"

singleRE_dirichlet_multinomial.so: singleRE_dirichlet_multinomial.cpp
	Rscript --vanilla -e "TMB::compile('singleRE_dirichlet_multinomial.cpp')"

clean:
	rm -rf *o

