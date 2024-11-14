# Test
rm(list = ls())
library(susieR);library(mr.ash.alpha)
source("New_Setting_10192024.R");source("susie_ash_joint_ELBO_v3 (1).R")
source("SuSiE-ASH/code/susie_versions/susie_inf.R")

X20 <- readRDS("G_block_20_match_ENCODE.RDS")

gen.dat = generate_eqtl_data(X20, seed = 54629110) #261568381)#, seed = 
for(i in 1:length(gen.dat)){assign(names(gen.dat)[i], gen.dat[[i]])}
threshold = 0.005
causal_SNPs <- is_causal(gen.dat, threshold)# which(proportion_var_explained > threshold)
# plot(gen.dat$beta)
# 
X = ori.X
y=  ori.y
v_threshold = threshold
true.beta = beta
res_susie_ash_cal = susie_ash(ori.X, ori.y, L = 10, est_var = "cal_v", true_var_res = NULL, v_threshold = threshold)
# res_susie_ash_mom = susie_ash(ori.X, ori.y, L = 10, est_var = "mom", true_var_res = NULL)
# res_susie=  susie(ori.X, ori.y, L = 10)
# res_susie_inf = susie_inf(X = scale(X,T,T), y= scale(y,T,F), L= 10)

### Variance Estimation:
gen.dat$var_epsilon
res_susie_ash_cal$sigma2
# res_susie_ash_mom$sigma2
# res_susie$sigma2
# res_susie_inf$sigmasq


### Check ELBO
plot(res_susie_ash_cal$elbo); res_susie_ash_cal$elbo;all(diff(res_susie_ash_cal$elbo) >= 0);diff(res_susie_ash_cal$elbo)
# plot(res_susie_ash_mom$elbo); res_susie_ash_mom$elbo;  
# plot(res_susie$elbo); res_susie$elbo;  all(diff(res_susie$elbo) >= 0) 

### calculate metrics:
compute_metrics(res_susie_ash_cal$sets$cs, causal_SNPs)
# compute_metrics(res_susie_ash_mom$sets$cs, causal_SNPs)
# compute_metrics(res_susie$sets$cs, causal_SNPs)
# compute_metrics(res_susie_inf$sets, causal_SNPs)


L = 10
scaled_prior_variance = 0.2;
residual_variance = NULL;
prior_weights = NULL;
null_weight = 0;
standardize = TRUE;
intercept = TRUE;
estimate_residual_variance = TRUE;
estimate_prior_variance = TRUE;
estimate_prior_method = "simple";
check_null_threshold = 0;
prior_tol = 1e-9;
residual_variance_upperbound = Inf;
s_init = NULL;
coverage = 0.95;
min_abs_corr = 0.5;
median_abs_corr = NULL;
compute_univariate_zscore = FALSE;
na.rm = FALSE;
max_iter = 100;
tol = 1e-3;
verbose = FALSE;
track_fit = FALSE;
residual_variance_lowerbound = var(drop(y))/1e4;
refine = FALSE;
n_purity = 100;
est_var = "cal_v";
true_var_res = NULL
v_threshold = threshold

