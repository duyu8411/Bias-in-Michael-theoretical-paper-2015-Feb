library(mvtnorm)
library(fBasics)

########################Function for STD design############################
MC_STD_BIAS_VCOV <- function(pie_1,kstar,vcov,uc,u1,u2,delta_1,delta_2,times){
        
        delta_c <- pie_1*delta_1+(1-pie_1)*delta_2
        ##convert variance covariance matrix to the standardized
        var_vec <- diag(vcov)
        n_row <- nrow(vcov)
        zvcov <- diag(1,n_row)
        for(i in 1:(n_row-1)){
                for(j in (i+1):n_row){
                        zvcov[i,j] <- vcov[i,j]/sqrt(var_vec[i]*var_vec[j])
                        zvcov[j,i] <- zvcov[i,j]
                }
        }
        se_c_k <- sqrt(var_vec[1:kstar])
        se_1_k <- sqrt(var_vec[(kstar+1):(kstar+5)])
        se_2_k <- sqrt(var_vec[(kstar+6):(kstar*2+5)])
        
        
        ##mean vector
        mean_c <- delta_c/se_c_k[1:kstar]
        mean_1 <- delta_1/se_1_k[1:5]
        mean_2 <- delta_2/se_2_k[1:kstar]
        
        mean_vec <- c(mean_c,mean_1,mean_2)
        set.seed(1228)
        dat <- rmvnorm(n=times,mean=mean_vec,sigma=zvcov,method="svd")
        dat <- as.matrix(dat)
        
        lc <- rep(0,kstar)
        lc[kstar] <- uc[kstar]
        l1 <- rep(0,5)
        l1[5] <- u1[5]
        l2 <- rep(0,kstar)
        l2[kstar] <- u2[kstar]
        
        Integrand <- function(vec){
                z_c <- vec[1:kstar]
                z_1 <- vec[(kstar+1):(kstar+5)]
                z_2 <- vec[(kstar+6):(kstar*2+5)]
                for(k in 1:3){
                        if(z_c[k]>uc[k] | z_c[k] <= lc[k]){
                                delta_c_star <- z_c[k]*se_c_k[k]
                                v_c <- abs(delta_c_star-delta_c)/(1.96*se_c_k[k])
                                break
                        }                        
                }
                
                for(k in 1:5){
                        if(z_1[k]>u1[k] | z_1[k] <= l1[k]){
                                delta_1_star <- z_1[k]*se_1_k[k]
                                v_1 <- abs(delta_1_star-delta_1)/(1.96*se_1_k[k])
                                break
                        }
                }
                
                for(k in 1:3){
                        if(z_2[k]>u2[k] | z_2[k] <= l2[k]){
                                delta_2_star <- z_2[k]*se_2_k[k]
                                v_2 <- abs(delta_2_star-delta_2)/(1.96*se_2_k[k])
                                break
                        }
                }
                return(c(delta_1_star,delta_2_star,delta_c_star,delta_1_star-1.96*se_1_k[k],delta_1_star+1.96*se_1_k[k],
                         delta_2_star-1.96*se_2_k[k],delta_2_star+1.96*se_2_k[k],delta_c_star-1.96*se_c_k[k],delta_c_star+1.96*se_c_k[k],v_1,v_2,v_c))
                
                
        }
        
        delta_star_sim <- t(apply(dat,1,Integrand))
        delta_star_sim_3 <- delta_star_sim[,1:3]
        delta <- c(delta_1,delta_2,delta_c)
        bias <- apply(delta_star_sim_3,2,mean)-delta
        se <- apply(delta_star_sim_3,2,sd)
        skewness <- apply(delta_star_sim_3,2,skewness)
        names(bias) <- c("bias_delta_1","bias_delta_2","bias_delta_c")
        names(se) <- c("se_delta_1","se_delta_2","se_delta_c")
        names(skewness) <- c("sk_delta_1","sk_delta_2","sk_delta_c")
        delta_star_mse_tmp <- data.frame(sqerr1=(delta_star_sim[,1]-delta_1)^2,sqerr2=(delta_star_sim[,2]-delta_2)^2,sqerr3=(delta_star_sim[,3]-delta_c)^2)
        mse <- apply(delta_star_mse_tmp,2,mean)
        names(mse) <- c("mse_delta_1","mse_delta_2","mse_delta_c")
        
        cp_1_tmp <- ifelse(delta_star_sim[,4]<=delta_1 & delta_star_sim[,5]>=delta_1,1,0)
        cp_2_tmp <- ifelse(delta_star_sim[,6]<=delta_2 & delta_star_sim[,7]>=delta_2,1,0)
        cp_c_tmp <- ifelse(delta_star_sim[,8]<=delta_c & delta_star_sim[,9]>=delta_c,1,0)
        
        cp_1 <- sum(cp_1_tmp)/length(cp_1_tmp)
        cp_2 <- sum(cp_2_tmp)/length(cp_2_tmp)
        cp_c <- sum(cp_c_tmp)/length(cp_c_tmp)
        
        names(cp_1) <- "cp_delta_1"
        names(cp_2) <- "cp_delta_2"
        names(cp_c) <- "cp_delta_c"
        
        infla_c1 <- quantile(delta_star_sim[,10],probs=0.95)
        infla_c2 <- quantile(delta_star_sim[,11],probs=0.95)
        infla_cc <- quantile(delta_star_sim[,12],probs=0.95)
        
        names(infla_c1) <- "inflac_delta_1"
        names(infla_c2) <- "inflac_delta_2"
        names(infla_cc) <- "inflac_delta_c"
        return(c(bias,se,skewness,mse,cp_1,cp_2,cp_c,infla_c1,infla_c2,infla_cc))
        
}

load("trial_result_with_updated_bdry_for_theoretical_paper.rda")
vcov_tmp <- result$cov_ltmle
vcov <- diag(11)
for(i in 1:11){
        if(i<=3 & i>=1){
                vcov[i,1:11] <- vcov_tmp[(3*i-2),c(1,4,7,2,5,8,11,14,3,6,9)] 
        }else if(i>=4 & i<=8){
                vcov[i,1:11] <- vcov_tmp[(3*i-10),c(1,4,7,2,5,8,11,14,3,6,9)]
        }else{
                vcov[i,1:11] <- vcov_tmp[(3*i-24),c(1,4,7,2,5,8,11,14,3,6,9)]
        }
}
rm(vcov_tmp)

uc <- c(3.470268,2.453850,2.003560)
u1 <- c(4.628951,3.273163,2.672526,2.314476,2.070130)
u2 <- c(3.469055,2.452993,2.002860)







############################ADJ-STD#################################
MC_STD_BIAS_VCOV(1/3,3,vcov,uc,u1,u2,0.125,0.125,10000)
MC_STD_BIAS_VCOV(1/3,3,vcov,uc,u1,u2,0.125,0,10000)
MC_STD_BIAS_VCOV(1/3,3,vcov,uc,u1,u2,0,0,10000)

d1 <- seq(-1,1,by=0.01)
d2 <- seq(-1,1,by=0.01)
delta <- expand.grid(d1,d2)
table_bias1 <- t(apply(delta,1,function(a) MC_STD_BIAS_VCOV(1/3,3,vcov,uc,u1,u2,a[1],a[2],10000)))
delta_STD1 <- cbind(delta,table_bias1)
save(delta_STD1,file="delta_Paper_STD_up_adjvc.RData")
##max bias 
apply(table_bias1,2,max)
delta_STD1[which.max(table_bias1[,1]),] ##delta_values for delta1 max bias
delta_STD1[which.max(table_bias1[,2]),] ##delta_values for delta2 max bias
delta_STD1[which.max(table_bias1[,3]),] ##delta_values for deltac max bias

apply(table_bias1,2,min)
delta_STD1[which.min(table_bias1[,1]),] ##delta_values for delta1 min bias
delta_STD1[which.min(table_bias1[,2]),] ##delta_values for delta2 min bias
delta_STD1[which.min(table_bias1[,3]),] ##delta_values for deltac min bias