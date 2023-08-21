# setwd("C:/Users/benny/Box/WMM_project/code/ARMM")
# setwd("~/Library/CloudStorage/Box-Box/WMM_project/code/ARMM")

source("0_functions.R")

library(ggplot2)
library(usmap)
library(data.table)
library(xtable)
# install.packages("remotes")
# remotes::install_github("CarolinaEuan/HMClust") # installed on 08-21-2023
library(HMClust) # must be installed with GitHub



# get COVID-19 data from only PA
# covid <- read.csv("us-counties.txt")
# dta_old <- covid[ covid$state %in% c("Pennsylvania"), ]
# write.csv(dta_old, "PA.txt", row.names=FALSE)

dta_old <- read.csv("PA.txt")
# fill in missing dates with zeros
for( state1 in unique(as.character(dta_old$state)) ){
  for( county1 in unique(as.character(dta_old$county)) ){
    county_dta <- dta_old[(dta_old$state==state1)&(dta_old$county==county1),]
    diff <- as.integer( max(as.Date(county_dta$date))-min(as.Date(county_dta$date)) )
    if(!is.na(diff) & nrow(county_dta)!=(1+diff)){
      print(county1)
      d_old <- 0
      for(d in 0:diff){
        if( nrow( county_dta[as.Date(county_dta$date)==(min(as.Date(county_dta$date))+d),] )==0 ){
          print(d)
          dta_old <- rbind(dta_old,county_dta[as.Date(county_dta$date)==(min(as.Date(county_dta$date))+d_old),])
        }else{
          d_old <- d
        }
      }
    }
  }
}
dta_old$date <- as.Date(dta_old$date)

# change data from cumulative to incremental
dta <- NULL
for(state in unique(dta_old$state)){
  for(county in unique(dta_old$county)){
    temp <- dta_old[(dta_old$county==county)&(dta_old$state==state),]
    temp <- temp[order(as.Date(temp$date)),]
    
    temp$cases[-1] <- temp$cases[-1] - temp$cases[-length(temp$cases)]
    temp$cases <- temp$cases*(temp$cases>0)
    
    temp$deaths[-1] <- temp$deaths[-1] - temp$deaths[-length(temp$deaths)]
    temp$deaths <- temp$deaths*(temp$deaths>0)
    
    dta <- rbind(dta,temp)
  }
}

# check there's no missing days
for( county in unique(dta$county) ){
  for( state in unique(dta$state) ){
    temp <- dta[(dta$county==county) & (dta$state==state),]
    n_rows <- nrow(temp)
    if( n_rows > 0){
      for( r in 1:(n_rows-1) ){
        if( ( as.Date(temp$date[r+1]) - as.Date(temp$date[r]) ) != 1 ){
          print(county)
          print(temp[c(r,r+1),])
          print(as.Date(temp$date[r]))
          print(as.Date(temp$date[r+1]))
        }
      }
    }
  }
}
dta <- dta[(dta$date>=as.Date("2020-10-01"))*(dta$date<as.Date("2021-03-01"))==1,]

# get counties mapping
counties <- read.csv("2015_counties_table.csv")
# replace Abb names with full state names
counties$USPS <- as.character( sapply( as.character(counties$USPS), function(x) state.name[state.abb==x]) )
# remove " County" from string
counties$NAME <- as.character( sapply( as.character(counties$NAME), function(x) strsplit(x, " County")[1] ) )
counties$state <- counties$USPS
counties$county <- counties$NAME
# rename lat and long columns
counties$y <- counties$INTPTLAT
counties$x <- counties$INTPTLONG
counties <- counties[c("state","county","GEOID","y","x")]
# add lat and long to PA and NJ data set
covid <- merge(dta, counties)

census <- read.csv("census_county_interpolated.txt", row.names = 1)
census <- census[census$year==2016,]
census <- census[ , colnames(census)[c(-2,-(12:14))] ]
colnames(census) <- c( "GEOID", colnames(census)[-1] )
dta <- merge(covid, census)
dta$cases <- dta$cases/dta$population

X <- c()
for(county in unique(dta$county) ){
  X <- cbind(X,dta[dta$county==county,]$cases-
               mean(dta[dta$county==county,]$cases) ) 
}
X <- data.frame(X)
colnames(X) <- unique(dta$county)

coeff_df <- NULL
matrix_list <- list() # will be used for clustering
for(i in 1:dim(X)[2]){
  matrix_list[[i]] <- dim(X)[1]*
    cov_to_matrix(acf(X[,i], lag.max = 7, plot = F, 
                      demean = F, type = "covariance")$acf)
  #
  coeff <- arima(X[,i],order = c(7,0,0),method = "ML",include.mean = F)$coef[1:7]
  coeff_df <- rbind(coeff_df, coeff)
}
n_vector <- rep(dim(X)[1],dim(X)[2])

Clust1 <- HMClust::HSM(X) # use to initialize groups

colnames(X)[HMClust::cutk(Clust1,3)[[1]]]
colnames(X)[HMClust::cutk(Clust1,3)[[2]]]
colnames(X)[HMClust::cutk(Clust1,3)[[3]]]

# calculate log-lik and BIC
ll <- c()
for(G in 1:10){
  tmp <- HMClust::cutk(Clust1,G)
  
  z_list_init <- list()
  for(g in 1:G){
    z_list_init[[g]] <- 1*((1:dim(X)[2]) %in% tmp[[g]])
  }
  
  out <- WMM(matrix_list=matrix_list, n_vector = n_vector, G = G,
             method = "EM1", upper=upper,
             burn_in = 10, max_iter = 1000,
             pi_list_init = NULL, sigma_list_init = NULL, 
             z_list_init = z_list_init, 
             lambda_list_init = NULL)
  #
  tmp <- logLik_ARMM(out$z_list,matrix_list,out$sigma_list,n_vector)
  ll <- c(ll,tmp)
  print(G)
}

m <- (1:10)*8-1
ll + 2*m
ll + m*log( sum(n_vector) )

which.min(ll + 2*m)
which.min(ll + m*log( sum(n_vector) ))

# plot AIC/BIC data
jpeg("BIC.jpg",width = 1500, height = 1000)
par(mfrow=c(1,1))
par(mar=c(5, 7, 2, 2))
plot(1:10,ll,type="l",
     col="black", lwd=3,
     xlab="",ylab="", pch=16,
     cex.lab=2, cex.axis=2, cex.main=3, cex=2)
lines(1:10,ll + m*log( sum(n_vector) ), 
      lwd=3, col="blue",lty=5)
mtext("Log-likelihood/BIC", cex=5, side=2, line=3)
mtext("G", cex=5, side=1, line=3)
legend(7, -163000, legend=c("Log-likelihood", "BIC"),
       lty=c(1,2,5), 
       col=c("black","blue"), 
       lwd = c(2,2,2), cex=3)
dev.off()

#
set.seed(1)
G=5
# GMM results
GMM <- Mclust(data = coeff_df, G=G)
colnames(X)[ apply(GMM$z, 1, which.max)==5 ]
tmp <- HMClust::cutk(Clust1,G)

z_list_init <- list()
for(g in 1:G){
  z_list_init[[g]] <- 1*((1:dim(X)[2]) %in% tmp[[g]])
}

# WMM output ####
out <- WMM(matrix_list=matrix_list, n_vector = n_vector, G = G,
           method = "EM1", upper=upper,
           burn_in = 10, max_iter = 1000,
           pi_list_init = NULL, sigma_list_init = NULL, 
           z_list_init = z_list_init, 
           lambda_list_init = NULL)
#
PA_groups <- data.frame( z_max(out$z_list) )
PA_groups <- data.frame( "z_WMM" = apply(PA_groups, MARGIN = 1, which.max) )
#PA_groups$z_WMM <- PA_groups$z_WMM - 1 
#PA_groups$z_WMM[PA_groups$z_WMM==0] <- 3
PA_groups$z_WMM <- as.factor( PA_groups$z_WMM )
PA_groups$county <- paste( colnames(X), "County" )
#
tmp <- data.frame( z_max(z_list_init) )
tmp <- apply(tmp, MARGIN = 1, which.max)
PA_groups$z_HSM <- tmp 
#PA_groups$z_HSM <- PA_groups$z_HSM - 1
#PA_groups$z_HSM[PA_groups$z_HSM==0] <- 3
PA_groups$z_HSM <- as.factor( PA_groups$z_HSM ) 

#

PA_counties <- countypop[countypop$abbr=="PA",]
PA_counties <- merge(PA_counties,PA_groups)

p1 <- plot_usmap(data = countypop, values = "pop_2015", "counties",
           include = c("PA")) + 
  scale_fill_continuous(
    low = "white", high = "black", name = "", label = scales::comma
  ) + theme(legend.position = "none",
            legend.text=element_text(size=7),
            title=element_text(size=40, face='bold'),
            plot.margin = margin(.5, .5, .5, .5, "cm"),
            plot.background = element_rect(
              colour = "black",
              size = 1
            )) +
  labs(title = "Population")

p2 <- plot_usmap(data = PA_counties, values = "z_WMM", "counties",
           include = c("PA")) + 
  scale_colour_steps() + labs(fill = "") +
  theme(legend.position = "right",
        legend.text=element_text(size=50),
        title=element_text(size=40, face='bold'),
        plot.margin = margin(1.5, .5, 1.5, .5, "cm"),
        plot.background = element_rect(
          colour = "black",
          size = 1
        )) +
  labs(title = "WMM")
p3 <- plot_usmap(data = PA_counties, values = "z_HSM", "counties",
           include = c("PA")) + 
  scale_colour_steps() + labs(fill = "") +
  theme(legend.position = "right",
        legend.text=element_text(size=50),
        title=element_text(size=40, face='bold'),
        plot.margin = margin(1.5, .5, 1.5, .5, "cm"),
        plot.background = element_rect(
          colour = "black",
          size = 1
        )) +
  labs(title = "HSM")

# plot PA map data ####
jpeg("PA_plot.jpg",width = 1500, height = 500)
cowplot::plot_grid(p1, p2, p3, nrow = 1)
dev.off()

table_latex <- NULL
for(k in 1:G){
  table_latex <- cbind(table_latex,
    paste0(round(as.numeric(tbl[[k]]$coeff),4),"(",
           round(diag( tbl[[k]]$var )**.5,4),")") )
}
print(xtable(table_latex, digits = 4, type = "latex"), file = "PA_table.tex")


# fit/plot forecast, minor revision 02/20/2022 ####
PA_groups[PA_groups$z_HSM==1,]
county_mean = mean(dta[dta$county=="Centre",]$cases)
case_rates = X$Centre 

dta_obs = data_matrix(case_rates,6)
dta_obs = dta_obs[complete.cases(dta_obs),]

pred_case_rates = as.matrix(dta_obs) %*% tbl[[3]]$coeff

mat_obs = matrix_list[[which(names(X)=="Centre")]]

sigma_hat = out$sigma_list[[3]]
schur = 1 - sigma_hat[2:7,1] %*% solve(sigma_hat[2:7,2:7]) %*% sigma_hat[1,2:7] / sigma_hat[1,1]
innovation_var = mat_obs[1,1]/dim(X)[1]*schur

se = sqrt( diag( dta_obs %*% tbl[[3]]$var %*% t(dta_obs) ) + 
             as.numeric(innovation_var) )

# fit stand alone ar model
yw_fit = ar.mle(case_rates,demean = F,aic = F,order.max = 7)
pred_case_rates_new = as.matrix(dta_obs) %*% yw_fit$ar

se_new = sqrt( diag( dta_obs %*% yw_fit$asy.var.coef %*% t(dta_obs) ) + 
                 yw_fit$var.pred )

U_fit = loess((pred_case_rates+county_mean+1.96*se)~c(1:nrow(dta_obs)),
              span = .07)
L_fit = loess((pred_case_rates+county_mean-1.96*se)~c(1:nrow(dta_obs)),
              span = .07)

U_fit_new = loess((pred_case_rates_new+county_mean+1.96*se_new)~c(1:nrow(dta_obs)),
                  span = .07)
L_fit_new = loess((pred_case_rates_new+county_mean-1.96*se_new)~c(1:nrow(dta_obs)),
                  span = .07)

par( oma = c(2,3,1,1), mfrow=c(1,2), mar = c(0.4,2,2,0.4) )
# plot WMM
plot(1:length(case_rates),case_rates+county_mean,type="l",
     xlab="", xaxt='n', main="WMM", cex.main=2,
     ylim=c(-0.00025,0.0015))
points(1:nrow(dta_obs)+7,pred_case_rates+county_mean,pch=20,col="red",cex=1)
lines(1:nrow(dta_obs)+7,U_fit$fitted,col="blue",lty=4,lwd=1.5)
lines(1:nrow(dta_obs)+7,L_fit$fitted,col="blue",lty=4,lwd=1.5)

# plot stand alone ar
plot(1:length(case_rates),case_rates+county_mean,type="l",
     xlab="", xaxt='n', main="MLE", cex.main=2,
     ylim=c(-0.00025,0.0015))
points(1:nrow(dta_obs)+7,pred_case_rates_new+county_mean,pch=20,col="red",cex=1)
lines(1:nrow(dta_obs)+7,U_fit_new$fitted,col="blue",lty=4,lwd=1.5)
lines(1:nrow(dta_obs)+7,L_fit_new$fitted,col="blue",lty=4,lwd=1.5)

# add y label
mtext("Incidence Rate",side=2,line=.5,outer=TRUE,cex=2,las=0)

# add legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n',ylim = c(0,0.1))
legend('bottom',
       legend = c("observed", "predicted", "prediction bounds"), 
       col = c("black","red", "blue"), 
       lty=c(1,NA,4),
       pch=c(NA,20,NA),
       lwd = 2, xpd = T, horiz = TRUE, cex = 1, bty = 'n')

