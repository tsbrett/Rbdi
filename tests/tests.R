# For namespace
#useDynLib(Rbdi)
#importFrom(Rcpp, sourceCpp)

Rcpp::sourceCpp("./src/bdi.cpp")

ll_test = bdi_likelihood(z=c(40,20), delta_t = 0.1, eta = 1, gamm = 1,
                         dR0 = c(0.0,0.01), R00 = 0.12*1)


bdi_likelihood_ratio_test(ll_test)

ts = read.csv("./tests/test.csv", row.names = 1)
ts = unlist(ts, use.names=FALSE)

test1 = bdi_likelihood(z=ts, delta_t = 6, eta = 1./7, gamm = 1./7,
                      dR0 = seq(0.0,0.0004,0.0004/20), R00 = seq(0,1,1.0/20))


test2 = bdi_lrt_moving_window(z=ts, delta_t = 6, eta = 1./7, gamm = 1./7,
                              dR0 = seq(0.0,0.0004,0.0004/20),
                              R00 = seq(0,1,1.0/20), window=length(ts))

plot(test2$i, test2$cox_delta)

ll2 = na.omit(ll_test)
ll2[ll2$LogLike == stats$ML_e, ]









P_BDIC(x=40,x0=20,v=1, lambda= 0.9, mu=1, t=10)
P_BDC(x=40,x0=20, lambda= 0.9, mu=1, t=10)
P_BDC(x=30,x0=20, lambda= 0.9, mu=1, t=10)
P_BDC(x=30,x0=20, lambda=0.9, mu=1, t=10)


norm = 0.0
df = data.frame(i=integer(),p=double(),norm=double())
for(i in seq(0,50,1)){
  p = P_BDC(x=i,x0=20, lambda= 0.9, mu=1, t=10)
  norm = norm + p
  df = rbind(df, list(i=i,p=p,norm=norm))
}
plot(df[,"i"],df[,"p"])






x=30
x0=20
P_BDC(x=x,x0=x0, lambda= 0.12, mu=1, t=0.1)
P_BDIC(x=x,x0=x0,v=1, lambda= 1.1, mu=1, t=10)
test = C_bdi_ll_linear_R0(x=c(x,x0), delta_t = 0.1, eta = 1, gamm = 1,
                          dR0 = 0.0, R00 = 0.12*1)
exp(test)

ll_test = bdi_likelihood(z=c(x,x0), delta_t = 0.1, eta = 1, gamm = 1,
                         dR0 = c(0.0,0.01), R00 = 0.12*1)


bdi_likelihood_ratio_test(ll_test)


exp(ll_test$LogLike)


norm = 0.0
for(i in 0:50){
  p = P_BDIC(x=i,x0=20,v=1, lambda= 0.9, mu=1, t=10)
  norm = norm + p
  print(c(i,p,norm))
}


print(c(p,norm))

#WTF is going on here

####

pow <- function(a,b){
  #if !((0 == 0) & (1 ==0))
  return(a^b)
}

P_BDIR <- function(x, x0, v,lambda, mu, t){

  A_BD = mu*(exp((lambda-mu)*t) -1)/(lambda*exp((lambda-mu)*t) - mu); #//tick
  B_BD = lambda*(exp((lambda-mu)*t) -1)/(lambda*exp((lambda-mu)*t) - mu);# //tick
  m_BDI = (v/(lambda-mu))*(exp((lambda-mu)*t) -1.0); #//tick

  r_BDI = v/lambda;
  p_BDI = m_BDI/(m_BDI +r_BDI);

  p_BDI = m_BDI/(m_BDI+v); #// this is not quite 1 - p
  R0 = lambda/mu; #// yes but unused
  rho = A_BD/mu;

  prob_BD = 0.0;
  prob_BDI0 = 0.0;
  prob_BDI = 0.0;
  min_x = 0;

  #for(i=0; i <= x; i++){
  for(i in 0:x){

    # //Calculate prob_BD -- contribution to prob_BDI from initial infected
    #// This is wront
    prob_BD = 0.0;
    if(x-i >0 && x0 >0) #// tick
    {
      #//Do I need to round?
        min_x = min(x-i,x0);
      #for (j = 0; j <= min_x; j++ ){
      for(j in 0:min_x){
        prob_BD = prob_BD + choose(x0,j)*choose(x0 + (x-i) -j - 1, x0 -1)* #//tick
        pow(A_BD,(x0-j))* #//tick
        pow(B_BD,((x-i)-j))* #//tick
        pow((1.0 - A_BD-B_BD),j); #//tick
      }
    }
    else prob_BD = pow(A_BD,x0); #// tick (no risk of overcounting as x-i <= 0 iff x==i)


    #//Calculate prob_BDI0 -- contribution to prob_BDI from new importations,
    #// given by a negative binomial distribution
    print(prob_BD)
    prob_BDI0 = choose(i+r_BDI-1,i)*pow(p_BDI,i)*pow((1-p_BDI),r_BDI); #// tick
    print(prob_BDI0)
    prob_BDI = prob_BDI + prob_BD*prob_BDI0; #//tick
  }
  return(prob_BDI)
}





P_BDIC(x=1,x0=0, v=1, lambda= 0.12, mu=1, t=0.1)
P_BDIR(x=1,x0=0, v=1, lambda= 0.12, mu=1, t=0.1)



