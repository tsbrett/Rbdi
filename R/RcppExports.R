# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

c_choose <- function(N, K) {
    .Call('_Rbdi_c_choose', PACKAGE = 'Rbdi', N, K)
}

P_BDC <- function(x, x0, lambda, mu, t) {
    .Call('_Rbdi_P_BDC', PACKAGE = 'Rbdi', x, x0, lambda, mu, t)
}

P_BDIC <- function(x, x0, v, lambda, mu, t) {
    .Call('_Rbdi_P_BDIC', PACKAGE = 'Rbdi', x, x0, v, lambda, mu, t)
}

C_bdi_ll_linear_R0 <- function(x, delta_t, eta, gamm, dR0, R00) {
    .Call('_Rbdi_C_bdi_ll_linear_R0', PACKAGE = 'Rbdi', x, delta_t, eta, gamm, dR0, R00)
}

