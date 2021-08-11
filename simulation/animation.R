#f_ppclp(train.Data, s_K, s_q, lambda_p, draw = TRUE);
f_Animeppslp(m_Y, s_K, s_q, pen=1,
             increment=0.1, maxit=60, s_lambdap=2,
             init_w = 6, init_Theta=NULL, tol=1e-3)
s_n <- 100L;
s_K <- 25L;
s_q <- 2L;
s_lambda_L = 0.1
noise <- runif(s_n, 0, 2*pi);
m_Y <- matrix(5*c(sin(noise), cos(noise)), s_n, s_q, byrow=FALSE)


s_lambdas = CV_ppclp(m_Y, k=5, range=3, s_K, s_q,
                     pen=1, increment=0.1, maxit=50);

m_Theta = f_ppclp(m_Y, s_K, s_q, s_lambdas,
                   pen=1, increment=0.1, maxit=50,
                   init_w = 6, init_Theta=NULL,
                   tol=1e-3, draw = TRUE)

f_plot_ppclp(m_Y, m_Theta, s_lambda_L, s_lambdas, maxit)
