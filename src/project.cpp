#include <Rcpp.h>
using namespace Rcpp;

//' Projection
//'
//' Compute scalar p_i for each y_i.
//' @param y data matrix y
//' @param theta weight matrix theta
//' @param m_Beta transform weight matrix m_Beta
//' @return A structure represent the projected points, including
//'         \code{m_B} the projected spline matrix
//'         \item{v_pi}{The projection index for each point}
//'         \item{dist}{The total squared distance from the curve}
//'         \item{dist_ind}{The squared distances from the curve to each of the respective points}
//'
//' @keywords Projection, Spline
//' @noRd
//' @examples
//' y <- matrix(c(13, -4, 2, -4, 11, -2, 2, -2, 8, 10), 5, 2, byrow = TRUE)
//' theta <- matrix(c(4.92, -2.67, 5.83, -1.76, 6.74, -0.85,
//'                   7.65, 0.054, 8.56, 0.96, 9.47, 1.87), 6, 2, byrow = TRUE)
//' m_Beta <- diff(theta)
//' out <- f_projection(y, theta, m_Beta)
//' proj_p <- out$m_B %*% theta
//' plot(y[, 1], y[, 2], type = "o")
// [[Rcpp::export]]
List project(
    const NumericMatrix& m_y,
    const NumericMatrix& m_theta,
    const NumericMatrix& m_beta
) {
    int s_n = m_y.nrow();
    int s_k = m_theta.nrow();
    int s_q = m_y.ncol();
    int s_seq = s_k - 1;

      // argument checks
    if (m_theta.ncol() != s_q) {
        stop("'y' and 'theta' must have an equal number of columns");
    }
    if (s_k < 2) {
        stop("'theta' must contain at least two rows.");
    }
    if (s_n == 0) {
        stop("'y' must contain at least one row.");
    }

    // preallocate variables
    double temp, u;

    
    //Precompute Values
    
    ////Define var
    double s_h = 1 / s_seq;
    NumericVector v_tk = no_init(s_seq);
    NumericVector v_denorm = no_init(s_seq);

    ////Initialize t_k = (0, 1, ..., s_k - 2) * s_h
    for (int k = 0; k < s_seq; ++k) {
        v_tk[k] = k * s_h;
    }

    ////Initialize denominator
    for (int k = 0; k < s_seq; ++k) {
        temp = 0;
        // OPTIMIZATION:
        //  v_denorm[k] = sum(pow(m_beta(k, _), 2));
        //  if (v_denorm[k] == 0){ v_denorm[k] = 1;}
        for (int j = 0; j < s_q; ++j) {
            temp += m_beta[j * s_seq + k] * m_beta[j * s_seq + k];
        }
        if (temp == 0) {
            temp = 1;
        }
        v_denorm[k] = temp;
    }
    
    //Allocate output data structure
    NumericVector v_proj = no_init(s_n);    // vector of projected values
    NumericVector v_dist = no_init(s_n);    // vector of distance from projected values
    NumericMatrix m_b = no_init(s_n, s_k);  // matrix B
    std::fill(m_b.begin(), m_b.end(), 0);

    //Allocate intermidiate data structure
    NumericVector y_i = no_init(s_q);
    NumericVector v_test = no_init(s_q);
    double a_k;
    double thresh;
    double s_pi;
    double dist;
    int seg;


    // outer loop (observations)
    for (int i = 0; i < s_n; ++i) {
        // reset best distance, best pi, best seg number
        s_pi = -1;
        a_k = -1;
        dist = R_PosInf;
        seg = -1;

        // copy current point to y_i
        for (int j = 0; j < s_q; ++j) {
            y_i[j] = m_y(i, j);
        }

        // inner loop (segments)
        for (int k = 0; k < s_seq; ++k) {
            // OPTIMIZATION:
            //  NumericVector diff = y_i - m_theta(k, _);
            //  double t = sum(m_beta(k, _) * diff) / v_denorm[k];
            thresh = 0;
            for (int j = 0; j < s_q; ++j) {
                thresh +=  m_beta(k, j) * (y_i[j] - m_theta(k, j));
            }

            thresh = thresh / v_denorm[k];

            if (thresh < 0){
                thresh = 0;
            } else if (thresh > 1) {
                thresh = 1;
            } 

            
            // Calculate the projected point on segment k
            // and the distance for this piece
            //  OPTIMIZATION:
            //      NumericVector test = m_theta(k, _) + thresh * m_beta(k, _);
            //      double dist = sum(pow(test - y_i, 2));
            temp = 0;
            for (int j = 0; j < s_q; ++j) {
                u = m_theta(k, j) + thresh * m_beta(k, j);
                temp += (u - y_i[j]) * (u - y_i[j]);
            }

            // if this is better than previous found, store it
            if (temp < dist) {
                dist = temp;
                a_k = thresh;
                s_pi = s_h * a_k + v_tk[k];
                seg = k;
            }
        }
        
        // save the best projection the output data structures
        v_proj[i] = s_pi;
        v_dist[i] = dist;
        m_b(i, seg) = 1 - a_k;
        m_b(i, seg + 1) = a_k;
    }

    dist = sum(v_dist);
    // return output
    List out;
    out["m_b"] = m_b;
//   out["v_dist"] = v_dist;
//    out["v_pi"] = v_proj;
    out["dist"] = dist;
    
    return out;
}


