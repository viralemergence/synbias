#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::cube ferrari(const arma::mat& initial_pop,
                   const arma::mat& interactions,
                   double beta_s,
                   int timesteps) {

  // Set up population list
  arma::cube state(initial_pop.n_rows, initial_pop.n_cols, timesteps);
  state.slice(0) = initial_pop;

  arma::mat current_pop = state.slice(1);

  // Print dimensions of variables after initialization
  if (current_pop.n_elem != initial_pop.n_elem) {
    Rcpp::stop("Error: Dimensions of current_pop and initial_pop do not match.");
  }
  arma::vec col_sums(initial_pop.n_cols);
  arma::vec total_strains(initial_pop.n_rows);
  arma::mat next_pop(initial_pop.n_rows, initial_pop.n_cols);
  arma::mat int_indiv;

  // Loop through timesteps
  for (int i = 0; i < (timesteps - 1); ++i) {
    current_pop = state.slice(i);
    col_sums.zeros();  // Ensure col_sums is initialized to zeros
    for (int col = 0; col < current_pop.n_cols; ++col) {
      col_sums(col) = arma::accu(current_pop.col(col));
    }

    if (arma::accu(col_sums)==current_pop.n_elem) {
      state.slice(i + 1) = current_pop;
    } else {
      total_strains = sum(current_pop, 1);

      for (int j = 0; j < current_pop.n_rows; ++j) {
        arma::rowvec indiv_prev = current_pop.row(j);
        // Convert indiv_prev to a boolean vector
        arma::uvec indices = find(indiv_prev > 0);
        if (indices.n_elem == 1) {
          // Only one TRUE value, select one row
          int_indiv = interactions.submat(indices(0), 0, indices(0), interactions.n_cols - 1);
        } else {
          // Multiple TRUE values, select multiple rows
          int_indiv = interactions.rows(indices);
        }

        arma::vec indiv_new(current_pop.n_cols);

        for (int k = 0; k < indiv_new.size(); ++k) {
          if (total_strains[j] > 1) {
            double rbinom_result = R::rbinom(col_sums[k], arma::prod(int_indiv.col(k)) * beta_s);
            indiv_new[k] = std::max(indiv_prev[k], std::min(rbinom_result, 1.0));
          } else if (total_strains[j] == 1) {
            double rbinom_result = R::rbinom(col_sums[k], int_indiv(0, k) * beta_s);
            indiv_new[k] = std::max(indiv_prev[k], std::min(rbinom_result, 1.0));
          } else {
            double rbinom_result = R::rbinom(col_sums[k], beta_s);
            indiv_new[k] = std::max(indiv_prev[k], std::min(rbinom_result, 1.0));
          }
        }

        next_pop.row(j) = arma::trans(indiv_new);
      }

      state.slice(i + 1) = next_pop;
    }
  }

  return state;
}
