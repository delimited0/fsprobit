# functions for multiplying multinomial probit constraint matrices

#' choice_idx is in original problem space, before identification
#' @param A n_choices-1 x n_choices-1 identified constraint matrix
#' @param v n_choices-1 x p vector
#' @return v relative to the choice index
vec_relative_choice = function(choice_idx, v) 
{
    vec = v[choice_idx] - v
    vec[choice_idx] = v[choice_idx]
    return(matrix(vec))
}

#' choice_idx is in original problem space, before identification
rows_relative_choice = function(choice_idx, M)
{
    row_vec = M[choice_idx, , drop = FALSE]
    mat = matrix(row_vec, nrow = nrow(M), ncol = ncol(M), byrow = TRUE) - M
    mat[choice_idx, ] = M[choice_idx, ]
    return(mat)
}

cols_relative_choice = function(choice_idx, M)
{
    col_vec = M[, choice_idx, drop = FALSE]
    mat = matrix(col_vec, nrow = nrow(M), ncol = ncol(M), byrow = FALSE) - M
    mat[, choice_idx] = M[, choice_idx]
    return(mat)
}

sandwich_choice = function(choice_idx, M)
{
    M_rel = rows_relative_choice(choice_idx, M)
    return(cols_relative_choice(choice_idx, M_rel))
}

