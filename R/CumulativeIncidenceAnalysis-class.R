#' @title Estimate cumulative incidence functions from competing risks data and
#' test equality across groups
#' @description
#' R6 class that performs different types of analyses.
#' @docType class
#' @import data.table
#' @import dplyr
#' @import dtplyr
#' @import parallel
#' @import tidyr
#' @export
CumulativeIncidenceAnalysis <- R6::R6Class( #nolint
  "CumulativeIncidenceAnalysis",
  inherit = Analysis,
  private = list(
    args = NULL,
    #' @description
    #' Dynamically builds a dplyr filter expression for the given columns.
    #'
    #' @param columns List of column names to filter.
    #' @param row Data.table with the values to use for each column.
    #' @returns Filter expression.
    build_filter_expression = function(columns, row) {
      if (is.null(ncol(row))) {
        col <- columns[[1]]

        return(
          rlang::call2("==", sym(col), row)
        )
      }

      results <- NULL

      for (col in columns) {
        eq_expr <- rlang::call2("==", sym(col), row[[col]])

        if (is.null(results)) {
          results <- eq_expr
        } else {
          results <- rlang::call2("&", results, eq_expr)
        }
      }

      results
    },
    #' @description
    #' Dynamically builds a dplyr mutate expression for the given columns.
    #'
    #' @param columns List of column names to filter.
    #' @param row Data.table with the values to use for each column.
    #' @returns Filter expression.
    build_mutate_expression = function(columns, row) {
      results <- list()

      if (is.null(ncol(row))) {
        col <- columns[[1]]

        results[[col]] <- row
      } else {

        for (col in columns) {
          results[[col]] <- row[[col]]
        }
      }

      rlang::call2("mutate", as.symbol("group_results"), !!!results)
    },
    #' @description
    #' Cause-1 cumulative incidence of one cohort: cmprsk's `cinc` recursion (Aalen-Johansen point,
    #' Aalen's delta-method variance) run on weighted counts. At unit weights it equals
    #' `cmprsk::cuminc` with the first row kept at each age.
    #'
    #' Rows are cmprsk's: age 0, each cause-1 age with the left limit F(age-) and its variance,
    #' and the last follow-up age when it has no cause-1 event, where F(age-) = F(age). Each
    #' variance term uses the weights' sum of squares over the squared weighted risk set, with
    #' cmprsk's tie correction 1 - (m - 1) / (n - 1) on people counts. Like cmprsk, it leaves out
    #' the covariance of cause-1 and competing events tied at one age (rwaples/epimight#8).
    #' `cases` counts cause-1 weight through each age, scaled by sum(w) / sum(w^2).
    #'
    #' @param trait_age Integer age at event or censoring.
    #' @param trait_status 0 censored, 1 the trait, any other value a competing event.
    #' @param weight Row weight, at least 0; a zero-weight row is the same as an absent row.
    #' @returns A data.table with `age`, `cif`, `se`, `l95`, `u95`, `var` and `cases`, or NULL
    #'   when no cause-1 event has a positive weight.
    wcinc = function(trait_age, trait_status, weight) {
      keep <- weight > 0
      if (!any(keep & trait_status == 1)) return(NULL)

      # Plain vectors and one rowsum: per-stratum data.tables cost more memory than cmprsk did.
      age    <- as.integer(trait_age[keep])
      status <- trait_status[keep]
      w      <- weight[keep]
      is_1   <- status == 1
      is_2   <- status != 0 & !is_1
      sums   <- rowsum(cbind(w, 1, w * is_1, w^2 * is_1, is_1, w * is_2, w^2 * is_2, is_2), age, reorder = TRUE)
      ages   <- as.integer(rownames(sums))
      at_risk <- rev(cumsum(rev(sums[, 1])))
      n_risk  <- rev(cumsum(rev(sums[, 2])))
      d1 <- sums[, 3]
      d2 <- sums[, 6]
      cases <- cumsum(d1) * sum(w) / sum(w^2)

      tie <- function(m, n) ifelse(m > 1, 1 - (m - 1) / (n - 1), 1)
      ev  <- which(d1 + d2 > 0)
      y   <- at_risk[ev]
      surv        <- cumprod(1 - (d1[ev] + d2[ev]) / y)
      surv_before <- c(1, head(surv, -1))
      cif  <- cumsum(surv_before * d1[ev] / y)
      t6_2 <- ifelse(d2[ev] > 0 & surv > 0, surv_before^2 * tie(sums[ev, 8], n_risk[ev]) * sums[ev, 7] / y^2, 0)
      t6_1 <- surv_before^2 * tie(sums[ev, 5], n_risk[ev]) * sums[ev, 4] / y^2
      inv  <- ifelse(surv > 0, 1 / surv, 0)
      a2   <- cif * inv
      a1   <- 1 + inv * cif
      v1   <- cumsum(t6_2 * a2^2 + t6_1 * a1^2)
      v2   <- cumsum((t6_2 * a2 + t6_1 * a1) * inv)
      v3   <- cumsum((t6_2 + t6_1) * inv^2)
      jump      <- d1[ev] > 0
      jump_age  <- ages[ev][jump]
      jump_cif  <- cif[jump]
      jump_var  <- pmax(v1[jump] + jump_cif^2 * v3[jump] - 2 * jump_cif * v2[jump], 0)

      # cmprsk's layout: age 0, the left limit at each cause-1 age, then the last follow-up age,
      # whose row holds the final CIF unless that age is itself a cause-1 age.
      out_age <- sort(unique(c(0L, jump_age, ages[length(ages)])))
      i       <- match(out_age, jump_age)
      at_jump <- !is.na(i)
      out_cif <- ifelse(at_jump, c(0, head(jump_cif, -1))[i], ifelse(out_age == 0L, 0, tail(jump_cif, 1)))
      out_var <- ifelse(at_jump, c(0, head(jump_var, -1))[i], ifelse(out_age == 0L, 0, tail(jump_var, 1)))
      index   <- findInterval(out_age, ages)

      data.table(
        age   = out_age,
        cif   = out_cif,
        se    = sqrt(out_var),
        l95   = out_cif - qnorm(0.975) * sqrt(out_var),
        u95   = out_cif + qnorm(0.975) * sqrt(out_var),
        var   = out_var,
        cases = ifelse(index > 0, cases[pmax(index, 1L)], 0)
      )
    },
    #' @description
    #' Runs CIF on the given TTE data as a single group. A cohort with no `weight` column has
    #' weight 1 throughout.
    #'
    #' @param tte Data.table of TTE data to use.
    #' @returns Risk estimations.
    run_single = function(tte) {
      if (!is.data.table(tte)) {
        stop("Given TTE was not a data.table")
      }

      weight <- if ("weight" %in% colnames(tte)) tte$weight else rep(1, nrow(tte))

      private$wcinc(tte$trait_age, tte$trait_status, weight)
    }
  ),
  public = list(
    initialize = function() {
    },
    #' @description
    #' Runs CIF on the given TTE data using the provided arguments.
    #'
    #' @returns Risk estimations.
    run = function(...) {
      validator <- ArgumentsValidator$new(
        tte = list(
          type     = "data.table",
          required = TRUE,
          columns  = list(
            person_id = list(
              type     = "string",
              required = TRUE
            ),
            trait_status = list(
              type     = "integer",
              required = TRUE,
              minimum  = 0
            ),
            trait_age = list(
              type     = "integer",
              required = TRUE,
              minimum  = 0
            ),
            weight = list(
              type    = "numeric",
              minimum = 0,
              maximum = 1
            )
          )
        ),
        stratify_columns = list(
          type  = "list",
          items = list(type = "string")
        )
      )

      args <- validator$run(...)

      self$assert_unique_individuals_tte(args$tte)

      if (!exists("stratify_columns", where = args) || length(args$stratify_columns) == 0) {
        return(
          private$run_single(args$tte)
        )
      }

      stratify_symbols <- rlang::syms(args$stratify_columns)
      permutations  <- args$tte |>
        select(!!!stratify_symbols) |>
        distinct(!!!stratify_symbols) |>
        arrange(!!!stratify_symbols) |>
        as.data.frame()

      runner <- function(idx) {
        filter_expr <- private$build_filter_expression(
          args$stratify_columns,
          permutations[idx, ]
        )

        group_tte <- args$tte[rlang::eval_tidy(filter_expr)]

        if (nrow(group_tte) == 0) return(NULL)

        group_results <- private$run_single(group_tte)

        if (is.null(group_results)) return(NULL)

        mutate_expr <- private$build_mutate_expression(
          args$stratify_columns, permutations[idx, ]
        )

        rlang::eval_tidy(mutate_expr) |> collect()
      }

      indexes       <- seq_len(nrow(permutations))
      group_results <- mclapply(indexes, runner)

      rbindlist(group_results)
    }
  )
)
