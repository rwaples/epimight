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

      rows   <- data.table(age = as.integer(trait_age[keep]), status = trait_status[keep], w = weight[keep])
      by_age <- rows[, .(
        w_all = sum(w),
        n_all = .N,
        d1    = sum(w[status == 1]),
        q1    = sum(w[status == 1]^2),
        m1    = sum(status == 1),
        d2    = sum(w[status != 0 & status != 1]),
        q2    = sum(w[status != 0 & status != 1]^2),
        m2    = sum(status != 0 & status != 1)
      ), keyby = age]
      by_age[, `:=`(at_risk = rev(cumsum(rev(w_all))), n_risk = rev(cumsum(rev(n_all))))]
      by_age[, cases := cumsum(d1) * sum(rows$w) / sum(rows$w^2)]

      tie <- function(m, n) ifelse(m > 1, 1 - (m - 1) / (n - 1), 1)
      ev  <- by_age[d1 + d2 > 0]
      ev[, surv := cumprod(1 - (d1 + d2) / at_risk)]
      ev[, surv_before := c(1, head(surv, -1))]
      ev[, cif := cumsum(surv_before * d1 / at_risk)]
      ev[, `:=`(
        t6_2 = ifelse(d2 > 0 & surv > 0, surv_before^2 * tie(m2, n_risk) * q2 / at_risk^2, 0),
        t6_1 = surv_before^2 * tie(m1, n_risk) * q1 / at_risk^2,
        inv  = ifelse(surv > 0, 1 / surv, 0)
      )]
      ev[, `:=`(a2 = cif * inv, a1 = 1 + inv * cif)]
      ev[, `:=`(
        v1 = cumsum(t6_2 * a2^2 + t6_1 * a1^2),
        v2 = cumsum((t6_2 * a2 + t6_1 * a1) * inv),
        v3 = cumsum((t6_2 + t6_1) * inv^2)
      )]
      jumps <- ev[d1 > 0, .(age, cif, var = pmax(v1 + cif^2 * v3 - 2 * cif * v2, 0))]

      # cmprsk's layout: age 0, the left limit at each cause-1 age, then the last follow-up age.
      out <- rbind(
        data.table(age = 0L, cif = 0, var = 0),
        data.table(age = jumps$age, cif = c(0, head(jumps$cif, -1)), var = c(0, head(jumps$var, -1))),
        data.table(age = by_age[.N, age], cif = tail(jumps$cif, 1), var = tail(jumps$var, 1))
      )[, head(.SD, 1), by = age]
      index <- findInterval(out$age, by_age$age)

      out[, .(
        age, cif,
        se    = sqrt(var),
        l95   = cif - qnorm(0.975) * sqrt(var),
        u95   = cif + qnorm(0.975) * sqrt(var),
        var,
        cases = ifelse(index > 0, by_age$cases[pmax(index, 1L)], 0)
      )]
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
