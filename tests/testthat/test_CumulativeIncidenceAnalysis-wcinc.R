library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(data.table, quietly = TRUE, warn.conflicts = FALSE)

#=================================================================================
# Preparation
#=================================================================================

wcinc <- CumulativeIncidenceAnalysis$new()$.__enclos_env__$private$wcinc

# |actual - expected| <= atol + rtol * |expected|, elementwise.
expect_close <- function(actual, expected, rtol = 1e-12, atol = 1e-15) {
  expect_equal(length(actual), length(expected))
  gap <- abs(actual - expected)
  expect_true(all(gap <= atol + rtol * abs(expected)), info = sprintf("max gap %.3g", max(gap)))
}

# cmprsk's cause-1 table as EPIMIGHT has always read it: the first row at each age.
cuminc_table <- function(age, status) {
  out <- cmprsk::cuminc(ftime = age, fstatus = status, cencode = 0)$`1 1`
  data.table(age = out$time, cif = out$est, var = out$var)[, head(.SD, 1), by = age]
}

# A random cohort of integer ages and statuses 0 (censored), 1 (trait) and 2 (competing).
random_cohort <- function(i) {
  n      <- sample(c(1:12, 20, 60), 1)
  max_age <- sample(c(0, 1, 3, 8, 20), 1)
  probs  <- switch(i %% 4 + 1, c(1, 1, 1), c(1, 2, 0), c(0, 1, 1), c(2, 1, 1))
  age    <- sample(0:max_age, n, replace = TRUE)
  status <- sample(0:2, n, replace = TRUE, prob = probs)
  if (!any(status == 1)) status[sample.int(n, 1)] <- 1L
  list(age = as.integer(age), status = as.integer(status))
}

# Cohorts that reach each boundary the recursion has.
edge_cohorts <- list(
  censored_at_0       = list(age = c(0L, 0L, 1L, 2L, 2L, 3L), status = c(0L, 1L, 1L, 2L, 1L, 0L)),
  censored_past_last  = list(age = c(1L, 2L, 2L, 5L, 7L, 9L), status = c(1L, 2L, 1L, 0L, 0L, 0L)),
  cause_1_at_last_age = list(age = c(1L, 2L, 3L, 3L, 4L, 4L), status = c(1L, 0L, 2L, 1L, 1L, 0L)),
  all_cause_1         = list(age = c(1L, 1L, 2L, 3L, 3L), status = c(1L, 1L, 1L, 1L, 1L)),
  exhausted_mixed     = list(age = c(1L, 2L, 4L, 4L, 4L), status = c(1L, 0L, 1L, 2L, 1L)),
  single_event_at_0   = list(age = c(0L, 2L, 5L), status = c(1L, 0L, 0L)),
  one_person          = list(age = 4L, status = 1L)
)

# The point from survival's weighted multi-state Aalen-Johansen, read as the left limit F(age-).
survfit_left_limit <- function(age, status, weight, at) {
  event <- factor(status, 0:2, c("censored", "trait", "competing"))
  fit   <- survival::survfit(survival::Surv(age, event) ~ 1, weights = weight, id = seq_along(age),
                             conf.type = "none")
  state <- which(fit$states == "trait")
  vapply(at, function(a) {
    k <- which(fit$time < a)
    if (length(k) == 0) 0 else fit$pstate[max(k), state]
  }, numeric(1))
}

# The variance as the direct quadratic form sum_j t6_j (a_j - F b_j)^2 over every event term j
# before the row's age, with F the row's CIF `at_cif`, written per age from the definitions
# rather than cinc's v1/v2/v3.
quadratic_form_var <- function(age, status, weight, at, at_cif) {
  keep <- weight > 0
  age <- age[keep]; status <- status[keep]; weight <- weight[keep]
  terms <- list()
  surv  <- 1
  cif   <- 0
  for (t in sort(unique(age))) {
    risk <- age >= t
    y    <- sum(weight[risk])
    n    <- sum(risk)
    hit1 <- age == t & status == 1
    hit2 <- age == t & status > 1
    d1   <- sum(weight[hit1])
    d2   <- sum(weight[hit2])
    surv_after <- surv * (1 - (d1 + d2) / y)
    cif        <- cif + surv * d1 / y
    tie        <- function(m) if (m > 1) 1 - (m - 1) / (n - 1) else 1
    if (d2 > 0 && surv_after > 0) {
      terms[[length(terms) + 1]] <- c(t, surv^2 * tie(sum(hit2)) * sum(weight[hit2]^2) / y^2,
                                      cif / surv_after, 1 / surv_after)
    }
    if (d1 > 0) {
      b <- if (surv_after > 0) 1 / surv_after else 0
      terms[[length(terms) + 1]] <- c(t, surv^2 * tie(sum(hit1)) * sum(weight[hit1]^2) / y^2, 1 + b * cif, b)
    }
    surv <- surv_after
  }
  terms <- do.call(rbind, terms)
  vapply(seq_along(at), function(i) {
    f    <- at_cif[i]
    used <- terms[terms[, 1] < at[i], , drop = FALSE]
    sum(used[, 2] * (used[, 3] - f * used[, 4])^2)
  }, numeric(1))
}

# Random fractional weights on a cohort, some of them zero.
fractional <- function(n) {
  sample(c(0, 0.1, 0.25, 1 / 3, 0.5, 2 / 3, 0.9, 1), n, replace = TRUE, prob = c(1, rep(2, 7)))
}

#=================================================================================
# Tests
#=================================================================================

describe("wcinc at unit weights", {
  skip_if_not_installed("cmprsk")

  it("equals cmprsk::cuminc's first row per age, to 1e-12, on random and boundary cohorts", {
    set.seed(11)
    cohorts <- c(edge_cohorts, lapply(seq_len(300), random_cohort))

    for (cohort in cohorts) {
      ours     <- wcinc(cohort$age, cohort$status, rep(1, length(cohort$age)))
      expected <- cuminc_table(cohort$age, cohort$status)

      expect_equal(ours$age, as.integer(expected$age))
      expect_close(ours$cif, expected$cif)
      expect_close(ours$var, expected$var)
    }
  })

  it("counts cases as the people with a cause-1 event through each age", {
    set.seed(12)
    for (cohort in lapply(seq_len(50), random_cohort)) {
      ours <- wcinc(cohort$age, cohort$status, rep(1, length(cohort$age)))

      expect_equal(ours$cases, vapply(ours$age, function(a) sum(cohort$status == 1 & cohort$age <= a), numeric(1)))
    }
  })

  it("reports integer ages and se, l95 and u95 from var", {
    ours <- wcinc(edge_cohorts$censored_at_0$age, edge_cohorts$censored_at_0$status, rep(1, 6))

    expect_type(ours$age, "integer")
    expect_equal(names(ours), c("age", "cif", "se", "l95", "u95", "var", "cases"))
    expect_equal(ours$se, sqrt(ours$var))
    expect_equal(ours$l95, ours$cif - qnorm(0.975) * ours$se)
    expect_equal(ours$u95, ours$cif + qnorm(0.975) * ours$se)
  })
})

describe("wcinc on hand-checked data", {
  it("gives em4's dataset C: F(1-) = 0, F(2-) = 1/8, F(3-) = 1/4, and F(4) = 11/20 on the trailing row", {
    ours <- wcinc(c(1, 2, 2, 3, 3, 3, 4, 4), c(1, 0, 1, 1, 1, 0, 2, 0), rep(1, 8))

    expect_equal(ours$age, 0:4)
    expect_close(ours$cif, c(0, 0, 1 / 8, 1 / 4, 11 / 20))
    expect_equal(ours$var[1:2], c(0, 0))
  })

  it("gives em5's variances: 1/9 at age 2 of (1:1, 2:1, 3:0) and 5/72 at age 3 of (1:2, 2:1, 3:1, 4:0)", {
    first  <- wcinc(c(1, 2, 3), c(1, 1, 0), rep(1, 3))
    second <- wcinc(c(1, 2, 3, 4), c(2, 1, 1, 0), rep(1, 4))

    expect_close(first[age == 2, c(cif, var)], c(1 / 3, 1 / 9))
    expect_close(second[age == 3, c(cif, var)], c(1 / 4, 5 / 72))
  })
})

describe("wcinc at fractional weights", {
  it("equals survival::survfit's weighted Aalen-Johansen point at every row, to 1e-12", {
    skip_if_not_installed("survival")
    set.seed(13)

    for (i in seq_len(200)) {
      cohort <- random_cohort(i)
      weight <- fractional(length(cohort$age))
      if (!any(weight > 0 & cohort$status == 1)) next
      ours <- wcinc(cohort$age, cohort$status, weight)
      keep <- weight > 0

      expect_close(ours$cif, survfit_left_limit(cohort$age[keep], cohort$status[keep], weight[keep], ours$age))
    }
  })

  it("gives 1/144 on Codex's counterexample (events w = 0.1, 0.1 at 1; censoring w = 1 at 2)", {
    # tie_variants.out prints 0.00694444 for the people-count tie correction; Kish counts give
    # -0.0198413 and no correction 0.0138889.
    ours <- wcinc(c(1, 1, 2), c(1, 1, 0), c(0.1, 0.1, 1))

    expect_equal(ours$age, 0:2)
    expect_close(ours[age == 2, var], 1 / 144)
  })

  it("equals the direct quadratic form sum_j t6_j (a_j - F b_j)^2 on random fractional cohorts", {
    set.seed(14)
    checked <- 0

    for (i in seq_len(300)) {
      cohort <- random_cohort(i)
      weight <- fractional(length(cohort$age))
      if (!any(weight > 0 & cohort$status == 1)) next
      ours   <- wcinc(cohort$age, cohort$status, weight)

      expect_close(ours$var, quadratic_form_var(cohort$age, cohort$status, weight, ours$age, ours$cif),
                   rtol = 1e-10)
      checked <- checked + 1
    }

    expect_gt(checked, 200)
  })

  it("is never negative", {
    set.seed(15)
    for (i in seq_len(300)) {
      cohort <- random_cohort(i)
      weight <- fractional(length(cohort$age))
      ours   <- wcinc(cohort$age, cohort$status, weight)
      if (!is.null(ours)) expect_true(all(ours$var >= 0))
    }
  })

  it("leaves cif, var and cases unchanged when every weight is divided by m", {
    set.seed(16)
    for (i in seq_len(40)) {
      cohort <- random_cohort(i)
      weight <- fractional(length(cohort$age))
      if (!any(weight > 0 & cohort$status == 1)) next
      ours <- wcinc(cohort$age, cohort$status, weight)

      for (m in c(2, 7)) {
        scaled <- wcinc(cohort$age, cohort$status, weight / m)

        expect_equal(scaled$age, ours$age)
        expect_close(scaled$cif, ours$cif, rtol = 1e-12, atol = 1e-14)
        expect_close(scaled$var, ours$var, rtol = 1e-10, atol = 1e-15)
        expect_close(scaled$cases, ours$cases, rtol = 1e-12)
      }
    }
  })

  it("keeps the point and doubles cases when every row is duplicated at half its weight", {
    set.seed(17)
    for (i in seq_len(40)) {
      cohort <- random_cohort(i)
      weight <- fractional(length(cohort$age))
      if (!any(weight > 0 & cohort$status == 1)) next
      ours  <- wcinc(cohort$age, cohort$status, weight)
      twice <- wcinc(rep(cohort$age, 2), rep(cohort$status, 2), rep(weight / 2, 2))

      expect_equal(twice$age, ours$age)
      expect_close(twice$cif, ours$cif, atol = 1e-14)
      expect_close(twice$cases, 2 * ours$cases)
    }
  })

  it("treats a zero-weight row as an absent row", {
    with_zero <- wcinc(c(1, 2, 3, 9, 9), c(1, 2, 1, 1, 0), c(0.5, 1, 0.25, 0, 0))
    without   <- wcinc(c(1, 2, 3), c(1, 2, 1), c(0.5, 1, 0.25))

    expect_equal(with_zero, without)
  })
})

describe("wcinc's known gap: tied cause-1 and competing events (rwaples/epimight#8)", {
  skip_if_not_installed("cmprsk")

  it("matches cuminc's 5.470899366e-05 on the balanced 4-outcome example, not the exact 6.25e-05", {
    # n = 4000: a quarter each of cause-1 and competing at age 1, cause-1 at 2, censored at 3.
    # F(3-) = 1/2 and the exact variance is p(1 - p) / n = 6.25e-05. cinc, and wcinc with it,
    # leave out the covariance of the tied cause-1 and competing counts at age 1
    # (notes/epimight-cif-paths/tied_cov.out). A fix for #8 changes this value on purpose.
    n      <- 4000
    age    <- rep(c(1, 1, 2, 3), each = n / 4)
    status <- rep(c(1, 2, 1, 0), each = n / 4)
    ours   <- wcinc(age, status, rep(1, n))

    expect_close(ours[age == 3, cif], 0.5)
    expect_equal(ours[age == 3, var], 5.470899366e-05, tolerance = 1e-9)
    expect_close(ours[age == 3, var], tail(cmprsk::cuminc(age, status, cencode = 0)$`1 1`$var, 1))
  })
})

describe("wcinc's layout", {
  it("has an age-0 row, a row at each cause-1 age, and a trailing row when the last age has no cause-1 event", {
    ours <- wcinc(c(2, 2, 5, 7, 7, 9), c(1, 0, 1, 2, 0, 0), c(0.5, 1, 1, 0.5, 1, 0.25))

    expect_equal(ours$age, c(0L, 2L, 5L, 9L))
    expect_close(ours$cif, c(0, 0, 0.5 / 4.25, 0.5 / 4.25 + (3.75 / 4.25) / 2.75))
    expect_equal(ours$var[1:2], c(0, 0))
  })

  it("keeps the pre-jump row when the last age has a cause-1 event", {
    ours <- wcinc(c(2, 4, 4), c(1, 1, 0), c(1, 0.5, 1))

    expect_equal(ours$age, c(0L, 2L, 4L))
    expect_close(ours$cif, c(0, 0, 1 / 2.5))
  })

  it("carries the last cause-1 age's post-jump CIF to the trailing row past a competing-only age", {
    ours <- wcinc(c(1, 2, 3, 6), c(1, 1, 2, 0), c(1, 0.5, 0.5, 1))

    expect_equal(ours$age, c(0L, 1L, 2L, 6L))
    expect_close(ours$cif, c(0, 0, 1 / 3, 1 / 2))
  })
})

describe("wcinc and the sandwich's plug-in CIF", {
  sandwich <- SandwichAnalysis$new()

  it("agree at a trailing row (no cause-1 event at the last age)", {
    age    <- c(1, 2, 2, 3, 5, 5, 8)
    status <- c(1, 2, 1, 1, 0, 2, 0)
    weight <- c(0.5, 1, 0.25, 1, 0.5, 1, 1)
    ours   <- wcinc(age, status, weight)

    expect_equal(tail(ours$age, 1), 8L)
    expect_close(sandwich$calculate_cif_influence(age, status, weight, ours$age)$cif, ours$cif)
  })

  it("agree at the pre-jump row of a cause-1 event at the last age", {
    age    <- c(1, 2, 2, 3, 5, 6, 6)
    status <- c(1, 2, 1, 1, 0, 1, 0)
    weight <- c(0.5, 1, 0.25, 1, 0.5, 1, 1)
    ours   <- wcinc(age, status, weight)

    expect_equal(tail(ours$age, 1), 6L)
    expect_close(sandwich$calculate_cif_influence(age, status, weight, ours$age)$cif, ours$cif)
  })
})

describe("wcinc's return contract", {
  it("returns NULL with no event, with only competing events, and with no positively weighted cause-1 event", {
    expect_null(wcinc(c(1, 2), c(0, 0), c(1, 1)))
    expect_null(wcinc(c(1, 2, 3), c(2, 0, 2), c(1, 1, 1)))
    expect_null(wcinc(c(1, 2, 3), c(1, 2, 0), c(0, 1, 1)))
    expect_null(wcinc(integer(0), integer(0), numeric(0)))
  })
})
