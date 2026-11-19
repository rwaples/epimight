library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(data.table, quietly = TRUE, warn.conflicts = FALSE)
library(parallel, quietly = TRUE, warn.conflicts = FALSE)

source("../utils.R")

#=================================================================================
# Preparation
#=================================================================================

h2_analysis  <- HeritabilityAnalysis$new()
cif_analysis <- CumulativeIncidenceAnalysis$new()

#=================================================================================
# Tests
#=================================================================================

describe("calculate_h2", {
  it("Produces same results regardless of vectorization or not", {
    estimates <- data.table(
      id = c(1, 2, 3),
      k1 = c(0.09653048, 0.1248009, 0.1060537),
      kr = c(0.2701333, 0.3139035, 0.3562753),
      a1 = c(6, 7, 6),
      ar = c(4, 5, 4),
      rc = c(0.25, 0.25, 0.25)
    )

    results1 <- h2_analysis$calculate_h2(
      estimates$id,
      estimates$k1,
      estimates$kr,
      estimates$a1,
      estimates$ar,
      estimates$rc
    )

    results2 <- NULL

    for (i in seq_len(nrow(estimates))) {
      est      <- estimates[i, ]
      curr     <- h2_analysis$calculate_h2(est$id, est$k1, est$kr, est$a1, est$ar, est$rc)
      results2 <- rbind(results2, curr)
    }

    expect_dataframe_equal(results1, results2)
  })

  it("produces h2 > 2 when all relatives have trait if the proband have the trait", {
    tte_pop <- tte_random_probands(10000) |>
      tte_random_trait("SCZ", 0.5, 20, 11) |>
      mutate(
        relatives_n_trait = ifelse(trait_status == 1, relatives_n, 0)
      )

    tte_fh  <- tte_pop |>
      filter(relatives_n_trait > 0) |>
      select(person_id, trait_status, trait_age) |> as.data.table()

    tte_pop <- tte_pop |> select(person_id, trait_status, trait_age) |> as.data.table()

    cif_pop <- cif_analysis$run(tte = tte_pop) |>
        select(age, cif, cases) |>
        rename(pop_cif = cif, pop_cases = cases)

    cif_fh <- cif_analysis$run(tte = tte_fh) |>
        select(age, cif, cases) |>
        rename(fh_cif = cif, fh_cases = cases)

    h2 <- h2_analysis$run(
      cif = inner_join(cif_pop, cif_fh, by = join_by(age)) |>
        arrange(desc(age)) |>
        filter(row_number() == 1) |>
        as.data.table(),
      relatedness = 0.5
    )

    expect_gt(h2$h2, 2)
  })
})

describe("run", {
  h2_rows <- function(pop_cif, fh_cif, cases = 5) {
    h2_analysis$run(
      cif = data.table(age = 10L, pop_cif = pop_cif, pop_cases = cases, fh_cif = fh_cif, fh_cases = cases),
      relatedness = 0.5
    )
  }

  it("drops a row whose h2 is zero up to round-off", {
    # Both CIFs are 5/7; cmprsk and the weighted estimator round it to adjacent doubles.
    expect_equal(nrow(h2_rows(0.71428571428571419, 0.71428571428571441)), 0)
    expect_equal(nrow(h2_rows(0.71428571428571441, 0.71428571428571419)), 0)
  })

  it("keeps a small h2 that is not round-off", {
    rows <- h2_rows(0.1, 0.1 + 1e-6)
    expect_equal(nrow(rows), 1)
    expect_gt(rows$h2, 1e-6)
  })
})
