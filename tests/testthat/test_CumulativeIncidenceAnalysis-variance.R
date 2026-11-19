library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(data.table, quietly = TRUE, warn.conflicts = FALSE)

source("../utils.R")

# The two CIF paths use different variance estimators (rwaples/epimight#5). The expected
# values are hand calculations, not outputs of the code under test.

analysis <- CumulativeIncidenceAnalysis$new()

tiny_tte <- function(age, status) {
  data.table(
    person_id    = as.character(seq_along(age)),
    trait_status = as.integer(status),
    trait_age    = as.integer(age)
  )
}

describe("run's CIF variance", {
  # Ages 1, 2, 3; statuses 1, 1, 0. The row at age 2 reports F(1) = 1/3, and nobody is
  # censored before age 3, so the sampling variance of F(1) is binomial: (1/3)(2/3)/3 = 2/27.
  tte <- tiny_tte(c(1, 2, 3), c(1, 1, 0))

  it("is Greenwood's formula on the weighted path, the binomial variance here", {
    # S(1)^2 * d1 / (n (n - d1)) = (2/3)^2 * 1 / (3 * 2) = 2/27
    weighted <- analysis$run(tte = copy(tte)[, weight := 1.0])

    expect_equal(weighted[age == 2]$cif, 1 / 3)
    expect_equal(weighted[age == 2]$var, 2 / 27)
  })

  it("is cmprsk's Aalen variance on the unweighted path, larger than binomial here", {
    # cmprsk's increment variance for one cause-1 event among 3 at risk is S(0)^2 d1 / n^2
    # = 1/9. With G = F(1) / S(1) = 1/2 it enters as (1 + G)^2 - 2 G (1 + G) + G^2 = 1.
    unweighted <- analysis$run(tte = copy(tte))

    expect_equal(unweighted[age == 2]$cif, 1 / 3)
    expect_equal(unweighted[age == 2]$var, 1 / 9)
  })

  it("has no competing-event term on the weighted path", {
    # Ages 1, 2, 3, 4; statuses 2, 1, 1, 0. The row at age 3 reports F(2) = 1/4.
    # Greenwood sums cause-1 terms only: S(2)^2 * 1 / (3 * 2) = (1/2)^2 / 6 = 1/24.
    # cmprsk: 1/16 from the cause-1 event at age 2, plus F(2)^2 * (S(0)^2 / 4^2) / S(1)^2
    # = (1/16)(1/16)(16/9) = 1/144 from the competing event at age 1; 1/16 + 1/144 = 5/72.
    # Nobody is censored before age 4, so the sampling variance is (1/4)(3/4)/4 = 3/64.
    tte <- tiny_tte(c(1, 2, 3, 4), c(2, 1, 1, 0))
    weighted   <- analysis$run(tte = copy(tte)[, weight := 1.0])
    unweighted <- analysis$run(tte = copy(tte))

    expect_equal(weighted[age == 3]$cif, 1 / 4)
    expect_equal(unweighted[age == 3]$cif, 1 / 4)
    expect_equal(weighted[age == 3]$var, 1 / 24)
    expect_equal(unweighted[age == 3]$var, 5 / 72)
  })
})
