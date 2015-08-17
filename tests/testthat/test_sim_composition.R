library(breedTools)
context("Testing the output of sim_composition()")


load("dat_for_breedTools_test.RData")



test_that("sim_composition stops since only one breed is provided", {
  expect_error(sim_composition(Y8K_sample, X8K, rep = 5, par1 = list("Yorkshire" = rownames(Y8K_sample))),
               "There must be at least two breeds to choose from")
})


test_that("sim_composition returns proper class", {
  expect_is(sim_composition(Y8K_sample, X8K, rep = 5,
                            par1 = list("Yorkshire" = rownames(Y8K_sample)[1:50],
                                        "Landrace" = rownames(Y8K_sample)[51:100])),
            "matrix")
})


