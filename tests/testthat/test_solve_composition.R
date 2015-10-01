library(breedTools)
context("Testing each of 2 main functions (in different ways) with a subset of Yorkshire data")


load("dat_for_breedTools_test.RData")




test <- solve_composition(Y8K_sample, X8K)

test_that("solve_composition returns proper class", {
  expect_is(test, "matrix")
})

test_that("solve_composition computes expected column means", {
  expect_equal(unname(round(colMeans(test)[1], 5)), 0.01225)
  expect_equal(unname(round(colMeans(test)[2], 5)), 0.01745)
  expect_equal(unname(round(colMeans(test)[3], 5)), 0.03085)
  expect_equal(unname(round(colMeans(test)[4], 5)), 0.93945)
  expect_equal(unname(round(colMeans(test)[5], 5)), 0.39806)
})




test_w_ped <- solve_composition(Y8K_sample, X8K,
                                ped = yorkshirePed)

test_that("sim_composition with ped option returns proper class", {
  expect_is(test_w_ped, "matrix")
})

test_that("solve_composition with ped option computes expected column means", {
  expect_equal(unname(round(colMeans(test_w_ped)[1], 5)), 0.04142)
  expect_equal(unname(round(colMeans(test_w_ped)[2], 5)), 0.00399)
  expect_equal(unname(round(colMeans(test_w_ped)[3], 5)), 0.02582)
  expect_equal(unname(round(colMeans(test_w_ped)[4], 5)), 0.92878)
  expect_equal(unname(round(colMeans(test_w_ped)[5], 5)), 0.27156)
})


test_w_ped_mia <- solve_composition(Y8K_sample, X8K,
                                    ped = yorkshirePed,
                                    mia = TRUE)

test_that("mia argument gives expected allele codings and output", {
  expect_equal(ncol(test_w_ped_mia), ncol(Y8K_sample))
  expect_true(names(table(test_w_ped_mia))[1] %in% c("?", "0", "0.5", "1"))
  expect_true(names(table(test_w_ped_mia))[2] %in% c("?", "0", "0.5", "1"))
  expect_true(names(table(test_w_ped_mia))[3] %in% c("?", "0", "0.5", "1"))
  expect_true(names(table(test_w_ped_mia))[4] %in% c("?", "0", "0.5", "1"))
})




test_w_group <- solve_composition(Y8K_sample, X8K, 
                                  groups = list("Yorkshire" = rownames(Y8K_sample)[1:50],
                                                "Landrace" = rownames(Y8K_sample)[51:100]))

test_that("sim_composition with group option returns proper class", {
  expect_is(test_w_group, "list")
})









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
            