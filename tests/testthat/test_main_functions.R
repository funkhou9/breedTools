library(breedTools)
context("Testing each of 2 main functions (in different ways) with a subset of Yorkshire data")


load("dat_for_breedTools_test.RData")




test <- solve_composition(Y8K_sample, X8K, 
                          names = c("Duroc", "Hampshire", "Landrace", "Yorkshire"))

test_that("solve_composition returns proper class", {
  expect_is(test, "matrix")
})

test_that("solve_composition computes expected column means", {
  expect_equal(unname(round(colMeans(test)["Duroc"], 5)), 0.01225)
  expect_equal(unname(round(colMeans(test)["Hampshire"], 5)), 0.01745)
  expect_equal(unname(round(colMeans(test)["Landrace"], 5)), 0.03085)
  expect_equal(unname(round(colMeans(test)["Yorkshire"], 5)), 0.93945)
  expect_equal(unname(round(colMeans(test)["R2"], 5)), 0.39806)
})




test_w_ped <- solve_composition(Y8K_sample, X8K,
                                names = c("Duroc", "Hampshire", "Landrace", "Yorkshire"),
                                ped = yorkshirePed)

test_that("sim_composition with ped option returns proper class", {
  expect_is(test_w_ped, "matrix")
})

test_that("solve_composition with ped option computes expected column means", {
  expect_equal(unname(round(colMeans(test_w_ped)["Duroc"], 5)), 0.04307)
  expect_equal(unname(round(colMeans(test_w_ped)["Hampshire"], 5)), 0.00416)
  expect_equal(unname(round(colMeans(test_w_ped)["Landrace"], 5)), 0.03003)
  expect_equal(unname(round(colMeans(test_w_ped)["Yorkshire"], 5)), 0.92274)
  expect_equal(unname(round(colMeans(test_w_ped)["R2"], 5)), 0.2677)
})




test_w_group <- solve_composition(Y8K_sample, X8K, 
                          names = c("Duroc", "Hampshire", "Landrace", "Yorkshire"),
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
            