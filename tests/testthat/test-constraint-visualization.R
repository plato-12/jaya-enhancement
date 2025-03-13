library(testthat)

test_that("Enhanced jaya_multi_enhanced function tracks constraints correctly", {
  skip_if_not_installed("Jaya")

  # Define test objectives
  objective1 <- function(x) sum(x^2)
  objective2 <- function(x) sum((x - 2)^2)

  # Define constraints
  constraint1 <- function(x) sum(x) - 3  # Sum of variables must be <= 3
  constraint2 <- function(x) -sum(x) + 1 # Sum of variables must be >= 1

  # Test with a small problem to keep the test quick
  set.seed(123)
  result <- jaya_multi_enhanced(
    objectives = list(objective1, objective2),
    lower = c(-2, -2),
    upper = c(2, 2),
    popSize = 10,
    maxiter = 10,
    n_var = 2,
    constraints = list(constraint1, constraint2)
  )

  # Check that the result is a jaya_multi object
  expect_s3_class(result, "jaya_multi")

  # Check that constraint info is stored
  expect_false(is.null(attr(result, "constraint_info")))

  # Check that constraint violations are tracked
  constraint_info <- attr(result, "constraint_info")
  expect_true("violations" %in% names(constraint_info))
  expect_true("feasible" %in% names(constraint_info))

  # Check dimensions of constraint violations
  expect_equal(nrow(constraint_info$violations), nrow(result$Solutions))
  expect_equal(ncol(constraint_info$violations), 2) # Two constraints

  # Check that feasibility flags match constraint violations
  for (i in 1:nrow(result$Solutions)) {
    x <- as.numeric(result$Solutions[i, 1:2])
    expect_equal(
      constraint_info$feasible[i],
      all(constraint_info$violations[i, ] <= 0)
    )
  }
})

test_that("plot_jaya_constraints function works correctly", {
  skip_if_not_installed("Jaya")
  skip_if_not_installed("ggplot2")

  # Define test objectives
  objective1 <- function(x) sum(x^2)
  objective2 <- function(x) sum((x - 2)^2)

  # Define constraints
  constraint1 <- function(x) sum(x) - 3  # Sum of variables must be <= 3

  # Create a test result
  set.seed(456)
  result <- jaya_multi_enhanced(
    objectives = list(objective1, objective2),
    lower = c(-2, -2),
    upper = c(2, 2),
    popSize = 10,
    maxiter = 5, # Keep small for test speed
    n_var = 2,
    constraints = list(constraint1)
  )

  # Test basic plotting functionality
  plot_obj <- plot_jaya_constraints(result)
  expect_s3_class(plot_obj, "ggplot")

  # Test with specific constraint indices
  plot_obj2 <- plot_jaya_constraints(result, constraints_to_show = 1)
  expect_s3_class(plot_obj2, "ggplot")

  # Test with invalid constraint indices
  expect_error(
    plot_jaya_constraints(result, constraints_to_show = 2),
    "Invalid constraint indices"
  )

  # Test with invalid objective indices
  expect_error(
    plot_jaya_constraints(result, objective_indices = c(1, 3)),
    "Invalid objective indices"
  )
})

test_that("plot_jaya_constraints handles non-enhanced objects gracefully", {
  skip_if_not_installed("Jaya")
  skip_if_not_installed("ggplot2")

  # Mock a regular jaya_multi object without constraint info
  mock_result <- list(
    Pareto_Front = data.frame(
      x1 = c(1, 2, 3),
      x2 = c(4, 5, 6),
      Objective_1 = c(1, 2, 3),
      Objective_2 = c(3, 2, 1)
    ),
    Solutions = data.frame(
      x1 = c(1, 2, 3, 4),
      x2 = c(4, 5, 6, 7),
      Objective_1 = c(1, 2, 3, 4),
      Objective_2 = c(3, 2, 1, 0)
    )
  )
  class(mock_result) <- "jaya_multi"

  # Function should warn but still produce a basic plot
  expect_warning(
    plot_obj <- plot_jaya_constraints(mock_result),
    "No constraint information found"
  )
  expect_s3_class(plot_obj, "ggplot")
})
