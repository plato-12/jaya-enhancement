#' Enhanced Multi-Objective Jaya Algorithm with Constraint Tracking
#'
#' This function extends the original jaya_multi function to track and store
#' information about constraint violations throughout the optimization process.
#'
#' @param objectives A list of objective functions to optimize.
#' @param lower Numeric vector specifying the lower bounds for variables.
#' @param upper Numeric vector specifying the upper bounds for variables.
#' @param popSize Population size. Default is 50.
#' @param maxiter Maximum number of iterations.
#' @param n_var Number of variables.
#' @param seed Random seed for reproducibility. Default is `NULL`.
#' @param suggestions Data frame of initial suggestions for starting population. Default is an empty data frame.
#' @param constraints A list of constraint functions. Each constraint should return a non-positive value if satisfied.
#' @param adaptive_pop Logical. Whether to adapt population size during optimization. Default is `FALSE`.
#' @param min_popSize Minimum population size if adaptive population is enabled. Default is 20.
#' @param max_popSize Maximum population size if adaptive population is enabled. Default is 100.
#' @param early_stopping Logical. Whether to stop early if no improvement is observed. Default is `FALSE`.
#' @param tolerance Numeric tolerance for early stopping. Default is 1e-6.
#' @param patience Number of iterations to wait for improvement before stopping. Default is 10.
#'
#' @return A list containing the same elements as jaya_multi, plus additional
#'         constraint information stored as attributes.
#'
#' @importFrom stats runif
#' @export
jaya_multi_enhanced <- function(objectives, lower, upper, popSize = 50, maxiter, n_var,
                                seed = NULL, suggestions = data.frame(),
                                constraints = list(), adaptive_pop = FALSE,
                                min_popSize = 20, max_popSize = 100,
                                early_stopping = FALSE, tolerance = 1e-6,
                                patience = 10) {
  # Set seed for reproducibility
  set.seed(seed)

  # Initialize population
  tab <- data.frame(matrix(stats::runif((popSize - nrow(suggestions)) * n_var,
                                        min = rep(lower, each = popSize - nrow(suggestions)),
                                        max = rep(upper, each = popSize - nrow(suggestions))),
                           ncol = n_var))
  colnames(tab) <- paste0("x", 1:n_var)
  tab <- rbind(suggestions, tab)

  # Initialize constraint tracking
  n_constraints <- length(constraints)
  constraint_violations <- matrix(NA, nrow = nrow(tab), ncol = n_constraints)
  is_feasible <- rep(TRUE, nrow(tab))

  # Evaluate initial population and store as a list of vectors
  func_values <- lapply(seq_len(nrow(tab)), function(i) {
    # Calculate objective values
    obj_vals <- sapply(objectives, function(f) f(as.numeric(tab[i, 1:n_var])))

    # Check constraints
    if (n_constraints > 0) {
      constraint_vals <- sapply(constraints, function(con) {
        con(as.numeric(tab[i, 1:n_var]))
      })
      constraint_violations[i, ] <<- constraint_vals
      is_feasible[i] <<- all(constraint_vals <= 0)
    }

    return(obj_vals)
  })

  func_values_df <- do.call(rbind, func_values)
  colnames(func_values_df) <- paste0("Objective_", 1:length(objectives))
  tab <- cbind(tab, func_values_df)

  # Non-dominated sorting function to include decision variables with objective values
  non_dominated_sort <- function(values, decision_vars) {
    pareto_front <- list()
    for (i in 1:nrow(values)) {
      is_dominated <- FALSE
      for (j in 1:nrow(values)) {
        if (all(values[j, ] <= values[i, ]) && any(values[j, ] < values[i, ])) {
          is_dominated <- TRUE
          break
        }
      }
      if (!is_dominated) {
        pareto_front[[length(pareto_front) + 1]] <- c(decision_vars[i, ], values[i, ])
      }
    }
    do.call(rbind, pareto_front)
  }

  # Initialize constraint info storage
  constraint_info <- list(
    constraints = constraints,
    violations = constraint_violations,
    feasible = is_feasible
  )

  # Initial Pareto front (only considering feasible solutions)
  feasible_indices <- which(is_feasible)
  if (length(feasible_indices) > 0) {
    pareto_front <- non_dominated_sort(func_values_df[feasible_indices, ],
                                       tab[feasible_indices, 1:n_var])
  } else {
    # If no feasible solutions, use all solutions (will be refined later)
    pareto_front <- non_dominated_sort(func_values_df, tab[, 1:n_var])
  }

  # Early stopping parameters
  no_improvement_count <- 0
  best_pareto_size <- nrow(pareto_front)

  # Main optimization loop
  for (i in 1:maxiter) {
    r1 <- stats::runif(n_var)
    r2 <- stats::runif(n_var)

    # Adaptive population adjustment
    if (adaptive_pop && i %% 10 == 0) {
      current_popSize <- max(min_popSize, min(nrow(tab) + ifelse(i %% 20 == 0, 5, -5), max_popSize))
      tab <- tab[1:current_popSize, ]
      is_feasible <- is_feasible[1:current_popSize]
      constraint_violations <- constraint_violations[1:current_popSize, , drop = FALSE]
    }

    for (can in 1:nrow(tab)) {
      candidate <- as.numeric(tab[can, 1:n_var])

      # Select random Pareto front solutions for updating candidate
      feasible_pareto_indices <- which(is_feasible)
      if (length(feasible_pareto_indices) > 0) {
        best_index <- sample(feasible_pareto_indices, 1)
      } else {
        best_index <- sample(seq_len(nrow(tab)), 1)
      }

      worst_index <- which.max(rowSums(func_values_df) * (!is_feasible * 10 + 1))
      best_solution <- as.numeric(tab[best_index, 1:n_var])
      worst_solution <- as.numeric(tab[worst_index, 1:n_var])

      # Update candidate position based on Pareto best and worst
      updated_candidate <- candidate + r1 * (best_solution - candidate) - r2 * (worst_solution - candidate)

      # Enforce boundaries
      updated_candidate <- pmax(pmin(updated_candidate, upper), lower)

      # Check constraints
      constraint_vals <- sapply(constraints, function(con) {
        con(updated_candidate)
      })
      feasible <- all(constraint_vals <= 0)

      # Evaluate updated candidate for each objective
      new_func_vals <- sapply(objectives, function(obj) obj(updated_candidate))

      # Only proceed if there are no NA values in new_func_vals
      if (!any(is.na(new_func_vals))) {
        # For feasible solutions or if both are infeasible, use normal domination rules
        if ((feasible && is_feasible[can]) || (!feasible && !is_feasible[can])) {
          is_better <- all(new_func_vals <= func_values[[can]]) &&
            any(new_func_vals < func_values[[can]])
        } else {
          # Feasible solutions always dominate infeasible ones
          is_better <- feasible
        }

        if (is_better) {
          tab[can, 1:n_var] <- updated_candidate
          func_values_df[can, ] <- new_func_vals
          func_values[[can]] <- new_func_vals
          constraint_violations[can, ] <- constraint_vals
          is_feasible[can] <- feasible
        }
      }
    }

    # Update Pareto front with decision variables included (considering only feasible solutions)
    feasible_indices <- which(is_feasible)
    if (length(feasible_indices) > 0) {
      pareto_result <- non_dominated_sort(func_values_df[feasible_indices, ],
                                          tab[feasible_indices, 1:n_var])
    } else {
      # If no feasible solutions, use all solutions
      pareto_result <- non_dominated_sort(func_values_df, tab[, 1:n_var])
    }

    # Combine with previous Pareto front
    if (!is.null(pareto_result) && !is.null(pareto_front)) {
      if (ncol(pareto_result) == ncol(pareto_front)) {
        pareto_front <- unique(rbind(pareto_front, pareto_result))
      } else if (is.null(pareto_front)) {
        pareto_front <- pareto_result
      }
    }

    # Early stopping check
    if (early_stopping) {
      current_pareto_size <- nrow(pareto_front)
      if (abs(current_pareto_size - best_pareto_size) < tolerance) {
        no_improvement_count <- no_improvement_count + 1
      } else {
        no_improvement_count <- 0
        best_pareto_size <- current_pareto_size
      }
      if (no_improvement_count >= patience) {
        message("Early stopping triggered due to lack of improvement.")
        break
      }
    }
  }

  # Prepare final output with decision variables and objective values in Pareto front
  pareto_front_df <- data.frame(pareto_front)
  colnames(pareto_front_df) <- c(paste0("x", 1:n_var), paste0("Objective_", 1:length(objectives)))

  # Update constraint info for final population
  constraint_info$violations <- constraint_violations
  constraint_info$feasible <- is_feasible

  # Return result with updated format
  result <- list("Pareto_Front" = pareto_front_df, "Solutions" = tab)

  # Set attributes for summary
  attr(result, "popSize") <- popSize
  attr(result, "maxiter") <- maxiter
  attr(result, "n_var") <- n_var
  attr(result, "lower") <- lower
  attr(result, "upper") <- upper
  attr(result, "objectives") <- objectives
  attr(result, "constraints") <- constraints
  attr(result, "constraint_info") <- constraint_info

  class(result) <- "jaya_multi"
  return(result)
}

#' Plot Constraints and Pareto Front for Multi-Objective Optimization
#'
#' Visualizes the relationship between constraints and Pareto-optimal solutions
#' in multi-objective optimization problems. This function shows both feasible and
#' infeasible solutions, highlighting the impact of constraints on the solution space.
#'
#' @param x An object of class \code{jaya_multi} containing optimization results.
#' @param constraints_to_show Vector of constraint indices to highlight in the plot.
#'        If NULL (default), all constraints are shown.
#' @param objective_indices Vector of length 2 indicating which objectives to plot.
#'        Default is c(1, 2) for the first two objectives.
#' @param interactive Logical; if TRUE, creates an interactive plot using plotly.
#'        Default is FALSE.
#' @param ... Additional graphical parameters passed to the underlying plotting functions.
#'
#' @return A ggplot2 object or plotly object (if interactive=TRUE) visualizing the
#'         constraints and Pareto front.
#'
#' @importFrom stats na.omit
#' @importFrom ggplot2 ggplot aes geom_point geom_line theme_minimal labs scale_color_manual scale_shape_manual
#' @importFrom rlang .data
#' @export
plot_jaya_constraints <- function(x, constraints_to_show = NULL,
                                  objective_indices = c(1, 2),
                                  interactive = FALSE, ...) {
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (interactive && !requireNamespace("plotly", quietly = TRUE)) {
    warning("Package 'plotly' is needed for interactive plots. Falling back to static plot.",
            call. = FALSE)
    interactive <- FALSE
  }

  # Check if input is a jaya_multi object
  if (!inherits(x, "jaya_multi")) {
    stop("Input must be an object of class 'jaya_multi'",
         call. = FALSE)
  }

  # Extract Pareto front and all solutions
  pareto_front <- x$Pareto_Front
  all_solutions <- x$Solutions

  # Get constraint information
  constraint_info <- attr(x, "constraint_info")
  if (is.null(constraint_info)) {
    warning("No constraint information found in the jaya_multi object. Was this optimization run with the enhanced version?",
            call. = FALSE)
    # Create a basic plot without constraint information
    p <- ggplot2::ggplot() +
      ggplot2::aes(x = pareto_front[, paste0("Objective_", objective_indices[1])],
                   y = pareto_front[, paste0("Objective_", objective_indices[2])]) +
      ggplot2::geom_point(color = "blue", size = 3) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Pareto Front",
        x = paste("Objective", objective_indices[1]),
        y = paste("Objective", objective_indices[2])
      )

    if (interactive) {
      return(plotly::ggplotly(p, ...))
    } else {
      return(p)
    }
  }

  # Extract objective columns
  obj_cols <- grep("^Objective_", names(all_solutions), value = TRUE)

  # Check if the specified objective indices are valid
  if (length(obj_cols) < max(objective_indices)) {
    stop("Invalid objective indices. The problem has only ", length(obj_cols), " objectives.",
         call. = FALSE)
  }

  # Select constraints to show
  n_constraints <- length(constraint_info$constraints)
  if (is.null(constraints_to_show)) {
    constraints_to_show <- seq_len(n_constraints)
  } else if (any(constraints_to_show > n_constraints)) {
    stop("Invalid constraint indices. The problem has only ", n_constraints, " constraints.",
         call. = FALSE)
  }

  # Convert to data frame for plotting
  plot_data <- data.frame(
    Objective1 = all_solutions[, obj_cols[objective_indices[1]]],
    Objective2 = all_solutions[, obj_cols[objective_indices[2]]],
    InPareto = FALSE,
    Feasible = constraint_info$feasible
  )

  # Mark points in the Pareto front
  for (i in seq_len(nrow(pareto_front))) {
    for (j in seq_len(nrow(all_solutions))) {
      if (all(all_solutions[j, obj_cols] == pareto_front[i, obj_cols])) {
        plot_data$InPareto[j] <- TRUE
        break
      }
    }
  }

  # Create the plot using explicit .data$ references to avoid NOTES
  p <- ggplot2::ggplot(plot_data) +
    ggplot2::aes(x = .data$Objective1, y = .data$Objective2,
                 color = .data$Feasible, shape = .data$InPareto) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue"),
                                name = "Feasible") +
    ggplot2::scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                                name = "In Pareto Front") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Constraints and Pareto Front Visualization",
      x = paste("Objective", objective_indices[1]),
      y = paste("Objective", objective_indices[2])
    )

  # Add constraint violation information if available
  if (!is.null(constraint_info$violations)) {
    for (i in constraints_to_show) {
      # Extract constraint violations for the current constraint
      constraint_data <- data.frame(
        Objective1 = all_solutions[, obj_cols[objective_indices[1]]],
        Objective2 = all_solutions[, obj_cols[objective_indices[2]]],
        Violation = constraint_info$violations[, i]
      )

      # Add contour or boundary if possible (simplistic approach)
      # This would be enhanced in a full implementation
      boundary_points <- constraint_data[abs(constraint_data$Violation) < 0.01, ]
      if (nrow(boundary_points) >= 2) {
        # Sort boundary points to create a somewhat continuous line
        boundary_points <- boundary_points[order(boundary_points$Objective1), ]
        p <- p + ggplot2::geom_line(data = boundary_points,
                                    ggplot2::aes(x = .data$Objective1, y = .data$Objective2),
                                    color = "green", linetype = "dashed", size = 1)
      }
    }
  }

  # Return the plot
  if (interactive) {
    if (requireNamespace("plotly", quietly = TRUE)) {
      return(plotly::ggplotly(p, ...))
    } else {
      warning("Package 'plotly' is needed for interactive plots but not available. Returning static plot.")
      return(p)
    }
  } else {
    return(p)
  }
}
