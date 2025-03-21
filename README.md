# Jaya Package Enhancement: Constraint Visualization

This repository contains an enhancement to the [Jaya](https://CRAN.R-project.org/package=Jaya) R package, adding new functionality for visualizing constraints in multi-objective optimization problems.

## Overview

The Jaya algorithm is a gradient-free, population-based optimization technique that can solve both single-objective and multi-objective problems without requiring algorithm-specific parameters. This enhancement adds the ability to visualize how constraints affect the Pareto front in multi-objective optimization scenarios.

## New Features

This enhancement includes:

1. **Enhanced Multi-Objective Optimization**: `jaya_multi_enhanced()` function that tracks constraint information throughout the optimization process
2. **Constraint Visualization**: `plot_jaya_constraints()` function that visualizes:
   - The Pareto front in objective space
   - Feasible and infeasible solutions
   - Constraint boundaries and their impact on the solution space

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("plato-12/jaya-enhancement")
```

## Usage

### Multi-Objective Optimization with Constraints

```r
library(Jaya)

# Define objectives
objective1 <- function(x) sum(x^2)
objective2 <- function(x) sum((x - 2)^2)

# Define constraints (must return <= 0 for feasibility)
constraint1 <- function(x) sum(x) - 3  # Sum of variables must be <= 3
constraint2 <- function(x) -sum(x) + 1 # Sum of variables must be >= 1

# Run optimization with constraint tracking
result <- jaya_multi_enhanced(
  objectives = list(objective1, objective2),
  lower = c(-5, -5),
  upper = c(5, 5),
  popSize = 50,
  maxiter = 100,
  n_var = 2,
  constraints = list(constraint1, constraint2)
)

# Visualize constraints and Pareto front
plot_jaya_constraints(result)

# Focus on specific constraints
plot_jaya_constraints(result, constraints_to_show = 1)

# Create interactive visualization
plot_jaya_constraints(result, interactive = TRUE)
```

## Examples

The constraint visualization shows:
- Blue points: Feasible solutions in the Pareto front
- Red points: Infeasible solutions that violate constraints
- Green dashed lines: Boundaries of constraints in the objective space

This visualization helps users understand:
- Which constraints are most restrictive
- How constraints shape the Pareto front
- Where the trade-offs exist between objectives under constraints

## Documentation

See the vignette "Visualizing Constraints in Multi-Objective Optimization" for detailed examples and explanations:

```r
vignette("Visualizing_Constraints_in_Multi-Objective_Optimization", package = "Jaya")
```

## Requirements

- R (>= 3.5.0)
- ggplot2
- plotly (optional, for interactive visualizations)

## License

MIT License

## Author

Priyanshu Tiwari (contributor)

Original Jaya package by Neeraj Bokde
