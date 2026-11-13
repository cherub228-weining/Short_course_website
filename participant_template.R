# Name / group: ______________________
# Run from the course folder root. Fill the TODO sections.
source("code/R/exercise_setup.R")
# Task 1: report units, dates, rows and state support.
# TODO
# Task 2: estimate impact and plot against the known truth.
impact <- estimate_horizon(0, Jset=4:10)
grid <- seq(-.8,.8,length.out=161)
curve <- infer_curve(impact, grid)
# TODO: plot estimate and 1-grid^2
# Task 3: contrast and uncertainty.
contrast <- infer_contrast(impact,zA=0,zB=.75)
print(contrast)
# TODO: compare pointwise and band widths at a grid point
# Task 4: sensitivity and interpretation.
fixed <- estimate_horizon(0,Jset=7)
# TODO: compare infer_contrast(fixed) with contrast
# Optional: repeat at horizon 4 and change HAC bandwidth.


# Optional E: endogenous state (still use its pre-shock value).
# endog <- simulate_panel(state="endogenous")
# efit <- estimate_horizon(0,Jset=4:10,data=endog)
# infer_contrast(efit)
# Optional F: compare finite-shift target and LP slope.
# nonlinear_shift_demo(z=0,delta=1,cubic=.25)
