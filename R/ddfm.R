#' @title Containment denominator degrees of freedom for a fitted lme model with crossed random effects
#'
#' @param lme.object A fitted lme object.
#' @param design A one-sided formula specifying the design structure.
#' The design structure must be fully specified: see details for info.
#' @details Nelder (1965) notes that crossing two factors is equivalent to
#'  nesting them at the same level (https://doi.org/10.1098/rspa.1965.0012).
#'  Dalgaard (2006), among others, have demonstrated how this can be leveraged to
#'  fit models with crossed random effects in the nlme::lme() function despite
#'  the positional notation which the package uses in specifying random effects
#'  which are assumed to be nested
#'  (https://staff.pubhealth.ku.dk/~pd/mixed-jan.2006/lme.pdf).
#'
#'  This enables the fitting of models with any arbitrary random effects structure, with
#'  patterned variance-covariance structures, unequal variances, and correlated
#'  residuals, which is not possible in other widely-used mixed modeling R
#'  packages (i.e. lme4). The available techniques for doing so, however,
#'  "trick" the algorithm the packages uses for calculating denominator degrees
#'  of freedom, (almost always) leading to inflated p-values. This function allows
#'  the calculation of the correct degrees of freedom based on the specified design
#'  structure.
#'
#'  The formula specifying the design structure must be fully specified, including
#'  a replication factor nested within treatment(s) for designs lacking a blocking
#'  structure.
#'
#'  The relevant operators are: "/" to indicate nesting, "*" to indicate crossing, and
#'  ":" to indicate an interaction (useful for factorial experiments where the experimental
#'  unit is at the level of an interaction among several factors).
#'
#'  Thus, for a 2-factor factorial experiment with replication number given by R:
#'
#'  CRD design: ~ A:B/R (experimental units are nested at the level of A:B interaction)
#'
#'  RCB design: ~ R/A:B (experimental units at the level of A:B interaction are nested in blocks)
#'
#'  Split-Plot RCB: ~ R/A/B
#'
#'  Split-Block RCB: ~R/(A*B)
#'
#' @return A suitably ordered vector of denominator degrees of freedom
#' corresponding to each fixed effect term in the lme.model, calculated using
#' the containment method.
#' @examples
#' # Adapted from page 425 of Steel, Torrie & Dickey (1997), an "Analysis of Variance for a Split-
#' # Plot in Space and Time":
#' data = expand.grid(R = factor(1:4), A = factor(1:5), B = factor(1:3), Y = factor(1:3))
#' keyout(trt.formula = ~ A*B*Y, des.formula = ~ R/(A/B*Y), data)

ddfm = function(lme.object, design){
  if (class(lme.object) != 'lme') {stop('\nlme.object must be model fitted with nlme:lme()')}
  if (class(design) != 'formula') {stop('\ndesign structure must be specified as a formula')}
  trt.form = formula(lme.object)
  key = keyout(formula(lme.object), des.formula = design, data = lme.object$data)
  for (n in (nrow(key)-1):1){
    key$ErrDF[n] = ifelse(is.na(key$ErrDF[n]), key$ErrDF[n+1], key$ErrDF[n])
  }
  term.names = names(lme.object$fixDF$terms)
  ddfm.new = rep(NA, length(term.names))
  for (n in 1:length(term.names)){
    temp = key$ErrDF[which(key$factors %in% term.names[n])]
    if (length(temp) == 0){
      ddfm.new[n] = key$ErrDF[which(is.na(key$factors)[1])]
    } else {
      ddfm.new[n] = temp
    }
  }
  ddfm.new
}
