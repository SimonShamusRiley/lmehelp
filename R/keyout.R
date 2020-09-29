#' @title Generate an ANOVA keyout for a given treatment and design structure
#'
#' @param trt.formula A one-side formula specifying the treatment structure. See
#' details for formula specification. Response variables are ignored if included.
#' @param des.formula A one-sided formula specifying the design structure.
#' The design structure must be fully specified to the level of experimental
#' unit. See details for information.
#' @param  data A data frame containing all factors and interactions included
#' in the treatment and design formula
#' @details Stroup (2012) has noted how the ANOVA keyout can be a powerful modeling tool
#'  (ISBN 9781439815120). The approach can also be useful in helping students to gain
#'  an intuition for the implications of various experimental designs on statistical power
#'  and the impact of replication (or lack there of) at various scales/error strata.
#'  This function therefore takes a treatment and design structure, and a data frame and
#'  produces an anova keyout, usefully ordered and grouped, with error strata and both
#'  numerator and denominator degrees of freedom. Note that where data is not yet available
#'  (because the study is hypothetical or in the design phase), the expand.grid() function
#'  can be used to quickly and easily simulate data sets.
#'
#'  The relevant operators are: "/" to indicate nesting, "*" to indicate crossing, and
#'  ":" to indicate an interaction (useful for factorial experiments where the experimental
#'  unit is at the level of an interaction among several factors).
#'
#'  Thus, for a 2-factor factorial experiment the treatment structure would be ~ A*B, whereas
#'  if the levels of the B factor are nested in A, the treatment structure would be ~ A/B.
#'
#'  Similarly, design structures for common (simple) designs are sepcified as follows:
#'
#'  CRD design: ~ A:B/R (experimental units are nested at the level of A:B interaction)
#'
#'  RCB design: ~ R/A:B (experimental units at the level of A:B interaction are nested in blocks)
#'
#'  Split-Plot RCB: ~ R/A/B
#'
#'  Split-Block RCB: ~R/(A*B)
#'
#' @return A data frame containing the fixed effects, levels, numerator DF, error strata and
#' denominator DF. For all fixed effects sharing the same experimental unit/error strata, only
#' the last in the list will have the error term or den DF, with others shown as NA.
#' @examples
#' library(agridat)
#' data("cox.stripsplit")
#' # fit the model in lme using "Dalgaard's technique":
#' fit = lme(fixed = yield~ soil*fert*calcium, data = cox.stripsplit,
#' random = list(rep = ~1, rep = pdIdent(form = ~soil-1), rep = pdIdent(form = ~ fert-1), calcium = ~1))
#' # denDF are wrong:
#' fit$fixDF$terms
#' # fix DF
#' fit$fixDF$terms = ddfm(fit, ~ rep/(fert/calcium * soil))
#' # denDF are correct:
#' fit$fixDF$terms
keyout = function(trt.formula, des.formula, data){
  if (class(trt.formula) != 'formula' | class(des.formula) != 'formula') {stop('\nTreatment and design structures must be specified as formulae')}
  if (!is.data.frame(data)) {stop('\nData must be provide as a data.frame')}
  trt.terms = attr(terms.formula(trt.formula), 'term.labels')
  trt.levels = sapply(trt.terms, nlevel, data, USE.NAMES = F)
  if (any(trt.levels == 1)) {warning('\nOne or more treatments has only a single level:\nVerify that all treatments specified in trt.formula are factors in the data frame')}
  trt.df = table(attr(model.matrix(trt.formula, data), 'assign'))[-1]

  trt.frame =data.frame(factors = trt.terms, levels = trt.levels, DF = as.numeric(trt.df))

  des.terms = attr(terms.formula(des.formula), 'term.labels')
  des.terms = des.terms[!as.character(des.terms) %in% trt.terms]
  if (length(des.terms) == 0) {stop('\nDesign formula contains no unique elements\nNote that the design structure must be fully specified,\nincluding a replication factor nested in treatment(s) if no blocking factor is present\nSee documentation for details.')}
  des.terms = ordered(des.terms, levels = des.terms)
  des.levels = sapply(des.terms, nlevel, data, USE.NAMES = F)
  if (any(des.levels == 1)) {warning('\nOne or more design components has only a single level:\nVerify that all replication/blocking factors specified in des.formula are factors in the data frame')}

  ErrDF = sapply(des.terms, nlevel, data, USE.NAMES = F)
  Err = data.frame(Error = des.terms, ErrDF)

  trt.frame$Error = match.eu(term = trt.frame$factors, des.terms)
  trt.frame = dplyr::right_join(trt.frame, Err)

  trt.frame = trt.frame[order(trt.frame$Error), ]
  trt.frame$ErrDF = trt.frame$ErrDF - 1
  trt.frame$Error = as.character(trt.frame$Error)

  des.term.list = strsplit(trt.frame$Error, ':')
  N = length(des.term.list)

  matches = lapply(des.term.list, match.chars, list = des.term.list)
  match.num = vector('list', N)
  for (n in 1:N){
    match.num[[n]] = which(sapply(matches, function(x, n) {n %in% x}, n=n))
  }

  trt.frame$ErrDF = ifelse(rev(duplicated(rev(trt.frame$Error))), NA, trt.frame$ErrDF)
  trt.frame$Error = ifelse(rev(duplicated(rev(trt.frame$Error))), '', as.character(trt.frame$Error))

  for (n in 2:N){
    if (is.na(trt.frame$ErrDF[n])) {next}
    randfs = match.num[[n]][match.num[[n]] < n]
    trtdf = match.num[[n]]
    trt.frame$ErrDF[n] = trt.frame$ErrDF[n] - sum(trt.frame$ErrDF[randfs], na.rm = T) - sum(trt.frame$DF[trtdf], na.rm = T)
  }
  if (any(trt.frame$ErrDF[!is.na(trt.frame$ErrDF)] == 0)) {warning('\nOne or more error DF is zero\nVerify that all replication/blocking factors specified in des.formula are factors in the data frame')}
  trt.frame
}
