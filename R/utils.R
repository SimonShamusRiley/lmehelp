# Count number of levels present in a data set for a factor or interaction term
# (for ease of use in -apply functions), appropriately ordered and grouped, with
# containment denomenator degrees of freedom
nlevel = function(x, data){sum(attr(model.matrix(formula(paste('~ 0 +', x)), data), 'assign'))}

# Syntactical matching of random effects and fixed effects terms
match.eu = function(term, des.terms) {
  term.split = strsplit(term, ':')
  matches = lapply(term.split, function(x) {sapply(x, function(x, des.terms) regexpr(x, des.terms) > 0, des.terms = des.terms)})
  match.num = vector('list', length(matches))

  for (n in 1:length(matches)){
    match.num[[n]] = ifelse(!is.array(matches[[n]]), which(mean(matches[[n]]) == 1), which(rowMeans(matches[[n]]) == 1))
  }

  out =  des.terms[sapply(match.num, min)]
  out
}

# Identify the list entry containing all individual terms in a vector
match.chars = function(terms, list){
  out = sapply(terms, function(x, list) {regexpr(x, list) > 0}, list = list)
  if (!is.array(out)){
    out = which(mean(out) == 1)
  } else {
    out = which(rowMeans(out) == 1)
  }
  out
}
