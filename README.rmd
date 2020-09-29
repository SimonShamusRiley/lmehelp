## {lmehelp}: Convenience functions for the design and analysis of complex experiments

This effort is still in its earliest stage, but its ultimate intention is to:

1) provide a suite of tools for students and researchers to explore how specific designs related to the design of their experiment will affect replication and statistical power,

2) provide tools for properly specifying the linear predictor in a linear mixed model,

3) enable the correction of the degrees of freedom in a fitted lme model when "Dalgaard's
technique" is required because crossed random effects are present, 

4) create a toy that will help students develop an intuition for the relationship(s) between
experimental design, analysis and inference

## Known bugs and limitations

The package uses a modification of the containment algorithem as specified in the SAS documentation, but the syntatic matching procedure *requires* a match for all terms (an "error strata" must be found for all fixed effects; nothing is assigned residual degrees of freedom) which, at present allows for any arbitrary design that can be expressed in terms of nesting, crossing and interaction *except* for latin squares and their ilk. This will be addressed in time. 

While some minimal sanity checks are in place, it is still possible to specify utterly non-sensical designs with the keyout() function. As a "toy", this is not only acceptable but perhaps desirable, as it confronts the user with a (hopefully) obvious problem but without simply throwing errors that can, for the novice (in either R or experimental design) quickly become frustrating and disheartening. In so far as it is used as a tool, this aspect of the keyout() function is clearly less desirable. In any case, and indeed any endeavor, "measure twice and cut once".

To report additional bugs, to request features/functions or to get involved in the effort, [please reach out](simon.riley@ufl.edu)
