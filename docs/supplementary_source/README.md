# Source of the supplementary appendix, part B

`docs/supplementary_appendix.pdf` has two parts. Part A (tables S1–S54, figures S1–S4, sections S.1–S.13) reproduces
the appendix of an earlier version of the thesis. Part B (tables S55 onwards, figure S5, sections S.14–S.25,
equations B1 onwards) holds the material moved out of Appendix A of the final thesis to shorten it.

Part B is built from `addendum.tex` and `addendum_body.tex` in the thesis's LaTeX project: compile from the
project root after the thesis (`pdflatex -jobname=addendum supplementary/addendum.tex`), with
`supplementary/thesis_bibcites.tex` holding the thesis's `\bibcite` lines so that citation numbers match the thesis,
and `xr-hyper` resolving references to the thesis through `thesis.aux`.
