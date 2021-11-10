(TeX-add-style-hook
 "padicl"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("hyperref" "pagebackref=true") ("algorithm2e" "ruled" "vlined") ("appendix" "toc" "page") ("tocbibind" "nottoc") ("footmisc" "stable") ("unicode-math" "math-style=ISO" "bold-style=ISO")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "amsmath"
    "amsthm"
    "amssymb"
    "amsfonts"
    "mathrsfs"
    "enumitem"
    "adjustbox"
    "caption"
    "setspace"
    "xspace"
    "hyperref"
    "algorithm2e"
    "xcolor"
    "bbm"
    "rotating"
    "array"
    "esint"
    "nicefrac"
    "cleveref"
    "tikz-cd"
    "tikz"
    "appendix"
    "natbib"
    "tocbibind"
    "float"
    "footmisc"
    "yfonts"
    "kpfonts"
    "fontspec"
    "unicode-math"
    "garamondx"
    "tocloft"
    "etoc"
    "physics")
   (TeX-add-symbols
    '("ls" 2)
    '("Norm" 1)
    '("pow" 1)
    '("ps" 1)
    '("Jap" 1)
    "A"
    "N"
    "Z"
    "Q"
    "R"
    "G"
    "C"
    "Ca"
    "F"
    "p"
    "Lip"
    "lav"
    "Lav"
    "loc"
    "nsub"
    "Set"
    "Ab"
    "Sc"
    "MHS"
    "Aff"
    "Coh"
    "Ring"
    "Comp"
    "dR"
    "shHom"
    "shExt"
    "defeq"
    "mbb"
    "mscr"
    "mc"
    "mbf"
    "mf")
   (LaTeX-add-labels
    "sec:introduction"
    "sec:background"
    "eq:8"
    "eq:9"
    "eq:10"
    "eq:12"
    "eq:11"
    "eq:13"
    "sec:overview"
    "sec:Theory"
    "eq:17"
    "eq:19"
    "eq:20"
    "eg:Jacobi-sum"
    "eq:24"
    "eq:25"
    "eq:27"
    "eq:18"
    "eq:21"
    "eq:Hecke-fn-eq"
    "eq:23"
    "eq:2"
    "eq:1"
    "eq:3"
    "sec:p-adic-interpolation"
    "sec:algorithm"
    "sec:lambda-invariants"
    "sec:quick-summary-ejv"
    "sec:ment-results-sant"
    "eq:14"
    "sec:implementation"
    "sec:quadr-forms-meth"
    "eq:4"
    "sec:comp-sets-nearly"
    "eq:5"
    "sec:optim-naive-algor"
    "sec:future-research")
   (LaTeX-add-bibliographies
    "/home/havard/Documents/bibliography/references")
   (LaTeX-add-amsthm-newtheorems
    "thm"
    "cor"
    "conj"
    "lemma"
    "prop"
    "mydef"
    "eg"
    "remark"
    "exercise")
   (LaTeX-add-xcolor-definecolors
    "phoen"
    "purpur"
    "skog"
    "hav"
    "himmel"))
 :latex)

