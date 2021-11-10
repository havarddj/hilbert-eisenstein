(TeX-add-style-hook
 "padicL_slides"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("beamer" "10pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("appendix" "toc" "page")))
   (add-to-list 'LaTeX-verbatim-environments-local "semiverbatim")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "beamer"
    "beamer10"
    "amsmath"
    "amsthm"
    "amssymb"
    "amsfonts"
    "mathrsfs"
    "mathabx"
    "epsfig"
    "cleveref"
    "xcolor"
    "bbm"
    "array"
    "esint"
    "nicefrac"
    "tikz-cd"
    "tikz"
    "enumitem"
    "subcaption"
    "wrapfig"
    "appendix"
    "physics"
    "natbib"
    "float")
   (TeX-add-symbols
    '("Pause" ["argument"] 0)
    '("ls" 2)
    '("Norm" 1)
    '("pow" 1)
    '("ps" 1)
    '("pur" 1)
    "A"
    "N"
    "Z"
    "Q"
    "R"
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
    "Aff"
    "Coh"
    "Ring"
    "Comp"
    "Mod"
    "defeq"
    "supp"
    "ord"
    "Nm"
    "Reg"
    "esssup"
    "spn"
    "Sl"
    "Gl"
    "GL"
    "Id"
    "Aut"
    "Hom"
    "Gal"
    "Pic"
    "Char"
    "AC"
    "OP"
    "pSh"
    "Sh"
    "Spec"
    "Fun"
    "im"
    "et"
    "Cl"
    "acts"
    "mbb"
    "mscr"
    "mc"
    "mf"
    "mbf"))
 :latex)

