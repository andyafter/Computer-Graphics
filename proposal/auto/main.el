(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "12pt")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "amscd"
    "amsmath"
    "amssymb"
    "amsthm"
    "cite"
    "epsfig"
    "verbatim"
    "graphicx"
    "color")
   (LaTeX-add-labels
    "fig:dipole"
    "fig:directional")
   (LaTeX-add-bibliographies
    "papers"))
 :latex)

