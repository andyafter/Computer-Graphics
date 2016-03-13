(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "dvips" "12pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("graphicx" "pdftex")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "graphicx"
    "url")
   (LaTeX-add-labels
    "m42")
   (LaTeX-add-bibitems
    "gonzalez2012")))

