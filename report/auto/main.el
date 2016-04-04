(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("beamer" "serif" "mathserif" "professionalfont" "hyperref={pdfpagelabels=false}")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("gnuplottex" "miktex") ("adjustbox" "export")))
   (add-to-list 'LaTeX-verbatim-environments-local "semiverbatim")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "beamer"
    "beamer10"
    "exscale"
    "amsmath"
    "graphicx"
    "helvet"
    "tcolorbox"
    "textpos"
    "xargs"
    "tikz"
    "amssymb"
    "algorithm2e"
    "pxfonts"
    "eulervm"
    "gnuplottex"
    "adjustbox"
    "cite"
    "mathtools")
   (TeX-add-symbols
    '("norm" 1)
    "localtextbulletone"
    "mynew"
    "gdw")))

