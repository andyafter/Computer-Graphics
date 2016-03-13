(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("beamer" "serif" "mathserif" "professionalfont" "hyperref={pdfpagelabels=false}")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("gnuplottex" "miktex") ("adjustbox" "export")))
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

