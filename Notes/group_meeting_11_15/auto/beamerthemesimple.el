(TeX-add-style-hook
 "beamerthemesimple"
 (lambda ()
   (TeX-run-style-hooks
    "xkeyval"
    "tikz")
   (TeX-add-symbols
    "setwatermark"))
 :latex)

