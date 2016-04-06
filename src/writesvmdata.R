writesvmdata <- function(model, cutoff) {
    fn <- "svmdata.h"

    # write license
    license <- c("/* Copyright 2015, 2016 Pol Cusco and Guillaume Filion",
                 "",
                 "   This file is part of Zerone.",
                 "",
                 "   Zerone is free software: you can redistribute it and/or modify",
                 "   it under the terms of the GNU General Public License as published by",
                 "   the Free Software Foundation, either version 3 of the License, or",
                 "   (at your option) any later version.",
                 "",
                 "   Zerone is distributed in the hope that it will be useful,",
                 "   but WITHOUT ANY WARRANTY; without even the implied warranty of",
                 "   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the",
                 "   GNU General Public License for more details.",
                 "",
                 "   You should have received a copy of the GNU General Public License",
                 "   along with Zerone. If not, see <http://www.gnu.org/licenses/>.",
                 "*/",
                 "")
    cat(license, file = fn, sep = "\n")

    # write defines
    defines <- c("#ifndef _SVM_DATA",
                 "#define _SVM_DATA",
                 "",
                 paste("#define NSV", model$tot.nSV),
                 paste("#define DIM", ncol(model$SV)),
                 paste("#define GAMMA", model$gamma),
                 paste("#define RHO", model$rho + cutoff),
                 "")
    cat(defines, file = fn, sep = "\n", append = TRUE)

    # write center values
    center <- c("const double CENTER[DIM] = {",
                paste(model$x.scale$`scaled:center`, collapse = ", "),
                "};",
                "")
    cat(center, file = fn, sep = "\n", append = TRUE)

    # write scale values
    scale <- c("const double SCALE[DIM] = {",
               paste(model$x.scale$`scaled:scale`, collapse = ", "),
               "};",
               "")
    cat(scale, file = fn, sep = "\n", append = TRUE)

    # write support vectors
    sv <- c("const double SV[NSV * DIM] = {",
            paste(model$SV, collapse = ", "),
            "};",
            "")
    cat(sv, file = fn, sep = "\n", append = TRUE)

    # write coefficients
    coefs <- c("const double COEFS[NSV] = {",
               paste(model$coefs, collapse = ", "),
               "};",
               "")
    cat(coefs, file = fn, sep = "\n", append = TRUE)

    # end define
    cat("#endif", file = fn, sep = "\n", append = TRUE)
}