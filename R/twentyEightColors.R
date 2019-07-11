# twentyEightColors ------------------------------------------------------
# I choose 28 colors built in R that can be used for scentific plotting.
# For examples, these colors can be used to represent different clusters
# in scRNA-seq.

library(gplot)

## Define 28 colors
twentyEightCol <- c("pink1", "violet", "mediumpurple1", "slateblue1", "purple", "purple3",
                    "turquoise2", "skyblue", "steelblue", "blue2", "navyblue",
                    "orange", "tomato", "coral2", "palevioletred", "violetred", "red2",
                    "springgreen2", "yellowgreen", "palegreen4",
                    "wheat2", "tan", "tan2", "tan3", "brown",
                    "grey70", "grey50", "grey30")

## Plot 28 colors
plot(1:28, rep(1, 28), type = "h", col = twentyEightCol,
     ylim = c(0, 1), lwd = 5)

## Convert colors to RGB
twentyEightColRGB <- col2rgb(twentyEightCol)

## Convert colors to HEX
twentyEightColHEX <- col2hex(twentyEightCol)