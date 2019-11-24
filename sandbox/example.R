library(rgdal)
library(spcomparison)
columbus <- readOGR(system.file("shapes/columbus.shp", package="spData")[1])
weights <- nb2listw(poly2nb(columbus))

formula <- CRIME ~ INC + HOVAL
elhorst(base = columbus,
        modele = CRIME ~ INC + HOVAL,
        matrice = weights,
        alpha = 0.05)

lesage_pace(data = columbus, formula = formula,
                 listw = weights, alpha = 0.05, criterion = "AIC")


