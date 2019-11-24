#' LeSage and Pace procedure
#'
#' @param formula a formula object (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data a data.frame object containing the variables in the model.
#' @param listw a listw object containing weights.
#' @param alpha a numeric value providing the confidence level between 0 and 1.
#' @param criterion either AIC (default) or BIC: in case if several models are selected, using information criterion to draw out the best model.
#'
#' @return returns a data.frame object containing likelihood ratio tests information.
#'
#' @references LeSage, J., & Pace, R. K. (2009). Introduction to spatial econometrics. Chapman and Hall/CRC.
#'
#' @importFrom stats AIC lm printCoefmat
#' @import spatialreg
#' @export
#'
#' @examples
#' require(spatialreg)
#' columbus <- rgdal::readOGR(system.file("shapes/columbus.shp", package="spData")[1])
#' weights <- spdep::nb2listw(spdep::poly2nb(columbus))
#' formula <- CRIME ~ INC + HOVAL
#'
#' lesage_pace(formula = formula, data = columbus,
#'             listw = weights, alpha = 0.05, criterion = "AIC")
#' lesage_pace(formula = INC ~ CRIME, data = columbus, listw = weights)

lesage_pace <- function(formula, data, listw, alpha = 0.05, criterion = "AIC") {

    ### ICI GESTION DES ERREURS DES INPUTS #############################################
    if (!inherits(formula, "formula")) stop("No formula given", call. = F)
    if (!inherits(listw, "listw")) stop("No neighbourhood list", call. = F)
    if (!inherits(alpha, "numeric")) stop("Confidence level alpha is not numeric", call. = F)
    if (nrow(data) != length(listw$neighbours)) stop("Input data and weights have different dimensions", call. = F)
    if (abs(alpha) > 1) stop("Confidence level alpha should be between 0 and 1", call. = F)
    if (!(criterion %in% c("AIC", "BIC"))) stop("Criterion should be 'AIC' or 'BIC'", call. = F)

    # Modele SDM
    sdm <- lagsarlm(formula = formula, data = data, listw = listw, type = "mixed")
    # Modele SAR
    sar <- lagsarlm(formula = formula, data = data, listw = listw)
    # Modele SEM
    sem <- errorsarlm(formula = formula, data = data, listw = listw)
    # Modele SLX (Modele a interaction exogene)
    slx <- lmSLX(formula = formula, data = data, listw = listw)
    # Modele OLS
    ols <- lm(formula = formula, data = data)

    df <- data.frame(vs = c("SDM vs SAR", "SAR vs OLS", "SDM vs SLX", "SLX vs OLS", "SDM vs SEM", "SEM vs OLS"),
                     stat = rep(0, 6),
                     pvalue = rep(0, 6))
    selected_models <- c() # RÉCUPÈRE LES MODÈLES QUI SONT RETENUS PAR BRANCHE

    # H0 : theta = 0
    # Modele non contraint : SDM
    # Modele contraint : SAR
    # Si on rejette H0 alors on garde le modele non contraint, ie le modele SDM
    sdm_vs_sar <- LR.sarlm(sdm, sar)
    row <- which(df$vs == "SDM vs SAR") # POUR EVITER DE RÉPÉTER DEUX FOIS
    df$stat[row] <- sdm_vs_sar$statistic
    df$pvalue[row] <- sdm_vs_sar$p.value

    if (sdm_vs_sar$p.value > alpha) {
        sar_vs_ols <- LR.sarlm(sar, ols)
        row <- which(df$vs == "SAR vs OLS")
        df$stat[row] <- sar_vs_ols$statistic
        df$pvalue[row] <- sar_vs_ols$p.value
        if (sar_vs_ols$p.value < alpha){
            selected_models <- c(selected_models, "SAR") # ON RECUPERE POUR PLUS TARD
        } else {
            selected_models <- c(selected_models, "OLS")
        }
    } else {
        df$stat[row] <- NA
        df$pvalue[row] <- NA
        selected_models <- c(selected_models, "SDM")
    }

    # H0 : rho=0, teta !=0 et teta + rho*beta != 0
    # Modele non contraint : SDM
    # Modele contraint : SLX
    # Si on rejette H0 alors on garde le modele non contraint, ie le modele SDM
    sdm_vs_slx <- LR.sarlm(sdm, slx)
    row <- which(df$vs == "SDM vs SLX")
    df$stat[row] <- sdm_vs_slx$statistic
    df$pvalue[row] <- sdm_vs_slx$p.value

    if (sdm_vs_slx$p.value > alpha) {
        slx_vs_ols <- LR.sarlm(slx, ols)
        row <- which(df$vs == "SLX vs OLS")
        df$stat[row] <- slx_vs_ols$statistic
        df$pvalue[row] <- slx_vs_ols$p.value
        if (slx_vs_ols$p.value < alpha){
            selected_models <- c(selected_models, "SLX")
        } else {
            selected_models <- c(selected_models, "OLS")
        }
    } else {
        df$stat[row] <- NA
        df$pvalue[row] <- NA
        selected_models <- c(selected_models, "SDM")
    }

    # H0 : teta + rho*beta= 0
    # Modele non contraint : SDM
    # Modele contraint : SEM
    # Si on rejette H0 alors on garde le modele non contraint, ie le modele SDM
    sdm_vs_sem <- LR.sarlm(sdm, sem)
    row <- which(df$vs == "SDM vs SEM")
    df$stat[row] <- sdm_vs_sem$statistic
    df$pvalue[row] <- sdm_vs_sem$p.value

    if (sdm_vs_sem$p.value > alpha) {
        sem_vs_ols <- LR.sarlm(sem, ols)
        row <- which(df$vs == "SEM vs OLS")
        df$stat[row] <- sem_vs_ols$statistic
        df$pvalue[row] <- sem_vs_ols$p.value
        if (sem_vs_ols$p.value < alpha){
            selected_models <- c(selected_models, "SEM")
        } else {
            selected_models <- c(selected_models, "OLS")
        }
    } else {
        df$stat[row] <- NA
        df$pvalue[row] <- NA
        selected_models <- c(selected_models, "SDM")
    }

    selected_models <- unique(selected_models) # ON RETIRE LES DOUBLONS (DEUX BRANCHES PEUVENT DONNER OLS PAR EXEMPLE)
    if (length(selected_models) == 1){
        bestmodel <- selected_models # SI Y A QU'UN SEUL MODELE, C'EST LE MEILLEUR
    } else { # SINON S'IL Y A PLUSIEURS, FAUT TESTER PAR AIC BIC
        if (criterion == "AIC"){
            info <- sapply(selected_models, function(x) { return(AIC(get(tolower(x))))})
        } else {
            info <- sapply(selected_models, function(x) { return(BIC(get(tolower(x))))})
        }
        bestmodel <- names(info[order(info)][1]) # ON PREND CELUI QUI MINIMISE
    }
    rownames(df) <- df$vs
    df["vs"] <- NULL # ON RETIRE LA COLONNE VS CAR C'EST PASSÉ EN ROWNAMES
    colnames(df) <- c("Likelihood ratio","p-value")

    cat("Formula:\n")
    print(formula)
    cat("\nWeights matrix:\n")
    cat(deparse(substitute(listw)))
    cat("\n\nLikelihood ratio tests:\n")
    printCoefmat(df, has.Pvalue = T) # C'EST CETTE FONCTION QUI MET LES TABLEAUX EN FORME AVEC LES ÉTOILES ET TOUT
    if (length(selected_models) == 1){
        cat("\n--------------------------------------------------------------\nOnly one model was selected:\n")
        cat(bestmodel, "\n")
    } else {
        cat("\n--------------------------------------------------------------\nSeveral models were selected:\n")
        cat(selected_models, sep = ", ")
        cat("\n\nSelected models ", criterion, ":\n", sep = "")
        print(info)
        cat("\nSelected model according to ", criterion, ":\n", sep = "")
        cat(bestmodel, "\n")
    }
    return(invisible(df)) # INVISIBLE PERMET DE NE PAS PRINT(DF) QUAND LE RÉSULTAT N'EST PAS ASSIGNÉ À UN OBJET
}

