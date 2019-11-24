#' Elhorst procedure
#'
#' @param formula a formula object (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data a data.frame object containing the variables in the model.
#' @param listw a listw object containing weights.
#' @param alpha a numeric value providing the confidence level between 0 and 1.
#' @param criterion either RLM (default), AIC or BIC: in case if several models are selected, using additional criterion to draw out the best model.
#'
#' @return returns a list of two dataframes of test results.
#'
#' @references Elhorst, J. P. (2010). Applied spatial econometrics: raising the bar. Spatial economic analysis, 5(1), 9-28.
#'
#' @importFrom stats AIC BIC lm printCoefmat
#' @import spatialreg
#' @importFrom spdep lm.LMtests
#' @export
#'
#' @examples
#' require(spatialreg)
#' require(spdep)
#' require(rgdal)
#' columbus <- readOGR(system.file("shapes/columbus.shp", package="spData")[1])
#' weights <- nb2listw(poly2nb(columbus))
#' formula <- CRIME ~ INC + HOVAL
#'
#' elhorst(formula = formula, data = columbus, listw = weights, alpha = 0.05)
#' elhorst(formula = formula, data = columbus,
#'         listw = weights, alpha = 0.05, criterion = "AIC")
#' elhorst(formula = CRIME ~ INC, data = columbus, listw = weights)

elhorst <- function(formula, data, listw, alpha = 0.05, criterion = "RLM") {

    ### VERIFICATION DES INPUTS
    if (!inherits(formula, "formula")) stop("No formula given", call. = F)
    if (!inherits(listw, "listw")) stop("No neighbourhood list", call. = F)
    if (!inherits(alpha, "numeric")) stop("Confidence level alpha is not numeric", call. = F)
    if (nrow(data) != length(listw$neighbours)) stop("Input data and weights have different dimensions", call. = F)
    if (abs(alpha) > 1) stop("Confidence level alpha should be between 0 and 1", call. = F)
    if (!(criterion %in% c("AIC", "BIC", "RLM"))) stop("Criterion should be 'AIC', 'BIC' or 'RLM'", call. = F)

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

    # DF CONTENANT LES RESULTATS DES TESTS (R)LM
    df_lm <- data.frame(test = c("LMlag", "LMerr", "RLMlag", "RLMerr"),
                        stat = rep(NA, 4),
                        pvalue = rep(NA, 4))

    # DF CONTENANT LES RESULTATS DES TESTS LR
    df_lr <- data.frame(vs = c("SDM vs SAR", "SLX vs OLS", "SDM vs SLX", "SDM vs SEM"),
                        stat = rep(NA, 4),
                        pvalue = rep(NA, 4))

    # CONTIENDRA LES MODELES RETENUS
    selected_models <- c(NA, NA)

    # Test LM-Lag : H0 : rho = 0
    # Si on accepte H0 alors rho est egale a 0
    rho <- lm.LMtests(ols, listw, test = "LMlag")
    row <- which(df_lm["test"] == "LMlag")
    df_lm$stat[row] <- rho$LMlag$statistic
    df_lm$pvalue[row] <- rho$LMlag$p.value

    robust_rho <- lm.LMtests(ols, listw, test = "RLMlag")
    row <- which(df_lm["test"] == "RLMlag")
    df_lm$stat[row] <- robust_rho$RLMlag$statistic
    df_lm$pvalue[row] <- robust_rho$RLMlag$p.value

    # Test LM-Error : H0 : lambda = 0
    # Si on accepte H0 alors lambda est egal a 0
    lambda <- lm.LMtests(ols, listw, test = "LMerr")
    row <- which(df_lm["test"] == "LMerr")
    df_lm$stat[row] <- lambda$LMerr$statistic
    df_lm$pvalue[row] <- lambda$LMerr$p.value

    robust_lambda <- lm.LMtests(ols, listw, test = "RLMerr")
    row <- which(df_lm["test"] == "RLMerr")
    df_lm$stat[row] <- robust_lambda$RLMerr$statistic
    df_lm$pvalue[row] <- robust_lambda$RLMerr$p.value

    if (rho$LMlag$p.value < alpha){ # CAS GAUCHE
        sdm_vs_sar <- LR.sarlm(sdm, sar)
        row <- which(df_lr$vs == "SDM vs SAR")
        df_lr$stat[row] <- sdm_vs_sar$statistic
        df_lr$pvalue[row] <- sdm_vs_sar$p.value
        if (sdm_vs_sar$p.value < alpha){
            selected_models[1] <- "SDM"
        } else {
            selected_models[1] <- "SAR"
        }
    }

    if (rho$LMlag$p.value > alpha & lambda$LMerr$p.value > alpha){ # CAS CENTRAL
        slx_vs_ols <- LR.sarlm(slx, ols)
        row <- which(df_lr$vs == "SLX vs OLS")
        df_lr$stat[row] <- slx_vs_ols$statistic
        df_lr$pvalue[row] <- slx_vs_ols$p.value

        if (slx_vs_ols$p.value < alpha){
            sdm_vs_slx <- LR.sarlm(sdm, slx)
            row <- which(df_lr$vs == "SDM vs SLX")
            df_lr$stat[row] <- sdm_vs_slx$statistic
            df_lr$pvalue[row] <- sdm_vs_slx$p.value
            if (sdm_vs_slx$p.value < alpha){
                selected_models[1:2] <- "SDM"
            } else {
                selected_models[1:2] <- "SLX"
            }
        } else {
            selected_models[1:2] <- "OLS"
        }
    }


    if (lambda$LMerr$p.value < alpha){ # CAS DROITE
        sdm_vs_sem <- LR.sarlm(sdm, sem)
        row <- which(df_lr$vs == "SDM vs SEM")
        df_lr$stat[row] <- sdm_vs_sem$statistic
        df_lr$pvalue[row] <- sdm_vs_sem$p.value
        if (sdm_vs_sem$p.value < alpha){
            selected_models[2] <- "SDM"
        } else {
            selected_models[2] <- "SEM"
        }
    }

    selected_models <- selected_models[!is.na(selected_models)] # IL EXISTE NA SI UN SEUL MODELE RETENU PARMI 2
    selected_models <- unique(selected_models) # RETIRER DOUBLONS
    if (length(selected_models) == 1){ # SI 1 SEUL MODEL, PAS BESOIN DE DEPARTAGER
        bestmodel <- selected_models
    } else { # SINON COMPARER AIC BIC OU RLM
        if (criterion == "AIC"){
            info <- sapply(selected_models, function(x) { return(AIC(get(tolower(x))))})
            bestmodel <- names(info[order(info)][1])
        } else if (criterion == "BIC") {
            info <- sapply(selected_models, function(x) { return(BIC(get(tolower(x))))})
            bestmodel <- names(info[order(info)][1])
        } else {
            if (robust_lambda$RLMerr$p.value > robust_rho$RLMlag$p.value) {
                bestmodel <- selected_models[1]
            } else {
                bestmodel <- selected_models[2]
            }
        }
    }
    rownames(df_lm) <- df_lm$test
    df_lm["test"] <- NULL
    colnames(df_lm) <- c("LM-stat", "p-value")

    rownames(df_lr) <- df_lr$vs
    df_lr["vs"] <- NULL
    colnames(df_lr) <- c("Likelihood ratio", "p-value")


    cat("Formula:\n")
    print(formula)
    cat("\nWeights matrix:\n")
    cat(deparse(substitute(listw)))
    cat("\n\n(Robust) Lagrange multiplier tests:\n")
    printCoefmat(df_lm, has.Pvalue = TRUE, signif.stars = TRUE, signif.legend = TRUE, na.print = "-")
    cat("\n--------------------------------------------------------------")
    cat("\nLikelihood ratio tests:\n")
    printCoefmat(df_lr, has.Pvalue = TRUE, signif.stars = TRUE, signif.legend = TRUE, na.print = "-")
    if (length(selected_models) == 1){
        cat("\n--------------------------------------------------------------\nOnly one model was selected:\n")
        cat(bestmodel, "\n")
    } else {
        cat("\n--------------------------------------------------------------\nSeveral models were selected:\n")
        cat(selected_models, sep = ", ")
        if (criterion %in% c("AIC", "BIC")){
            cat("\n\nSelected models ", criterion, ":\n", sep = "")
            print(info)
            cat("\nSelected model according to ", criterion, ":\n", sep = "")
            cat(bestmodel, "\n")
        } else {
            cat("\n\nSelected model according to robust LM tests:\n")
            cat(bestmodel, "\n")
        }
    }

    return(invisible(list(df_lm, df_lr)))
}
