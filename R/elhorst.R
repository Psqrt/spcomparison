#' Elhorst procedure
#'
#' @param base data
#' @param modele formula
#' @param matrice weights matrix
#' @param alpha confidence level between 0 and 1
#'
#' @return returns a dataframe of results
#' @importFrom spdep lm.LMtests
#' @export
#'
#' @examples

elhorst <- function(base,modele,matrice,alpha) {

    # Modele SDM
    sdm<-lagsarlm(modele, data=base, matrice, type="mixed")
    # Modele SAR
    sar<-lagsarlm(modele, data=base, matrice)
    # Modele SEM
    sem<-errorsarlm(modele, data=base, matrice)
    # Modele SLX (Modele a interaction exogene)
    slx<-lmSLX(modele, data=base, matrice)
    # Modele OLS
    ols<-lm(modele,data = base)

    vs <-c()
    stat <- c()
    pvalue <- c()

    # Test LM-Error : H0 : lambda = 0
    # Si on accepte H0 alors lambda est egal a 0
    lambda<-lm.LMtests(ols,matrice,test="RLMerr")
    vs <- c(vs, "Test LM-Error")
    stat <- c(stat, round(lambda$RLMerr$statistic,2))
    pvalue <- c(pvalue, round(lambda$RLMerr$p.value,4))

    # Test LM-Lag : H0 : rho = 0
    # Si on accepte H0 alors rho est egale a 0
    rho<-lm.LMtests(ols,matrice,test="RLMlag")
    vs <- c(vs, "Test LM-Lag")
    stat <- c(stat, round(rho$RLMlag$statistic,2))
    if (rho$RLMlag$p.value < 0.001) {pvalue <- c(pvalue, "< 0.001")} else {pvalue <- c(pvalue, round(rho$RLMlag$p.value,4))}


    if (rho$RLMlag$p.value <= alpha & lambda$RLMerr$p.value > alpha) {
        sdm_vs_sar<-LR.sarlm(sdm,sar)
        vs <- c(vs, "SDM vs SAR")
        stat <- c(stat, round(sdm_vs_sar$statistic,2))
        if (sdm_vs_sar$p.value < 0.001) {pvalue <- c(pvalue, "< 0.001")} else {pvalue <- c(pvalue, round(sdm_vs_sar$p.value,4))}
        if (sdm_vs_sar$p.value > alpha) {
            vs <- c(vs, "Modele retenu : SAR")
            stat <- c(stat, "")
            pvalue <- c(pvalue, "")
        } else {
            vs <- c(vs, "Modele retenu : SDM")
            stat <- c(stat, "")
            pvalue <- c(pvalue, "")
        }
    } else if (rho$RLMlag$p.value > alpha & lambda$RLMerr$p.value > alpha) {
        slx_vs_ols<-LR.sarlm(slx,ols)
        vs <- c(vs, "SLX vs OLS")
        stat <- c(stat, round(slx_vs_ols$statistic,2))
        if (slx_vs_ols$p.value < 0.001) {pvalue <- c(pvalue, "< 0.001")} else {pvalue <- c(pvalue, round(slx_vs_ols$p.value,4))}
        if (slx_vs_ols$p.value > alpha) {
            vs <- c(vs, "Modele retenu : OLS")
            stat <- c(stat, "")
            pvalue <- c(pvalue, "")
        } else {
            sdm_vs_slx<-LR.sarlm(sdm,slx)
            vs <- c(vs, "SDM vs SLX")
            stat <- c(stat, round(sdm_vs_slx$statistic,2))
            if (sdm_vs_slx$p.value < 0.001) {pvalue <- c(pvalue, "< 0.001")} else {pvalue <- c(pvalue, round(sdm_vs_slx$p.value,4))}
            if (sdm_vs_slx$p.value > alpha) {
                vs <- c(vs, "Modele retenu : SLX")
                stat <- c(stat, "")
                pvalue <- c(pvalue, "")
            } else {
                vs <- c(vs, "Modele retenu : SDM")
                stat <- c(stat, "")
                pvalue <- c(pvalue, "")
            }
        }
    } else if (rho$RLMlag$p.value > alpha & lambda$RLMerr$p.value <= alpha) {
        sdm_vs_sem<-LR.sarlm(sdm,sem)
        vs <- c(vs, "SDM vs SEM")
        stat <- c(stat, round(sdm_vs_sem$statistic,2))
        if (sdm_vs_sem$p.value < 0.001) {pvalue <- c(pvalue, "< 0.001")} else {pvalue <- c(pvalue, round(sdm_vs_sem$p.value,4))}
        if (sdm_vs_sem$p.value > alpha) {
            vs <- c(vs, "Modele retenu : SEM")
            stat <- c(stat, "")
            pvalue <- c(pvalue, "")
        } else {
            vs <- c(vs, "Modele retenu : SDM")
            stat <- c(stat, "")
            pvalue <- c(pvalue, "")
        }
    }

    df<-data.frame(vs, stat, pvalue)
    colnames(df)<-c("Procedure d'Elhorst","Statistique de test","P-value")

    return(df)
}
