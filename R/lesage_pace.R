#' LeSage and Pace procedure
#'
#' @param base data
#' @param modele formula
#' @param matrice weights matrix
#' @param alpha confidence level between 0 and 1
#'
#' @return returns a dataframe of results
#'
#' @importFrom stats AIC lm
#' @import spatialreg
#' @export
#'
#' @examples

lesage_pace <- function(base,modele,matrice,alpha) {
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

    df<-data.frame(vs=c("SDM vs SAR","SAR vs OLS","SDM vs SLX","SLX vs OLS","SDM vs SEM","SEM vs OLS"),
                   stat=rep(0,6),
                   pvalue=rep(0,6))


    # H0 : teta = 0
    # Modele non contraint : SDM
    # Modele contraint : SAR
    # Si on rejette H0 alors on garde le modele non contraint, ie le modele SDM
    sdm_vs_sar<-LR.sarlm(sdm,sar)
    df$stat[which(df$vs=="SDM vs SAR")]<-round(sdm_vs_sar$statistic ,2)
    if (sdm_vs_sar$p.value<0.001) {df$pvalue[which(df$vs=="SDM vs SAR")]="< 0.001"} else {df$pvalue[which(df$vs=="SDM vs SAR")]=round(sdm_vs_sar$p.value,4)}

    if (sdm_vs_sar$p.value > alpha) {
        sar_vs_ols<-LR.sarlm(sar,ols)
        df$stat[which(df$vs=="SAR vs OLS")]<-round(sar_vs_ols$statistic ,2)
        if (sar_vs_ols$p.value<0.001) {df$pvalue[which(df$vs=="SAR vs OLS")]="< 0.001"} else {df$pvalue[which(df$vs=="SAR vs OLS")]=round(sar_vs_ols$p.value,4)}
    } else {
        df$stat[which(df$vs=="SAR vs OLS")]<-""
        df$pvalue[which(df$vs=="SAR vs OLS")]<-""
    }

    # H0 : rho=0, teta !=0 et teta + rho*beta != 0
    # Modele non contraint : SDM
    # Modele contraint : SLX
    # Si on rejette H0 alors on garde le modele non contraint, ie le modele SDM
    sdm_vs_slx<-LR.sarlm(sdm,slx)
    df$stat[which(df$vs=="SDM vs SLX")]<-round(sdm_vs_slx$statistic ,2)
    if (sdm_vs_slx$p.value<0.001) {df$pvalue[which(df$vs=="SDM vs SLX")]="< 0.001"} else {df$pvalue[which(df$vs=="SDM vs SLX")]=round(sdm_vs_slx$p.value,4)}

    if (sdm_vs_slx$p.value > alpha) {
        slx_vs_ols<-LR.sarlm(slx,ols)
        df$stat[which(df$vs=="SLX vs OLS")]<-round(slx_vs_ols$statistic ,2)
        if (slx_vs_ols$p.value<0.001) {df$pvalue[which(df$vs=="SLX vs OLS")]="< 0.001"} else {df$pvalue[which(df$vs=="SLX vs OLS")]=round(slx_vs_ols$p.value,4)}
    } else {
        df$stat[which(df$vs=="SLX vs OLS")]<-""
        df$pvalue[which(df$vs=="SLX vs OLS")]<-""
    }

    # H0 : teta + rho*beta= 0
    # Modele non contraint : SDM
    # Modele contraint : SEM
    # Si on rejette H0 alors on garde le modele non contraint, ie le modele SDM
    sdm_vs_sem<-LR.sarlm(sdm,sem)
    df$stat[which(df$vs=="SDM vs SEM")]<-round(sdm_vs_sem$statistic ,2)
    if (sdm_vs_sem$p.value<0.001) {df$pvalue[which(df$vs=="SDM vs SEM")]="< 0.001"} else {df$pvalue[which(df$vs=="SDM vs SEM")]=round(sdm_vs_sem$p.value,4)}

    if (sdm_vs_sem$p.value > alpha) {
        sem_vs_ols<-LR.sarlm(sem,ols)
        df$stat[which(df$vs=="SEM vs OLS")]<-round(sem_vs_ols$statistic ,2)
        if (sem_vs_ols$p.value<0.001) {df$pvalue[which(df$vs=="SEM vs OLS")]="< 0.001"} else {df$pvalue[which(df$vs=="SEM vs OLS")]=round(sem_vs_ols$p.value,4)}
    } else {
        df$stat[which(df$vs=="SEM vs OLS")]<-""
        df$pvalue[which(df$vs=="SEM vs OLS")]<-""
    }

    if (sdm_vs_sar$p.value <= alpha & sdm_vs_slx$p.value <= alpha & sdm_vs_sem$p.value <= alpha) { ##SAR_1 SLX_1 SEM_1
        mod<-data.frame(vs="Modele retenu : SDM", stat="", pvalue="")
    } else if (sdm_vs_sar$p.value > alpha & sdm_vs_slx$p.value <= alpha & sdm_vs_sem$p.value <= alpha) { ##SAR_0 SLX_1 SEM_1
        if (sar_vs_ols$p.value > alpha) {
            mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")
        } else {mod<-data.frame(vs="Modele retenu : SAR", stat="", pvalue="")}
    } else if (sdm_vs_sar$p.value <= alpha & sdm_vs_slx$p.value > alpha & sdm_vs_sem$p.value <= alpha) { ##SAR_1 SLX_0 SEM_1
        if (slx_vs_ols$p.value > alpha) {
            mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")
        } else {mod<-data.frame(vs="Modele retenu : SLX", stat="", pvalue="")}
    } else if (sdm_vs_sar$p.value <= alpha & sdm_vs_slx$p.value <= alpha & sdm_vs_sem$p.value > alpha) { ##SAR_1 SLX_1 SEM_0
        if (sem_vs_ols$p.value > alpha) {
            mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")
        } else {mod<-data.frame(vs="Modele retenu : SEM", stat="", pvalue="")}
    } else if (sdm_vs_sar$p.value > alpha & sdm_vs_slx$p.value > alpha & sdm_vs_sem$p.value <= alpha) { ##SAR_0 SLX_0 SEM_1
        if (sar_vs_ols$p.value > alpha & slx_vs_ols$pvalue > alpha) { #SAR_0 SLX_0
            mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")
        } else if (sar_vs_ols$p.value <= alpha & slx_vs_ols$pvalue > alpha) { #SAR_1 SLX_0
            if (AIC(sar) < AIC(ols)) {
                mod<-data.frame(vs="Modele retenu : SAR", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")}
        } else if (sar_vs_ols$p.value > alpha & slx_vs_ols$pvalue <= alpha) { #SAR_0 SLX_1
            if (AIC(slx) < AIC(ols)) {
                mod<-data.frame(vs="Modele retenu : SLX", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")}
        } else if (sar_vs_ols$p.value <= alpha & slx_vs_ols$pvalue <= alpha) { #SAR_1 SLX_1
            if (AIC(sar) < AIC(slx)) {
                mod<-data.frame(vs="Modele retenu : SAR", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : SLX", stat="", pvalue="")}
        }
    } else if (sdm_vs_sar$p.value > alpha & sdm_vs_slx$p.value <= alpha & sdm_vs_sem$p.value > alpha) { ##SAR_0 SLX_1 SEM_0
        if (sar_vs_ols$p.value > alpha & sem_vs_ols$p.value > alpha) { #SAR_0 SEM_0
            mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")
        } else if (sar_vs_ols$p.value <= alpha & sem_vs_ols$p.value > alpha) { #SAR_1 SEM_0
            if (AIC(sar) < AIC(ols)) {
                mod<-data.frame(vs="Modele retenu : SAR", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")}
        } else if (sar_vs_ols$p.value > alpha & sem_vs_ols$p.value <= alpha) { #SAR_0 SEM_1
            if (AIC(sem) < AIC(ols)) {
                mod<-data.frame(vs="Modele retenu : SEM", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")}
        } else if (sar_vs_ols$p.value <= alpha & sem_vs_ols$p.value <= alpha) { #SAR_1 SEM_1
            if (AIC(sar) < AIC(sem)) {
                mod<-data.frame(vs="Modele retenu : SAR", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : SEM", stat="", pvalue="")}
        }
    } else if (sdm_vs_sar$p.value <= alpha & sdm_vs_slx$p.value > alpha & sdm_vs_sem$p.value > alpha) { ##SAR_1 SLX_0 SEM_0
        if (slx_vs_ols$pvalue > alpha & sem_vs_ols$p.value > alpha) { #SLX_0 SEM_0
            mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")
        } else if (slx_vs_ols$pvalue <= alpha & sem_vs_ols$p.value > alpha) { #SLX_1 SEM_0
            if (AIC(slx) < AIC(ols)) {
                mod<-data.frame(vs="Modele retenu : SLX", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")}
        } else if (slx_vs_ols$pvalue > alpha & sem_vs_ols$p.value <= alpha) { #SLX_0 SEM_1
            if (AIC(sem) < AIC(ols)) {
                mod<-data.frame(vs="Modele retenu : SEM", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")}
        } else if (slx_vs_ols$pvalue <= alpha & sem_vs_ols$p.value <= alpha) { #SLX_1 SEM_1
            if (AIC(slx) < AIC(sem)) {
                mod<-data.frame(vs="Modele retenu : SLX", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : SEM", stat="", pvalue="")}
        }
    } else if (sdm_vs_sar$p.value <= alpha & sdm_vs_slx$p.value <= alpha & sdm_vs_sem$p.value <= alpha) { ##SAR_1 SLX_1 SEM_1
        if (sar_vs_ols$p.value > alpha & slx_vs_ols > alpha & sem_vs_ols > alpha) { # SAR_0 SLX_0 SEM_0
            mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")
        } else if (sar_vs_ols$p.value <= alpha & slx_vs_ols > alpha & sem_vs_ols > alpha) { # SAR_1 SLX_0 SEM_0
            if (AIC(sar) < AIC(ols)) {
                mod<-data.frame(vs="Modele retenu : SAR", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")}
        } else if (sar_vs_ols$p.value > alpha & slx_vs_ols <= alpha & sem_vs_ols > alpha) { # SAR_0 SLX_1 SEM_0
            if (AIC(slx) < AIC(ols)) {
                mod<-data.frame(vs="Modele retenu : SLX", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")}
        } else if (sar_vs_ols$p.value <= alpha & slx_vs_ols > alpha & sem_vs_ols > alpha) { # SAR_0 SLX_0 SEM_1
            if (AIC(sem) < AIC(ols)) {
                mod<-data.frame(vs="Modele retenu : SEM", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : OLS", stat="", pvalue="")}
        } else if (sar_vs_ols$p.value <= alpha & slx_vs_ols <= alpha & sem_vs_ols > alpha) { # SAR_1 SLX_1 SEM_0
            if (AIC(sar) < AIC(slx)) {
                mod<-data.frame(vs="Modele retenu : SAR", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : SLX", stat="", pvalue="")}
        } else if (sar_vs_ols$p.value <= alpha & slx_vs_ols > alpha & sem_vs_ols <= alpha) { # SAR_1 SLX_0 SEM_1
            if (AIC(sar) < AIC(sem)) {
                mod<-data.frame(vs="Modele retenu : SAR", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : SEM", stat="", pvalue="")}
        } else if (sar_vs_ols$p.value > alpha & slx_vs_ols <= alpha & sem_vs_ols <= alpha) { # SAR_0 SLX_1 SEM_1
            if (AIC(slx) < AIC(sem)) {
                mod<-data.frame(vs="Modele retenu : SLX", stat="", pvalue="")
            } else {mod<-data.frame(vs="Modele retenu : SEM", stat="", pvalue="")}
        } else if (sar_vs_ols$p.value <= alpha & slx_vs_ols <= alpha & sem_vs_ols <= alpha) { # SAR_1 SLX_1 SEM_1
            if (AIC(sar) < AIC(slx) & AIC(sar) < AIC(sem)) {
                mod<-data.frame(vs="Modele retenu : SAR", stat="", pvalue="")
            } else if  (AIC(slx) < AIC(sar) & AIC(slx) < AIC(sem)) {
                mod<-data.frame(vs="Modele retenu : SLX", stat="", pvalue="")
            } else if (AIC(sem) < AIC(sar) & AIC(sem) < AIC(slx)) {
                mod<-data.frame(vs="Modele retenu : SEM", stat="", pvalue="")
            }
        }
    }

    df<-rbind.data.frame(df,mod)

    colnames(df)<-c("Procedure de Lesage Pace","Statistique de test","P-value")
    return(df)
}
