#This is the script that I am using to analyze my corn kernel qPCR assay
#10112023 MLF


###### CUSTOM FUNCTIONS ######

#To calculate ng per test sample --------
abs_quant <- function(C,m,b){
  10^((C-b)/m)
}


# Bootstrapping function -------
boot.fn <- function(x, N=5000) {
  Int.1 <- replicate(N, mean(sample(x, size= length(x), replace=T)))
  Int.CI <- quantile(Int.1, probs=c(0.025,0.975))
  Int.CI
}


###### PREPPING WORKSPACE ######

x <- c("dplyr","ggpubr", "ggplot2", "data.table","gridExtra", "nlme") 
lapply(x, FUN = function(X) {do.call("library", list(X))}) #loading libraries

setwd("~/Downloads")

data <- read.csv("Vip_qPCR_2ndDer_10172023.csv", header = T)


#sanity check for data entry - all tests should have 4 except for one refuge sample, two per std, 4 negs 
table(data$Name)

#Getting NA samples, which were not detectable
data_NoAmp <- data[is.na(data$Cp),]
table(data_NoAmp$Name)# includes negatives

data_clean <- data[,c(4:8,11)]
data_clean <- subset(data_clean, Name != "neg")

#separating standard curve from test samples
data_SC <- subset(data_clean, Treat == "Std")
data_SC$Name <- as.numeric(data_SC$Name)

data_test <- subset(data_clean, Treat != "Std")
str(data_test)

#updating data structure
data_test$Primers <- as.factor(data_test$Primers)
data_test$Treat <- replace(data_test$Treat, data_test$Treat == "Vip", "AVip")#force reordering factor levels for analysis
data_test$Treat <- as.factor(data_test$Treat)
data_test$Name <- as.factor(data_test$Name)


###### STANDARD CURVE ANALYSES ######

#separating standard curves by primer pair
data_SC_hmg <- subset(data_SC, Primers == "hmg")
data_SC_Vip <- subset(data_SC, Primers == "Vip")


#model for hmg - the control gene
hmg_reg <- lm(Cp ~ log10(Name), data = data_SC_hmg)
summary(hmg_reg)

hmg_coef <- as.matrix(hmg_reg$coefficients)#pull m & b values from here.

plot(data_SC_hmg$Cp ~ log10(data_SC_hmg$Name), pch = 16, 
     cex = 1.3, col = "blue", xlab = "Log10 DNA Conc (ng)", ylab = "Cp value")
abline(lm(Cp ~ log10(Name), data = data_SC_hmg))


#model for Vip - the test gene
Vip_reg <- lm(Cp ~ log10(Name), data = data_SC_Vip)
summary(Vip_reg)

Vip_coef <- as.matrix(Vip_reg$coefficients)#pull m & b values from here.

plot(data_SC_Vip$Cp ~ log10(data_SC_Vip$Name), pch = 16, 
     cex = 1.3, col = "red", xlab = "Log10 DNA Conc (ng)", ylab = "Cp value")
abline(lm(Cp ~ log10(Name), data = data_SC_Vip))


###### ANALYSIS OF COPY NUMBER ACCORDING TO STEIBEL ET AL. 2009 ######

#direct analysis Cp values

data_test$CpNA <- data_test$Cp
data_test$CpNA <- ifelse(is.na(data_test$CpNA), 30, data_test$CpNA)#replacing high NAs with 30 (late cycle call)
fm1 <- lme(CpNA ~ Treat*Primers, random = ~Tech_rep|Name, data = data_test)
summary(fm1)

anova(fm1)


#analysis of estimated init quantities of DNA

#hmg primers
data_test_hmg <- subset(data_test, Primers == "hmg")
data_test_hmg$quant <- abs_quant(data_test_hmg$Cp, hmg_coef[2], hmg_coef[1])

#vip primers
data_test_Vip <- subset(data_test, Primers == "Vip")

#separating no amp vip from amped vip
data_Vip_NoAmp <- data_test_Vip[is.na(data_test_Vip$Cp),]
data_Vip_NoAmp$quant <- 0

data_test_Vip <- data_test_Vip[!is.na(data_test_Vip$Cp),]
data_test_Vip$quant <- abs_quant(data_test_Vip$Cp, Vip_coef[2], Vip_coef[1])

data_test_full <- rbind(data_test_hmg, data_test_Vip, data_Vip_NoAmp)#putting datasets back together

fm2 <- lme(quant ~ Treat*Primers, random = ~Tech_rep|Name, data = data_test_full)
summary(fm2)

anova(fm2)


###### EDA PLOTS OF TEST SAMPLES ######
data_test_full <- data_test_full[,c(1:4,8)]#removing CpNA

#calculating copy number

data_test_hmg <- subset(data_test_full, Primers == "hmg")

hmg_reshape <- reshape(data_test_hmg, idvar = c("Name", "Primers", "Treat"), 
                       timevar = "Tech_rep", direction = "wide")

hmg_reshape$meanQuant <- rowMeans(hmg_reshape[,c('quant.1', 'quant.2')], na.rm=FALSE)
hmg_reshape$Copies <- (hmg_reshape$meanQuant/0.002725)

hmg_meansByTreat <- tapply(hmg_reshape$Copies, hmg_reshape$Treat, mean)
hmg_meansByTreat

hmg_bootByTreat <- tapply(hmg_reshape$Copies, hmg_reshape$Treat, boot.fn)
hmg_bootByTreat


full_vip <- subset(data_test_full, Primers == "Vip")
Vip_reshape <- reshape(full_vip, idvar = c("Name", "Primers", "Treat"), 
                       timevar = "Tech_rep", direction = "wide")


Vip_reshape$meanQuant <- rowMeans(Vip_reshape[,c('quant.1', 'quant.2')], na.rm=FALSE)
Vip_reshape$Copies <- (Vip_reshape$meanQuant/0.002725)

Vip_meansByTreat <- tapply(Vip_reshape$Copies, Vip_reshape$Treat, mean)
Vip_meansByTreat

Vip_bootByTreat <- tapply(Vip_reshape$Copies, Vip_reshape$Treat, boot.fn)
Vip_bootByTreat



###### COMPARATIVE PLOTS NUM COPIES BY TREATMENT ######

#assuming haploid genome weight of maize to be 1 copy = 2.725 pg (from Liang et al. 2014 pp. 2605)

pdf("Vip_and_Control_DNA_Value_boxplot.pdf", width = 11, height = 8)

par(mfrow = c(1,2))

boxplot(Copies~Treat,
        data=Vip_reshape,
        main=substitute(paste(italic('vip3Aa20'))),
        xaxt = "n",
        xlab="Treatment",
        ylab="Estimated DNA Copies",
        col="orange",
        border="brown")
axis(1, at=1:3, labels=c("Vip", "Close Refuge", "Far Refuge"))

boxplot(Copies~Treat,
          data=hmg_reshape,
          main=substitute(paste(italic('hmg'))),
          xaxt = "n",
          ylim=c(0,14000),
          xlab="Treatment",
          ylab="Estimated DNA Copies",
          col="lightblue",
          border="navy")
axis(1, at=1:3, labels=c("Vip", "Close Refuge", "Far Refuge"))


dev.off()








