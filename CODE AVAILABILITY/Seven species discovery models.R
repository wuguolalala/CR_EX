

########################################CRanimals##############################################################################
################################################################################################################################
data <- read.csv("D:\\other study\\model5\\CRanimalsmodel.csv")

library(nlme)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(broom)
library(minpack.lm)  


start.list <- list(A = max(data$seq.cum), x0 = mean(data$Year), B = sd(data$Year))  # 初始值列表
model <- nlsLM(seq.cum ~ A / (1 + exp(-(Year - x0)/B)), data = data, start = start.list)  # nlsLM函数


bic <- BIC(model)


print(paste("BIC: ", bic))



new_data <- data.frame(Year = c(1993:2022, 2023:2042))
new_data$predicted_seq_cum <- predict(model, newdata = new_data)
new_data$predicted_seq_cum <- round(new_data$predicted_seq_cum,0)


observed <- data$seq.cum
predicted <- predict(model)
fit <- cor(observed, predicted)^2


model_eq <- paste("seq.cum =", round(coef(model)[1], 3), "/ (1 + exp(-(", round(coef(model)[2], 3), "* (Year - ", round(coef(model)[3], 3), "))))")


ggplot(data, aes(x = Year, y = seq.cum)) +
  geom_point(color = "blue") +
  geom_line(data = new_data, aes(x = Year, y = predicted_seq_cum), color = "red") +
  annotate("text", x = 2025, y = 80000, label = paste("R^2 =", round(fit, 3)), color = "black", hjust = 4, vjust = 1) +
  annotate("text", x = 2025, y = 70000, label = model_eq, color = "black", hjust = 0.8, vjust = -1) +
  labs(title = "CR animals", x = "Year", y = "seq.cum")


write.csv(new_data, file = "D:\\other study\\model5\\CRanimalsN2.csv", row.names = FALSE)  
############################################################################################################################################
# ModelI-a


model1 <- gnls(sp.per ~ (a + b * sp.cum) * (3797 - sp.cum),
               data = data,
               start = list(a = 1e-02, b = 1e-04),
               weights = varPower(),
               verbose = TRUE,
               control = gnlsControl(returnObject = TRUE, minScale = 1e-500,
                                     tolerance = 0.001, nlsMaxIter = 3))


parameters <- coef(model1)
confidence_intervals <- confint(model1)
summary(model1)


yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[29+i])
  sp.cum.predict <- predict(model1, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model1)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval1 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval1)

################################################################################################################################################
#ModelI-b
#k(t)=a+b×S+c×S^2

model2=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2)*(3797-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))



parameters <- coef(model2)
confidence_intervals <- confint(model2)
summary(model2)



yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[29+i])
  sp.cum.predict <- predict(model2, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model2)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval2 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval2)

############################################################################################################################################
# ModelI-c

model3 <- gnls(sp.per ~ (a + b * sp.cum + c * log(seq.cum)) * (3797 - sp.cum),
               data = data,
               start = list(a = 1e-02, b = 1e-04, c = 1e-04),
               weights = varPower(),
               verbose = TRUE,
               control = gnlsControl(returnObject = TRUE, minScale = 1e-500, tolerance = 0.001, nlsMaxIter = 3))


parameters <- coef(model3)
confidence_intervals <- confint(model3)
summary(model3)

# Preparing the initial value of sp.cum.t
sp.cum.t <- tail(data$sp.cum, 1)


yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[29+i])
  sp.cum.predict <- predict(model3, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model3)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval3 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval3)
################################################################################################
#ModelI-d
#k(t)=a+b×S+c×log(N)+d×S×log(N)

model4=gnls(sp.per~(a+b*sp.cum+c*log(seq.cum)+d*sp.cum*log(seq.cum))*(3797-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))



parameters <- coef(model4)
confidence_intervals <- confint(model4)
summary(model4)


yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[29+i])
  sp.cum.predict <- predict(model4, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model4)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval4 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval4)

####################################################################################################################

# Model I-e
#k(t)=a+b×S+c×S^2+d×log(N)

model5=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2+d*log(seq.cum))*(3797-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))




parameters <- coef(model5)
confidence_intervals <- confint(model5)
summary(model5)

yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[29+i])
  sp.cum.predict <- predict(model5, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model5)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval5 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval5)

###########################################################################################################################
# ModelI-f
#k(t)=a+b×S+c×log(N)+d×log(N)^2
# Parameter estimation 


model6=gnls(sp.per~(a+b*sp.cum+c*log(seq.cum)+d*log(seq.cum)^2)*(3797-sp.cum),
             data=data,
             start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
             weights=varPower(),verbose=T,
             control=gnlsControl(returnObject=T,minScale=1e-500,
                                 tolerance=0.001,nlsMaxIter=3))




parameters <- coef(model6)
confidence_intervals <- confint(model6)
summary(model6)


yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[29+i])
  sp.cum.predict <- predict(model6, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model6)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval6 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval6)

##########################################################################################################################

# ModelI-g
#k(t)=a+b×S+c×S^2+d×log(N)+e×log(N)^2
# Parameter estimation 


model7=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2+d*log(seq.cum)+e*log(seq.cum)^2)*(3797-sp.cum),
             data=data,
             start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01, e=0.01),
             weights=varPower(),verbose=T,
             control=gnlsControl(returnObject=T,minScale=1e-500,
                                 tolerance=0.001,nlsMaxIter=3))


parameters <- coef(model7)
confidence_intervals <- confint(model7)
summary(model7)


yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[29+i])
  sp.cum.predict <- predict(model7, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model7)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval7 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval7)

#######################################################################################################################

combined_results <- cbind(prediction_interval1, prediction_interval2, prediction_interval3, prediction_interval4, prediction_interval5, prediction_interval6,prediction_interval7,prediction_interval8,prediction_interval9,prediction_interval10,prediction_interval11,prediction_interval12)
colnames(combined_results) <- c("Year", "Model1", "Model2", "Model3", "Model4", "Model5","Model6","Model7")
# 导出为CSV文件
write.csv(combined_results, file = "D:\\other study\\model5\\CRanimalsresult.csv", row.names = FALSE)




########################################CRplants################################################################################
################################################################################################################################
data <- read.csv("D:\\other study\\model5\\CRplantsmodel.csv")

library(nlme)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(broom)
library(minpack.lm)  


start.list <- list(A = max(data$seq.cum), x0 = mean(data$Year), B = sd(data$Year))  # 初始值列表
model <- nlsLM(seq.cum ~ A / (1 + exp(-(Year - x0)/B)), data = data, start = start.list)  # nlsLM函数


bic <- BIC(model)


print(paste("BIC: ", bic))



new_data <- data.frame(Year = c(1993:2022, 2023:2042))
new_data$predicted_seq_cum <- predict(model, newdata = new_data)
new_data$predicted_seq_cum <- round(new_data$predicted_seq_cum,0)


observed <- data$seq.cum
predicted <- predict(model)
fit <- cor(observed, predicted)^2


model_eq <- paste("seq.cum =", round(coef(model)[1], 3), "/ (1 + exp(-(", round(coef(model)[2], 3), "* (Year - ", round(coef(model)[3], 3), "))))")


ggplot(data, aes(x = Year, y = seq.cum)) +
  geom_point(color = "blue") +
  geom_line(data = new_data, aes(x = Year, y = predicted_seq_cum), color = "red") +
  annotate("text", x = 2025, y = 30000, label = paste("R^2 =", round(fit, 3)), color = "black", hjust = 4, vjust = 1) +
  annotate("text", x = 2025, y = 20000, label = model_eq, color = "black", hjust = 0.8, vjust = -10) +
  labs(title = "CR plants", x = "Year", y = "seq.cum")
write.csv(new_data, file = "D:\\other study\\濒危物种研究\\物种发现数预测\\model5\\CRplantsN2逻辑.csv", row.names = FALSE) 

############################################################################################################################################
# ModelI-a
#sp.per~(a+b*sp.cum)*(5232-sp.cum
# Input data

model1 <- gnls(sp.per ~ (a + b * sp.cum) * (5232 - sp.cum),
               data = data,
               start = list(a = 1e-02, b = 1e-04),
               weights = varPower(),
               verbose = TRUE,
               control = gnlsControl(returnObject = TRUE, minScale = 1e-500,
                                     tolerance = 0.001, nlsMaxIter = 3))


parameters <- coef(model1)
confidence_intervals <- confint(model1)
summary(model1)



yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[28+i])
  sp.cum.predict <- predict(model1, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model1)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval1 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval1)

################################################################################################################################################
#ModelI-b
#k(t)=a+b×S+c×S^2

model2=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2)*(5232-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))



parameters <- coef(model2)
confidence_intervals <- confint(model2)
summary(model2)



yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[28+i])
  sp.cum.predict <- predict(model2, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model2)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval2 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval2)

############################################################################################################################################
# ModelI-c

model3 <- gnls(sp.per ~ (a + b * sp.cum + c * log(seq.cum)) * (5232 - sp.cum),
               data = data,
               start = list(a = 1e-02, b = 1e-04, c = 1e-04),
               weights = varPower(),
               verbose = TRUE,
               control = gnlsControl(returnObject = TRUE, minScale = 1e-500, tolerance = 0.001, nlsMaxIter = 3))

parameters <- coef(model3)
confidence_intervals <- confint(model3)
summary(model3)

# Preparing the initial value of sp.cum.t
sp.cum.t <- tail(data$sp.cum, 1)


yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[28+i])
  sp.cum.predict <- predict(model3, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model3)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval3 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval3)
################################################################################################
#ModelI-d
#k(t)=a+b×S+c×log(N)+d×S×log(N)

model4=gnls(sp.per~(a+b*sp.cum+c*log(seq.cum)+d*sp.cum*log(seq.cum))*(5232-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))




parameters <- coef(model4)
confidence_intervals <- confint(model4)
summary(model4)


yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[28+i])
  sp.cum.predict <- predict(model4, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model4)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval4 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval4)

####################################################################################################################

# ModelI-e
#k(t)=a+b×S+c×S^2+d×log(N)
# Parameter estimation 

model5=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2+d*log(seq.cum))*(5232-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))




parameters <- coef(model5)
confidence_intervals <- confint(model5)
summary(model5)


yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[28+i])
  sp.cum.predict <- predict(model5, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model5)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval5 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval5)

###########################################################################################################################

# ModelI-f
#k(t)=a+b×S+c×log(N)+d×log(N)^2
# Parameter estimation 


model6=gnls(sp.per~(a+b*sp.cum+c*log(seq.cum)+d*log(seq.cum)^2)*(5232-sp.cum),
             data=data,
             start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
             weights=varPower(),verbose=T,
             control=gnlsControl(returnObject=T,minScale=1e-500,
                                 tolerance=0.001,nlsMaxIter=3))




parameters <- coef(model6)
confidence_intervals <- confint(model6)
summary(model6)


yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[28+i])
  sp.cum.predict <- predict(model6, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model6)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval6 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval6)

##########################################################################################################################
# ModelI-g
#k(t)=a+b×S+c×S^2+d×log(N)+e×log(N)^2
# Parameter estimation 


model7=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2+d*log(seq.cum)+e*log(seq.cum)^2)*(5232-sp.cum),
             data=data,
             start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01, e=0.01),
             weights=varPower(),verbose=T,
             control=gnlsControl(returnObject=T,minScale=1e-500,
                                 tolerance=0.001,nlsMaxIter=3))



parameters <- coef(model7)
confidence_intervals <- confint(model7)
summary(model7)



yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[28+i])
  sp.cum.predict <- predict(model7, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model7)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval7 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval7)

#######################################################################################################################

combined_results <- cbind(prediction_interval1, prediction_interval2, prediction_interval3, prediction_interval4, prediction_interval5, prediction_interval6,prediction_interval7,prediction_interval8,prediction_interval9,prediction_interval10,prediction_interval11,prediction_interval12)
colnames(combined_results) <- c("Year", "Model1", "Model2", "Model3", "Model4", "Model5","Model6","Model7")

write.csv(combined_results, file = "D:\\other study\\model5\\CRplantsresult.csv", row.names = FALSE)

########################################EXanimals################################################################################
################################################################################################################################
data <- read.csv("D:\\other study\\濒危物种研究\\物种发现数预测\\model5\\EXanimalsmodel.csv")


library(nlme)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(broom)
library(minpack.lm)  


start.list <- list(A = max(data$seq.cum), x0 = mean(data$Year), B = sd(data$Year))  # 初始值列表
model <- nlsLM(seq.cum ~ A / (1 + exp(-(Year - x0)/B)), data = data, start = start.list)  # nlsLM函数


bic <- BIC(model)


print(paste("BIC: ", bic))



new_data <- data.frame(Year = c(1993:2022, 2023:2042))
new_data$predicted_seq_cum <- predict(model, newdata = new_data)
new_data$predicted_seq_cum <- round(new_data$predicted_seq_cum,0)


observed <- data$seq.cum
predicted <- predict(model)
fit <- cor(observed, predicted)^2


model_eq <- paste("seq.cum =", round(coef(model)[1], 3), "/ (1 + exp(-(", round(coef(model)[2], 3), "* (Year - ", round(coef(model)[3], 3), "))))")


ggplot(data, aes(x = Year, y = seq.cum)) +
  geom_point(color = "blue") +
  geom_line(data = new_data, aes(x = Year, y = predicted_seq_cum), color = "red") +
  annotate("text", x = 2025, y = 2500, label = paste("R^2 =", round(fit, 3)), color = "black", hjust = 4, vjust = 1) +
  annotate("text", x = 2025, y = 2400, label = model_eq, color = "black", hjust = 0.8, vjust = 1) +
  labs(title = "EX animals", x = "Year", y = "seq.cum")

write.csv(new_data, file = "D:\\other study\\model5\\EXanimalsN.csv", row.names = FALSE) 


############################################################################################################################################
# ModelI-a
#sp.per~(a+b*sp.cum)*(3797-sp.cum
# Input data


model1 <- gnls(sp.per ~ (a + b * sp.cum) * (684 - sp.cum),
               data = data,
               start = list(a = 1e-02, b = 1e-04),
               weights = varPower(),
               verbose = TRUE,
               control = gnlsControl(returnObject = TRUE, minScale = 1e-500,
                                     tolerance = 0.001, nlsMaxIter = 3))


parameters <- coef(model1)
confidence_intervals <- confint(model1)
summary(model1)


yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[26+i])
  sp.cum.predict <- predict(model1, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model1)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval1 <- data.frame(
  Year = 2022:(2021 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval1)

################################################################################################################################################
#ModelI-b
#k(t)=a+b×S+c×S^2

model2=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2)*(684-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))


parameters <- coef(model2)
confidence_intervals <- confint(model2)
summary(model2)


yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[26+i])
  sp.cum.predict <- predict(model2, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model2)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval2 <- data.frame(
  Year = 2022:(2021 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval2)

############################################################################################################################################
# ModelI-c

model3 <- gnls(sp.per ~ (a + b * sp.cum + c * log(seq.cum)) * (684 - sp.cum),
               data = data,
               start = list(a = 1e-02, b = 1e-04, c = 1e-04),
               weights = varPower(),
               verbose = TRUE,
               control = gnlsControl(returnObject = TRUE, minScale = 1e-500, tolerance = 0.001, nlsMaxIter = 3))

parameters <- coef(model3)
confidence_intervals <- confint(model3)
summary(model3)

# Preparing the initial value of sp.cum.t
sp.cum.t <- tail(data$sp.cum, 1)


yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[26+i])
  sp.cum.predict <- predict(model3, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model3)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval3 <- data.frame(
  Year = 2022:(2021 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval3)
################################################################################################
#ModelI-d
#k(t)=a+b×S+c×log(N)+d×S×log(N)

model4=gnls(sp.per~(a+b*sp.cum+c*log(seq.cum)+d*sp.cum*log(seq.cum))*(684-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))




parameters <- coef(model4)
confidence_intervals <- confint(model4)
summary(model4)


yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[26+i])
  sp.cum.predict <- predict(model4, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model4)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval4 <- data.frame(
  Year = 2022:(2021 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval4)

####################################################################################################################

# ModelI-e
#k(t)=a+b×S+c×S^2+d×log(N)
# Parameter estimation 

model5=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2+d*log(seq.cum))*(684-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))




parameters <- coef(model5)
confidence_intervals <- confint(model5)
summary(model5)


yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[26+i])
  sp.cum.predict <- predict(model5, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model5)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval5 <- data.frame(
  Year = 2022:(2021 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval5)

###########################################################################################################################

# ModelI-f
#k(t)=a+b×S+c×log(N)+d×log(N)^2
# Parameter estimation 
model6=gnls(sp.per~(a+b*sp.cum+c*log(seq.cum)+d*log(seq.cum)^2)*(684-sp.cum),
             data=data,
             start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
             weights=varPower(),verbose=T,
             control=gnlsControl(returnObject=T,minScale=1e-500,
                                 tolerance=0.001,nlsMaxIter=3))



parameters <- coef(model6)
confidence_intervals <- confint(model6)
summary(model6)


yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[26+i])
  sp.cum.predict <- predict(model6, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model6)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval6 <- data.frame(
  Year = 2022:(2021 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval6)

##########################################################################################################################

# ModelI-g
#k(t)=a+b×S+c×S^2+d×log(N)+e×log(N)^2
# Parameter estimation 


model7=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2+d*log(seq.cum)+e*log(seq.cum)^2)*(684-sp.cum),
             data=data,
             start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01, e=0.01),
             weights=varPower(),verbose=T,
             control=gnlsControl(returnObject=T,minScale=1e-500,
                                 tolerance=0.001,nlsMaxIter=3))



parameters <- coef(model7)
confidence_intervals <- confint(model7)
summary(model7)



yr <- 21
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[26+i])
  sp.cum.predict <- predict(model7, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model7)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval7 <- data.frame(
  Year = 2022:(2021 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval7)

#######################################################################################################################

combined_results <- cbind(prediction_interval1, prediction_interval2, prediction_interval3, prediction_interval4, prediction_interval5, prediction_interval6,prediction_interval7,prediction_interval8,prediction_interval9,prediction_interval10,prediction_interval11,prediction_interval12)
colnames(combined_results) <- c("Year", "Model1", "Model2", "Model3", "Model4", "Model5","Model6","Model7")


write.csv(combined_results, file = "D:\\other study\\model5\\EXanimalsresult.csv", row.names = FALSE)

########################################EXplants################################################################################
################################################################################################################################
data <- read.csv("D:\\other study\\model5\\EXplantsmodel.csv")

library(nlme)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(broom)
library(minpack.lm)  


start.list <- list(A = max(data$seq.cum), x0 = mean(data$Year), B = sd(data$Year))  
model <- nlsLM(seq.cum ~ A / (1 + exp(-(Year - x0)/B)), data = data, start = start.list)  


bic <- BIC(model)


print(paste("BIC: ", bic))



new_data <- data.frame(Year = c(1993:2022, 2023:2042))
new_data$predicted_seq_cum <- predict(model, newdata = new_data)
new_data$predicted_seq_cum <- round(new_data$predicted_seq_cum,0)


observed <- data$seq.cum
predicted <- predict(model)
fit <- cor(observed, predicted)^2


model_eq <- paste("seq.cum =", round(coef(model)[1], 3), "/ (1 + exp(-(", round(coef(model)[2], 3), "* (Year - ", round(coef(model)[3], 3), "))))")


ggplot(data, aes(x = Year, y = seq.cum)) +
  geom_point(color = "blue") +
  geom_line(data = new_data, aes(x = Year, y = predicted_seq_cum), color = "red") +
  annotate("text", x = 2025, y = 1000, label = paste("R^2 =", round(fit, 3)), color = "black", hjust = 4, vjust = 1) +
  annotate("text", x = 2025, y = 900, label = model_eq, color = "black", hjust = 0.8, vjust = -1) +
  labs(title = "EX plants", x = "Year", y = "seq.cum")

write.csv(new_data, file = "D:\\other study\\model5\\EXplantsN2.csv", row.names = FALSE) 



############################################################################################################################################
# ModelI-a
#sp.per~(a+b*sp.cum)*(3797-sp.cum
# Input data

model1 <- gnls(sp.per ~ (a + b * sp.cum) * (163 - sp.cum),
               data = data,
               start = list(a = 1e-02, b = 1e-04),
               weights = varPower(),
               verbose = TRUE,
               control = gnlsControl(returnObject = TRUE, minScale = 1e-500,
                                     tolerance = 0.001, nlsMaxIter = 3))


parameters <- coef(model1)
confidence_intervals <- confint(model1)
summary(model1)


yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[21+i])
  sp.cum.predict <- predict(model1, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model1)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval1 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval1)

################################################################################################################################################
#ModelI-b
#k(t)=a+b×S+c×S^2

model2=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2)*(163-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))


parameters <- coef(model2)
confidence_intervals <- confint(model2)
summary(model2)



yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[21+i])
  sp.cum.predict <- predict(model2, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model2)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval2 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval2)

############################################################################################################################################
# ModelI-c

model3 <- gnls(sp.per ~ (a + b * sp.cum + c * log(seq.cum)) * (163 - sp.cum),
               data = data,
               start = list(a = 1e-02, b = 1e-04, c = 1e-04),
               weights = varPower(),
               verbose = TRUE,
               control = gnlsControl(returnObject = TRUE, minScale = 1e-500, tolerance = 0.001, nlsMaxIter = 3))

parameters <- coef(model3)
confidence_intervals <- confint(model3)
summary(model3)

# Preparing the initial value of sp.cum.t
sp.cum.t <- tail(data$sp.cum, 1)


yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[21+i])
  sp.cum.predict <- predict(model3, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model3)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval3 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval3)
################################################################################################
#ModelI-d
#k(t)=a+b×S+c×log(N)+d×S×log(N)

model4=gnls(sp.per~(a+b*sp.cum+c*log(seq.cum)+d*sp.cum*log(seq.cum))*(163-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))



parameters <- coef(model4)
confidence_intervals <- confint(model4)
summary(model4)


yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[21+i])
  sp.cum.predict <- predict(mode4, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model4)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval4 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval6)

####################################################################################################################

# ModelI-e
#k(t)=a+b×S+c×S^2+d×log(N)
# Parameter estimation 

model5=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2+d*log(seq.cum))*(163-sp.cum),
            data=data,
            start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
            weights=varPower(),verbose=T,
            control=gnlsControl(returnObject=T,minScale=1e-500,
                                tolerance=0.001,nlsMaxIter=3))



parameters <- coef(model5)
confidence_intervals <- confint(model5)
summary(model5)


yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[21+i])
  sp.cum.predict <- predict(model5, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model5)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval5 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval5)

###########################################################################################################################

# ModelI-f
#k(t)=a+b×S+c×log(N)+d×log(N)^2
# Parameter estimation 


model6=gnls(sp.per~(a+b*sp.cum+c*log(seq.cum)+d*log(seq.cum)^2)*(163-sp.cum),
             data=data,
             start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01),
             weights=varPower(),verbose=T,
             control=gnlsControl(returnObject=T,minScale=1e-500,
                                 tolerance=0.001,nlsMaxIter=3))




parameters <- coef(model6)
confidence_intervals <- confint(model6)
summary(model6)


yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[21+i])
  sp.cum.predict <- predict(model6, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model6)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval6 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval6)

##########################################################################################################################

# ModelI-g
#k(t)=a+b×S+c×S^2+d×log(N)+e×log(N)^2
# Parameter estimation 


model7=gnls(sp.per~(a+b*sp.cum+c*sp.cum^2+d*log(seq.cum)+e*log(seq.cum)^2)*(163-sp.cum),
             data=data,
             start=list(a=1e-02,b=1e-04,c=1e-04,d=0.01, e=0.01),
             weights=varPower(),verbose=T,
             control=gnlsControl(returnObject=T,minScale=1e-500,
                                 tolerance=0.001,nlsMaxIter=3))



parameters <- coef(model7)
confidence_intervals <- confint(model7)
summary(model7)



yr <- 20
sp.cum.predict.final <- NULL
lower_quantile <- numeric(yr)
upper_quantile <- numeric(yr)
sp.cum.t <- tail(data$sp.cum, 1)

for (i in 1:yr) {
  x_pred <- data.frame(sp.cum = sp.cum.t, seq.cum = new_data$predicted_seq_cum[21+i])
  sp.cum.predict <- predict(model7, newdata = x_pred)
  sp.cum.t <- sp.cum.predict + sp.cum.t  # Add the prediction to the cumulative sum
  
  # Calculate confidence interval manually
  se <- sigma(model7)
  t_value <- qt(0.975, df = nrow(data) - length(parameters))  # T-distribution quantile for 95% confidence interval
  margin_error <- t_value * se
  
  lower_quantile[i] <- sp.cum.t - margin_error
  upper_quantile[i] <- sp.cum.t + margin_error
  sp.cum.predict.final <- c(sp.cum.predict.final,sp.cum.t)
}

prediction_interval7 <- data.frame(
  Year = 2023:(2022 + yr),
  Prediction = sp.cum.predict.final,
  Lower_Quantile = lower_quantile,
  Upper_Quantile = upper_quantile
)

print(prediction_interval7)

#######################################################################################################################

combined_results <- cbind(prediction_interval1, prediction_interval2, prediction_interval3, prediction_interval4, prediction_interval5, prediction_interval6,prediction_interval7,prediction_interval8,prediction_interval9,prediction_interval10,prediction_interval11,prediction_interval12)
colnames(combined_results) <- c("Year", "Model1", "Model2", "Model3", "Model4", "Model5","Model6","Model7")

write.csv(combined_results, file = "D:\\other study\\濒危物种研究\\物种发现数预测\\model5\\EXplantsresult.csv", row.names = FALSE)
