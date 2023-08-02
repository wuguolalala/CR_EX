#########################################################################
#CRanimals
###########################################################################
library(ggplot2)
library(forecast) 
library(tseries) 
library(dplyr)
library(readxl)
library(tidyr)
library(ggfortify)
library(gridExtra)
library(nlme)
library(tidyverse)
library(broom)
library(minpack.lm) 

data <- read.csv("D:\\other study\\model5\\CRanimalsmodel.csv")  


sp_cum <- ts(data$sp.cum,frequency = 1, start= 1993)
seq_cum <- ts(data$seq.cum,frequency = 1, start= 1993)


p1 <- autoplot(sp_cum) + labs(x = "Year", y = "sp_cum")
p2 <- autoplot(seq_cum) + labs(x = "Year", y = "seq_cum")
gridExtra::grid.arrange(p1, p2, nrow = 2)



adf.test(sp_cum) 
adf.test(seq_cum) 



arimax <- auto.arima(y=sp_cum,xreg = seq_cum)
arimax
ARIMAXmod <- Arima(sp_cum,order = c(0,2,1),xreg = seq_cum)
summary(ARIMAXmod)

价
residuals <- residuals(ARIMAXmod)
qqnorm(residuals)
qqline(residuals)
Box.test(residuals, type = "Ljung-Box")


start.list <- list(A = max(data$seq.cum), x0 = mean(data$Year), B = sd(data$Year))  # 初始值列表
model <- nlsLM(seq.cum ~ A / (1 + exp(-(Year - x0)/B)), data = data, start = start.list)  # nlsLM函数


bic <- BIC(model)


print(paste("BIC: ", bic))



new_data <- data.frame(Year = c(2023:2042))
new_data$predicted_seq_cum <- predict(model, newdata = new_data)
new_data$predicted_seq_cum <- round(new_data$predicted_seq_cum,0)


future_seq_cum <- ts(new_data$predicted_seq_cum, frequency = 1, start = 2023, end = 2042)  
pred <- forecast(ARIMAXmod, h = 20, xreg = future_seq_cum)
print(pred)


plot(pred, xlab = "Year", ylab = "sp_cum")
lines(sp_cum, col = "black")
write.csv(pred, file = "D:\\other study\\CRanimals.csv", row.names = FALSE) 


write.csv(pred, "D://other study//CRanimals.csv", row.names = TRUE)


#########################################################################
#CRplants
###########################################################################
library(ggplot2)
library(forecast) 
library(tseries) 
library(dplyr)
library(readxl)
library(tidyr)
library(ggfortify)
library(gridExtra)
library(nlme)
library(tidyverse)
library(broom)
library(minpack.lm) 

data <- read.csv("D:\\other study\\model5\\CRplantsmodel.csv")  


sp_cum <- ts(data$sp.cum,frequency = 1, start= 1994)
seq_cum <- ts(data$seq.cum,frequency = 1, start= 1994)


p1 <- autoplot(sp_cum) + labs(x = "Year", y = "sp_cum")
p2 <- autoplot(seq_cum) + labs(x = "Year", y = "seq_cum")
gridExtra::grid.arrange(p1, p2, nrow = 2)



adf.test(sp_cum)  
adf.test(seq_cum) 



arimax <- auto.arima(y=sp_cum,xreg = seq_cum)
arimax
ARIMAXmod <- Arima(sp_cum,order = c(0,2,2),xreg = seq_cum)
summary(ARIMAXmod)


residuals <- residuals(ARIMAXmod)
qqnorm(residuals)
qqline(residuals)
Box.test(residuals, type = "Ljung-Box")


start.list <- list(A = max(data$seq.cum), x0 = mean(data$Year), B = sd(data$Year)) 
model <- nlsLM(seq.cum ~ A / (1 + exp(-(Year - x0)/B)), data = data, start = start.list) 


bic <- BIC(model)


print(paste("BIC: ", bic))



new_data <- data.frame(Year = c(2023:2042))
new_data$predicted_seq_cum <- predict(model, newdata = new_data)
new_data$predicted_seq_cum <- round(new_data$predicted_seq_cum,0)


future_seq_cum <- ts(new_data$predicted_seq_cum, frequency = 1, start = 2023, end = 2042)  # 未来外部变量
pred <- forecast(ARIMAXmod, h = 20, xreg = future_seq_cum)
print(pred)


plot(pred, xlab = "Year", ylab = "sp_cum")
lines(sp_cum, col = "black")



write.csv(pred, "D://other study//CRplants.csv", row.names = TRUE)


#########################################################################
#EXanimals
###########################################################################
library(ggplot2)
library(forecast) 
library(tseries) 
library(dplyr)
library(readxl)
library(tidyr)
library(ggfortify)
library(gridExtra)
library(nlme)
library(tidyverse)
library(broom)
library(minpack.lm) 

data <- read.csv("D:\\other study\\model5\\EXanimalsmodel.csv")  


sp_cum <- ts(data$sp.cum,frequency = 1, start= 1993)
seq_cum <- ts(data$seq.cum,frequency = 1, start= 1993)


p1 <- autoplot(sp_cum) + labs(x = "Year", y = "sp_cum")
p2 <- autoplot(seq_cum) + labs(x = "Year", y = "seq_cum")
gridExtra::grid.arrange(p1, p2, nrow = 2)


性
adf.test(sp_cum)  
adf.test(seq_cum)  



arimax <- auto.arima(y=sp_cum,xreg = seq_cum)
arimax
ARIMAXmod <- Arima(sp_cum,order = c(0,1,0),xreg = seq_cum)
summary(ARIMAXmod)


residuals <- residuals(ARIMAXmod)
qqnorm(residuals)
qqline(residuals)
Box.test(residuals, type = "Ljung-Box")


start.list <- list(A = max(data$seq.cum), x0 = mean(data$Year), B = sd(data$Year))  # 初始值列表
model <- nlsLM(seq.cum ~ A / (1 + exp(-(Year - x0)/B)), data = data, start = start.list)  # nlsLM函数


bic <- BIC(model)


print(paste("BIC: ", bic))



new_data <- data.frame(Year = c(2021:2042))
new_data$predicted_seq_cum <- predict(model, newdata = new_data)
new_data$predicted_seq_cum <- round(new_data$predicted_seq_cum,0)


future_seq_cum <- ts(new_data$predicted_seq_cum, frequency = 1, start = 2021, end = 2042)  # 未来外部变量
pred <- forecast(ARIMAXmod, h = 20, xreg = future_seq_cum)
print(pred)


plot(pred, xlab = "Year", ylab = "sp_cum")
lines(sp_cum, col = "black")



write.csv(pred, "D://other study//EXanimals.csv", row.names = TRUE)

#########################################################################
#EXplants
###########################################################################
library(ggplot2)
library(forecast) 
library(tseries) 
library(dplyr)
library(readxl)
library(tidyr)
library(ggfortify)
library(gridExtra)
library(nlme)
library(tidyverse)
library(broom)
library(minpack.lm) 

data <- read.csv("D:\\other study\\model5\\EXplantsmodel.csv") 


sp_cum <- ts(data$sp.cum,frequency = 1, start= 1)
seq_cum <- ts(data$seq.cum,frequency = 1, start= 1)


p1 <- autoplot(sp_cum) + labs(x = "Year", y = "sp_cum")
p2 <- autoplot(seq_cum) + labs(x = "Year", y = "seq_cum")
gridExtra::grid.arrange(p1, p2, nrow = 2)



adf.test(sp_cum) 
adf.test(seq_cum) 



arimax <- auto.arima(y=sp_cum,xreg = seq_cum)
arimax
ARIMAXmod <- Arima(sp_cum,order = c(2,2,0),xreg = seq_cum)
summary(ARIMAXmod)


residuals <- residuals(ARIMAXmod)
qqnorm(residuals)
qqline(residuals)
Box.test(residuals, type = "Ljung-Box")


start.list <- list(A = max(data$seq.cum), x0 = mean(data$Year), B = sd(data$Year)) 
model <- nlsLM(seq.cum ~ A / (1 + exp(-(Year - x0)/B)), data = data, start = start.list)  


bic <- BIC(model)


print(paste("BIC: ", bic))



new_data <- data.frame(Year = c(2023:2042))
new_data$predicted_seq_cum <- predict(model, newdata = new_data)
new_data$predicted_seq_cum <- round(new_data$predicted_seq_cum,0)


future_seq_cum <- ts(new_data$predicted_seq_cum, frequency = 1, start = 2023, end = 2042) 
pred <- forecast(ARIMAXmod, h = 20, xreg = future_seq_cum)
print(pred)


plot(pred, xlab = "Year", ylab = "sp_cum")
lines(sp_cum, col = "black")


write.csv(pred, "D://other study//EXplants.csv", row.names = TRUE)

