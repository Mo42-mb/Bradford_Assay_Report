### Bradford assay ###
library(readxl)
library(knitr)
library(dplyr)

###### Load palettes and data ######
# UPDATE working directory and data files. Knit/run the rest.
setwd("/Users")
dat_full <- read_excel("example/example_samples_Bradford.xlsx", sheet = "standards")
dat <- dat_full %>%
  select(protein_concentration_mg_per_mL, absorbance_corrected)
samples_full <- read_excel("example/example_samples_Bradford.xlsx", sheet = "samples")
# Make concise & exclude multiplication information until the end
samples <- samples_full %>%
  select(sample, absorbance_corrected)

##### first: LINEAR regression #####
linear_model_whole <- lm(protein_concentration_mg_per_mL ~ absorbance_corrected, data = dat)

#linear model should only be from absorbance range 0.2 to 1.0
linear_model <- lm(protein_concentration_mg_per_mL ~ absorbance_corrected, data = subset(dat, protein_concentration_mg_per_mL >= 0.2 & protein_concentration_mg_per_mL <= 1.0))
linear_model
summary(linear_model)

# save equation variables
# y = mx + b
b <- coef(linear_model)[1]
m <- coef(linear_model)[2]
R2adj_linear <- summary(linear_model)$adj.r.squared

# linear function
predict_y <- function(x) {
  y <- m * x + b
  return(y)
}

## CALCULATION
(linear_predicted_y <- predict_y(samples$absorbance_corrected))

#save data frame
df_predicted_l <- cbind(samples, linear_predicted_y)
colnames(df_predicted_l)[2:3] <- c("absorbance_corrected", "protein_concentration_mg_per_mL")

print(paste("Model: y = ", round(m, digits=3), "*x +", round(b, digits=3)))

# Basic plot in Base R
plot(protein_concentration_mg_per_mL ~ absorbance_corrected, data = dat)
points(df_predicted_l[,2], df_predicted_l[,3], col = "black", pch = 16, cex = 1.5)
abline(linear_model, col = "cyan4")
abline(linear_model_whole, col = "orange")
abline(h = 0.2, col = "cyan2", lty = 2)
abline(h = 1, col = "cyan2", lty = 2)
title("Linear Model")
text(0.2,0, labels = paste("R^2 (linear) = ", round(R2adj_linear, digits = 3)))
legend("topleft", legend = c("Model: Linear dataset", "Model: All data", "~Linear Region"),
       col = c("cyan4", "orange", "cyan2"), lwd = 1.5, lty = c(1,1,2))


##### second: quadratic regression #####
quad_model <- lm(protein_concentration_mg_per_mL ~ absorbance_corrected + I(absorbance_corrected^2), data = dat)
quad_model
summary(quad_model)

# save equation variables
# y = ax2 + bx + c
c <- coef(quad_model)[1]
b <- coef(quad_model)[2]
a <- coef(quad_model)[3]

R2adj_quad <- summary(quad_model)$adj.r.squared

# linear function
quad_predict_y <- function(x) {
  y <- a * x^2 + b * x + c
  return(y)
}

#for plot
quad_fun <- function(x) {
  a * x^2 + b * x + c
}

## CALCULATION
(quad_predicted_y <- quad_predict_y(samples$absorbance_corrected))

#save data frame
df_predicted_q <- cbind(samples, quad_predicted_y)
colnames(df_predicted_q)[2:3] <- c("absorbance_corrected", "protein_concentration_mg_per_mL")

print(paste("Model: y = ", round(a, digits = 3), "*x^2 +", round(b, digits =3), "*x +", round(c, digits=3)))

# Basic plot in Base R
plot(protein_concentration_mg_per_mL ~ absorbance_corrected, data = dat)
points(df_predicted_q[,2], df_predicted_q[,3], col = "black", pch = 16, cex = 1.5)
curve(quad_fun, from = min(dat$absorbance_corrected), to = max(dat$absorbance_corrected), col = "purple", lwd = 2, add = TRUE)
title("Quadratic Model")
text(0.2,0, labels = paste("R^2 = ", round(R2adj_quad, digits = 3)))

##### Finally - compare models #####
plot(protein_concentration_mg_per_mL ~ absorbance_corrected, data = dat)
points(df_predicted_l[,2], df_predicted_l[,3], col = "cyan4", pch = 16, cex = 1.5)
points(df_predicted_q[,2], df_predicted_q[,3], col = "purple", pch = 16, cex = 1.5)
abline(linear_model, col = "cyan4", lwd = 2)
curve(quad_fun, from = min(dat$absorbance_corrected), to = max(dat$absorbance_corrected), col = "purple", lwd = 2, add = TRUE)
title("Compare Linear and Quadratic Models")
legend("topleft", legend = c("Linear Model", "Quadratic Model"),
       col = c("cyan4", "purple"), lwd = 2)

# print summary
# print sample, linear-predicted value, quadratic-predicted value, and R2 from both equations
summary <- cbind(samples_full, df_predicted_l[,3], df_predicted_q[,3])
summary$linear_prediction <- summary$multiplication_factor*summary$`df_predicted_l[, 3]`
summary$quadratic_prediction <- summary$multiplication_factor*summary$`df_predicted_q[, 3]`

table_summary <- summary %>% 
  select(sample, linear_prediction, quadratic_prediction)

table_summary[,2] <- round(table_summary[,2], digits = 2)
table_summary[,3] <- round(table_summary[,3], digits = 2)

kable(table_summary, caption = "Protein concentration (mg/mL)")
print(paste("Linear Model R2 adjusted =", round(R2adj_linear, digits = 3)))
print(paste("Quadratic Model R2 adjusted =", round(R2adj_quad, digits = 3)))

# save as markdown
knitr::spin("Bradford_analysis.R", knit = FALSE)



