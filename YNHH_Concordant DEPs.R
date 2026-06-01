##########################################################################################################################################
#ALL SAMPLES ANALYSIS (YNHH cohort)
##########################################################################################################################################

#Load Packages
library(tidyverse)
library(plyr)
library(dplyr)
library(readxl)
library(writexl)
library(gplots)
library(ggrepel)
library(factoextra)

#Proteomics data file upload
excel_sheets("YNHH_LVO_data.xlsx")
YNHH_LVO_data= excel_sheets("YNHH_LVO_data.xlsx") %>% map(~read_xlsx("YNHH_LVO_data.xlsx",.))
YNHH_LVO_data

####################################################
#Visualize the raw data on a PCA plot (YNHH cohort)
####################################################

#PCA plot of raw LVO data
R.PCA<- prcomp(YNHH_LVO_data[[2]][,3:7303], scale=TRUE)
R.PCA
X11()
pca_plot= fviz_pca_ind(R.PCA, col.ind=YNHH_LVO_data[[2]]$Outcome, title= "PCA Plot of Raw LVO Data (YNHH)", addEllipses = FALSE, 
                       label= "none", pointsize= 4)
# Define colors for the groups
colors = c("LVO" = "red", "noLVO" = "blue")
# Add custom colors
pca_plot + scale_color_manual(values = colors)

#######################################################################################################################################

#####################################
#Log transformation of Raw YNHH Data
#####################################

#Select columns and log transform the data
dat_log <- YNHH_LVO_data[[1]] %>%
  dplyr::mutate(across(-Protein_ID, log2))

View(dat_log)

write.csv(dat_log, "Log_YNHH_data.csv")

#######################################################################################################################################

#################################################
# NORMALITY ASSESSMENT AFTER LOG2 TRANSFORMATION
#################################################

log_mat <- YNHH_LVO_data[[3]] %>%
  dplyr::select(-Protein_ID) %>%
  as.data.frame()

# Keep only numeric columns
log_mat <- log_mat[, sapply(log_mat, is.numeric)]

# Remove sample columns with too few values or zero variance
log_mat <- log_mat[, sapply(log_mat, function(x) sum(is.finite(x)) > 3)]
log_mat <- log_mat[, sapply(log_mat, function(x) sd(x, na.rm = TRUE) > 0)]

# GLOBAL VALUES

all_values <- unlist(log_mat)
all_values <- all_values[is.finite(all_values)]

hist_df <- data.frame(value = all_values)

# Calculate global mean and median
global_mean <- mean(all_values, na.rm = TRUE)
global_median <- median(all_values, na.rm = TRUE)

# Print values
cat("Global Mean:", global_mean, "\n")
cat("Global Median:", global_median, "\n")

# GLOBAL HISTOGRAM OF ALL LOG2 VALUES WITH MEAN AND MEDIAN

p_hist <- ggplot(hist_df, aes(x = value)) +
  geom_histogram(bins = 60, fill = "grey75", color = "black") +
  
  # Mean (dashed line)
  geom_vline(xintercept = global_mean, linetype = "dashed", size = 1) +
  
  # Median (dotted line)
  geom_vline(xintercept = global_median, linetype = "dotted", size = 1) +
  
  theme_classic() +
  labs(
    title = "Global Distribution of Log2-Transformed Intensities",
    subtitle = paste0(
      "Mean = ", round(global_mean, 2),
      " | Median = ", round(global_median, 2)
    ),
    x = "Log2 intensity",
    y = "Frequency"
  )

print(p_hist)

ggsave("Histogram_Log2_Global_Mean_Median_YNHH.png", p_hist,
       width = 8, height = 6, dpi = 300, bg = "white")

ggsave("Histogram_Log2_Global_Mean_Median_YNHH.pdf", p_hist,
       width = 8, height = 6)

#######################################################################################################################################

#########################################
# Distribution of unnormalized YNHH data
#########################################

unnorm_data <- read.csv("YNHH_Expression_Data.csv", row.names = 1)
meta_data <- read.csv("YNHH_Meta_Data.csv", row.names = 1)

# Create a vector of colors
colors <- c("red", "green", "blue", "purple", "orange", "yellow", "pink", "cyan")  # Add as many colors as needed

# Assign a unique color to each condition
conditionNames <- unique(meta_data$Outcome)
colorVector <- setNames(colors[1:length(conditionNames)], conditionNames)

# Create a color vector for the samples
statusCol <- colorVector[meta_data$Outcome]

# Set margin sizes (bottom, left, top, right)
par(mar = c(5, 4, 4, 2) + 0.1)

# Check distributions of samples using boxplots
boxplot(unnorm_data, 
        xlab="", 
        ylab="Expression levels",
        las=2,
        col=statusCol,
        main="Box plot of unnormalized data (YNHH)")

# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(as.matrix(unnorm_data)), col="blue")

############################################################################################################################################

###########################################################################
#YNHH dataset- Differential abundance analysis between LVO and non-LVO AIS
###########################################################################

#Visualizing the Data
dat = YNHH_LVO_data[[3]]
View(dat)

#Creating a T-test function for multiple experiments
t_test <- function(dt,grp1,grp2){
  # Subset Total Stroke Case group and convert to numeric
  x <- dt[grp1] %>% unlist %>% as.numeric()
  # Subset Healthy Control group and convert to numeric
  y <- dt[grp2] %>% unlist %>% as.numeric()
  # Perform t-test using the mean of x and y
  result <- t.test(x, y)
  # Extract p-values from the results
  p_vals <- tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
} 

#Apply t-test function to data using plyr adply
#.margins = 1, slice by rows, .fun = t_test plus t_test arguments
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:35), grp2 = c(36:65)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:35]), as.numeric(dat[1,36:65]))$p.value

#Bind columns to create transformed data frame
dat_combine = bind_cols(dat, dat_pvals[,66])
View (dat_combine)

#Calculating log-fold change
dat_fc = dat_combine %>% 
  group_by(Protein_ID) %>% 
  dplyr::mutate(mean_LVO_case = mean(c(LVO1,	LVO2,	LVO3,	LVO4,	LVO5,	LVO6,	LVO7,	LVO8,	LVO9,	LVO10,	LVO11,	LVO12,	LVO13,	LVO14,	
                                       LVO15,	LVO16,	LVO17,	LVO18, LVO19,	LVO20,	LVO21,	LVO22,	LVO23,	LVO24,	LVO25,	LVO26,	
                                       LVO27,	LVO28,	LVO29,	LVO30,	LVO31,	LVO32,	LVO33,	LVO34)),
                mean_noLVO_case= mean(c(noLVO1,	noLVO2,	noLVO3,	noLVO4,	noLVO5,	noLVO6,	noLVO7,	noLVO8,	noLVO9,	noLVO10,	noLVO11,	noLVO12,	
                                        noLVO13,	noLVO14,	noLVO15,	noLVO16,	noLVO17,	noLVO18,	noLVO19,	noLVO20,	noLVO21,	noLVO22, 
                                        noLVO23,	noLVO24,	noLVO25,	noLVO26,	noLVO27,	noLVO28,	noLVO29,	noLVO30)),
                log_fc = mean_LVO_case - mean_noLVO_case,
                log_pval = -1*log10(p_val))
View(dat_fc)

#Save final data with list of final data in csv file
write.csv(dat_fc, "Final_LVO_data_YNHH.csv")

#Volcano plot of log-fold change on x-axis and log p-value on y-axis
dat_fc %>% ggplot(aes(log_fc,log_pval)) + geom_point()

#Volcano plot at FC 1.2 and p<0.05

VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val))) + geom_point() + theme_minimal()
VP 
#Add vertical lines for Log2 FC and a horizontal line for p-value threshold
VP2= VP + geom_vline(xintercept = c(-0.26, 0.26), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Add a column of NAs
dat_fc$diffexpressed= "NO"
#Set Log2 FC and p-value cut-offs in the new column
dat_fc$diffexpressed[dat_fc$log_fc>0.26 & dat_fc$p_val<0.05] <- "UP"
dat_fc$diffexpressed[dat_fc$log_fc< -0.26 & dat_fc$p_val<0.05] <- "DOWN"
#Re-plot but this time color the points with "diffexpressed"
VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val), col= diffexpressed)) + geom_point() + theme_minimal()
VP
#Add lines as before..
VP2= VP + geom_vline(xintercept = c(-0.26, 0.26), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Change point colors
VP3= VP2 + scale_color_manual(values= c("blue", "black", "red"))
mycolors= c("blue", "red", "black")  
names(mycolors) = c("DOWN", "UP", "NO")
VP3= VP2 + scale_color_manual(values=mycolors)
#Create a new column "proteinlabel" that will contain names of differentially expressed protein IDs
dat_fc$proteinlabel= NA
dat_fc$proteinlabel[dat_fc$diffexpressed != "NO"] <- dat_fc$Protein_ID[dat_fc$diffexpressed != "NO"]
ggplot(data=dat_fc, aes(x= log_fc, y= -log10(p_val), col= diffexpressed, label=proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text()
View(dat_fc)

#Plot the Volcano plot using all layers used so far
X11()
ggplot(data= dat_fc, aes(x=log_fc, y= -log10(p_val), col= diffexpressed, label= proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline (xintercept = c(-0.26, 0.26), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")

#Proteins with significant observations
final_data<-dat_fc %>%
  #Filter for significant observations
  filter(log_pval >= 1.3 & (log_fc >= 0.26 | log_fc <= -0.26)) %>% 
  #Ungroup the data
  ungroup() %>% 
  #Select columns of interest
  select(Protein_ID, LVO1:LVO34, noLVO1:noLVO30, mean_LVO_case, mean_noLVO_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in csv file
write.csv(final_data, "Final_diffproteins_LVO_YNHH.csv")

############################################################################################################################################

###############################################################################################
# Bivariable Logistic Regression analysis of Concordant DEPs across two datasets (YNHH cohort)
###############################################################################################

library(readxl)
library(dplyr)
library(purrr)
library(openxlsx)

# ------------ USER INPUTS ------------
input_file  <- "YNHH_LVO_Data.xlsx"
sheet       <- 4

group_col   <- "Group"

meta_vars <- c("Age", "Sex", "DM", "HTN", "CAD", "AFIB", "Smoking", "Prior_Stroke", "CE", "EVT", "NIHSS0", "GCS")

protein_start_col <- 15
# -------------------------------------

df <- read_excel(input_file, sheet = sheet)

# Ensure binary coding: controls = noLVO (0), cases = LVO (1)
df[[group_col]] <- factor(df[[group_col]], levels = c("noLVO", "LVO"))

protein_cols <- names(df)[protein_start_col:ncol(df)]
protein_cols

# ---------- MODEL FUNCTION ----------
run_model <- function(data, protein, covariate, outcome = "Group") {
  
  tmp <- data %>%
    select(all_of(c(outcome, protein, covariate))) %>%
    na.omit()
  
  if (length(unique(tmp[[outcome]])) < 2)
    return(c(OR=NA, CI_low=NA, CI_high=NA, p=NA))
  
  # Wrap variable names in backticks for safe formula
  formula_text <- paste0("`", outcome, "` ~ `", covariate, "` + `", protein, "`")
  fit <- glm(
    formula = as.formula(formula_text),
    data = tmp,
    family = binomial()
  )
  
  s <- summary(fit)$coefficients
  protein_row <- paste0("`", protein, "`")
  
  if (!(protein_row %in% rownames(s)))
    return(c(OR=NA, CI_low=NA, CI_high=NA, p=NA))
  
  est <- s[protein_row, "Estimate"]
  se  <- s[protein_row, "Std. Error"]
  p   <- s[protein_row, "Pr(>|z|)"]
  
  OR      <- exp(est)
  CI_low  <- exp(est - 1.96 * se)
  CI_high <- exp(est + 1.96 * se)
  
  c(OR=OR, CI_low=CI_low, CI_high=CI_high, p=p)
}

# -------- RUN MODELS ACROSS ALL PROTEINS × COVARIATES --------

results <- map_df(protein_cols, function(prot) {
  map_df(meta_vars, function(meta) {
    out <- run_model(df, protein = prot, covariate = meta)
    
    dplyr::tibble(
      Protein   = prot,
      Covariate = meta,
      OR        = out["OR"],
      CI_low    = out["CI_low"],
      CI_high   = out["CI_high"],
      p_value   = out["p"]
    )
  })
})

# -------- SAVE AS CSV --------
write.csv(results, "Bivariable_logistic_results_GMH.csv", row.names = FALSE)

###########################################################################################################################################

#################################################
# Prediction models to classify LVO from non-LVO
#################################################

library(dplyr)
library(purrr)
library(pROC)
library(tidyr)
library(readxl)

# ------------ USER INPUTS ------------
input_file  <- "YNHH_LVO_Data.xlsx"
sheet       <- 4

group_col   <- "Group"
protein_start_col <- 15

# ------------ READ DATA ------------
df <- read_excel(input_file, sheet = sheet)

# ------------ Extract concordant DAPs ------------
concordant_DAPs <- colnames(df)[protein_start_col:ncol(df)]
print(concordant_DAPs)

# ------------ Subset model data ------------
model_data <- df %>%
  select(all_of(c(group_col, "NIHSS0", concordant_DAPs)))

# Ensure correct factor levels (controls first, cases second)
model_data$Group <- factor(model_data$Group, levels = c("noLVO","LVO"))

# ------------ Helper function for ROC + AUC ------------
run_model_auc <- function(formula_text, data) {
  
  vars <- all.vars(as.formula(formula_text))
  tmp <- data %>% select(all_of(vars)) %>% na.omit()
  
  if(length(unique(tmp$Group)) < 2){
    warning(paste("Skipping model", formula_text, "- only one class present"))
    return(list(model = NULL, roc = NULL, auc = NA, ci = c(NA,NA,NA)))
  }
  
  fit <- glm(as.formula(formula_text), data = tmp, family = binomial)
  preds <- predict(fit, newdata = tmp, type = "response")
  
  roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"), direction = "<")
  auc_val <- auc(roc_obj)
  ci_val  <- ci.auc(roc_obj)
  
  list(model = fit, roc = roc_obj, auc = auc_val, ci = ci_val)
}

# ------------ Helper function for PDF ROC plot only ------------
plot_roc_stata_clean <- function(roc_obj, auc_val, ci_val, title, filename_base){
  
  pdf(paste0(filename_base, ".pdf"), width = 6, height = 6)
  
  plot(roc_obj,
       col = "blue",
       lwd = 2,
       legacy.axes = TRUE,
       main = title,
       print.auc = FALSE)
  
  abline(a = 0, b = 1, lty = 2, col = "grey")
  
  legend("bottomright", 
         legend = paste0(
           "AUC = ", round(auc_val, 3),
           " (95% CI: ", round(ci_val[1], 3), "-", round(ci_val[3], 3), ")"
         ),
         bty = "n")
  
  dev.off()
}

# ------------ MODEL 1: NIHSS0 only ------------
model1 <- run_model_auc("Group ~ NIHSS0", model_data)

# ------------ MODEL 2: individual proteins ------------
model2_list <- map(concordant_DAPs, function(prot) {
  run_model_auc(paste0("Group ~ `", prot, "`"), model_data)
})
names(model2_list) <- concordant_DAPs

# ------------ MODEL 3: protein + NIHSS0 ------------
model3_list <- map(concordant_DAPs, function(prot) {
  run_model_auc(paste0("Group ~ NIHSS0 + `", prot, "`"), model_data)
})
names(model3_list) <- concordant_DAPs

# ------------------- Extract significant concordant DAPs from Model 3 ------------------- 
sig_proteins_m3 <- keep(names(model3_list), function(prot) { 
  fit <- model3_list[[prot]]$model
  if(is.null(fit)) return(FALSE)  # skip if model didn't run
  
  coefs <- summary(fit)$coefficients
  # Match protein term using grep to handle backticks
  prot_row <- grep(prot, rownames(coefs), fixed = TRUE)
  if(length(prot_row) == 0) return(FALSE)  # protein not found
  
  p_val <- coefs[prot_row, "Pr(>|z|)"]
  !is.na(p_val) && p_val < 0.05
})

sig_proteins_m3

# ------------------- 4) Concordant DAP pairs + NIHSS0 (filtered from Model 3) ------------------- 
model4_list <- list()

if(length(sig_proteins_m3) >= 2){ 
  combos <- combn(sig_proteins_m3, 2, simplify = FALSE) 
  
  for(combo in combos){ 
    # wrap protein names in backticks
    combo_backtick <- paste0("`", combo, "`")
    formula4 <- paste("Group ~ NIHSS0 +", paste(combo_backtick, collapse = " + "))
    
    vars <- c("Group", "NIHSS0", combo) 
    tmp <- model_data %>% select(all_of(vars)) %>% na.omit() 
    
    if(length(unique(tmp$Group)) < 2) next 
    
    fit <- glm(as.formula(formula4), data = tmp, family = binomial) 
    coefs <- summary(fit)$coefficients 
    
    # Extract p-values using grep to match backtick names
    protein_rows <- sapply(combo, function(p) grep(p, rownames(coefs), fixed = TRUE))
    protein_p <- coefs[protein_rows, "Pr(>|z|)"]
    
    if(all(!is.na(protein_p)) && all(protein_p < 0.05)){ 
      preds <- predict(fit, newdata = tmp, type = "response") 
      roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO")) 
      
      model4_list[[paste(combo, collapse = "_")]] <- list( 
        model = fit, 
        roc = roc_obj, 
        auc = auc(roc_obj), 
        ci = ci.auc(roc_obj) 
      ) 
    } 
  } 
}

model4_list
  
# ------------ SAVE AUC RESULTS ------------

# Model 1
auc_results <- tibble(
  Model = "NIHSS0",
  AUC = as.numeric(model1$auc),
  CI_lower = as.numeric(model1$ci[1]),
  CI_upper = as.numeric(model1$ci[3])
)

# Model 2
m2_df <- tibble(
  Model = paste0("Protein ", names(model2_list)),
  AUC = sapply(model2_list, function(x) as.numeric(x$auc)),
  CI_lower = sapply(model2_list, function(x) as.numeric(x$ci[1])),
  CI_upper = sapply(model2_list, function(x) as.numeric(x$ci[3]))
) %>% mutate(across(everything(), as.vector))

auc_results <- bind_rows(auc_results, m2_df)

# Model 3
m3_df <- tibble(
  Model = paste0(names(model3_list), " + NIHSS0"),
  AUC = sapply(model3_list, function(x) as.numeric(x$auc)),
  CI_lower = sapply(model3_list, function(x) as.numeric(x$ci[1])),
  CI_upper = sapply(model3_list, function(x) as.numeric(x$ci[3]))
) %>% mutate(across(everything(), as.vector))

auc_results <- bind_rows(auc_results, m3_df)

# Model 4: two proteins + NIHSS0 
if(length(model4_list) > 0){ 
  m4_df <- tibble( 
    Model = paste0(names(model4_list), " + NIHSS0"), 
    AUC = sapply(model4_list, function(x) if(!is.null(x$auc)) 
      as.numeric(x$auc) else NA_real_), 
    CI_lower = sapply(model4_list, function(x) if(!is.null(x$ci)) 
      as.numeric(x$ci[1]) else NA_real_), 
    CI_upper = sapply(model4_list, function(x) if(!is.null(x$ci)) 
      as.numeric(x$ci[3]) else NA_real_) 
  ) 
  m4_df <- m4_df %>% mutate(across(everything(), ~ as.vector(.))) 
  auc_results <- bind_rows(auc_results, m4_df) 
  }

# ------------ DeLong test: Model 1 vs all other models ------------

get_delong_p <- function(roc_ref, roc_comp){
  if(is.null(roc_ref) || is.null(roc_comp)) return(NA_real_)
  
  # Ensure both ROC objects have valid responses
  if(length(roc_ref$response) == 0 || length(roc_comp$response) == 0) return(NA_real_)
  
  # Try DeLong test safely
  out <- tryCatch({
    test <- roc.test(roc_ref, roc_comp, method = "delong")
    as.numeric(test$p.value)
  }, error = function(e) NA_real_)
  
  return(out)
}

# Initialize vector
delong_pvals <- c()

# Model 1 (reference) → p-value = NA
delong_pvals <- c(delong_pvals, NA_real_)

# Model 2
delong_pvals <- c(
  delong_pvals,
  sapply(model2_list, function(x) get_delong_p(model1$roc, x$roc))
)

# Model 3
delong_pvals <- c(
  delong_pvals,
  sapply(model3_list, function(x) get_delong_p(model1$roc, x$roc))
)

# Model 4
if(length(model4_list) > 0){
  delong_pvals <- c(
    delong_pvals,
    sapply(model4_list, function(x) get_delong_p(model1$roc, x$roc))
  )
}

# Add to results
auc_results$DeLong_p_value <- delong_pvals

# ============================================================
# Extract Odds Ratios (OR), 95% CI, and p-values
# for all models and save combined results
# ============================================================

# ------------ Helper function to extract OR results ------------
extract_or_results <- function(fit, model_name){
  
  if(is.null(fit)){
    return(NULL)
  }
  
  coefs <- summary(fit)$coefficients
  
  # Remove intercept
  coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
  
  # OR and 95% CI
  OR  <- exp(coefs[, "Estimate"])
  LCL <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
  UCL <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
  
  # p-values
  pval <- coefs[, "Pr(>|z|)"]
  
  tibble(
    Model = model_name,
    Variable = rownames(coefs),
    OR = OR,
    CI_lower = LCL,
    CI_upper = UCL,
    P_value = pval
  )
}

# ============================================================
# Extract OR results from all models
# ============================================================

or_results <- list()

# ------------ Model 1 ------------
or_results[["Model1"]] <- extract_or_results(
  model1$model,
  "NIHSS0"
)

# ------------ Model 2 ------------
m2_or <- map2(
  model2_list,
  names(model2_list),
  ~ extract_or_results(.x$model, paste0("Protein ", .y))
)

or_results[["Model2"]] <- bind_rows(m2_or)

# ------------ Model 3 ------------
m3_or <- map2(
  model3_list,
  names(model3_list),
  ~ extract_or_results(.x$model, paste0(.y, " + NIHSS0"))
)

or_results[["Model3"]] <- bind_rows(m3_or)

# ------------ Model 4 ------------
if(length(model4_list) > 0){
  
  m4_or <- map2(
    model4_list,
    names(model4_list),
    ~ extract_or_results(.x$model, paste0(.y, " + NIHSS0"))
  )
  
  or_results[["Model4"]] <- bind_rows(m4_or)
}

# ============================================================
# Combine all OR results
# ============================================================

final_or_results <- bind_rows(or_results)

# ============================================================
# Merge with AUC results
# ============================================================

final_combined_results <- final_or_results %>%
  left_join(
    auc_results,
    by = "Model"
  )

# ============================================================
# View results
# ============================================================

View(final_combined_results)

# ============================================================
# Save combined results to CSV
# ============================================================

write.csv(
  final_combined_results,
  "LVO_prediction_models_OR_AUC_YNHH.csv",
  row.names = FALSE
)

# ------------ PLOT ROC CURVES ------------

# Model 1
if(!is.null(model1$roc)){
  plot_roc_stata_clean(
    roc_obj = model1$roc,
    auc_val = as.numeric(model1$auc),
    ci_val  = as.numeric(model1$ci),
    title = "Model 1: NIHSS0 only",
    filename_base = "ROC_Model1_NIHSS0"
  )
}

# Model 2
for(prot in names(model2_list)){
  m2 <- model2_list[[prot]]
  if(!is.null(m2$roc)){
    plot_roc_stata_clean(
      roc_obj = m2$roc,
      auc_val = as.numeric(m2$auc),
      ci_val  = as.numeric(m2$ci),
      title = paste0("Model 2: Protein ", prot),
      filename_base = paste0("ROC_Model2_", prot)
    )
  }
}

# Model 3
for(prot in names(model3_list)){
  m3 <- model3_list[[prot]]
  if(!is.null(m3$roc)){
    plot_roc_stata_clean(
      roc_obj = m3$roc,
      auc_val = as.numeric(m3$auc),
      ci_val  = as.numeric(m3$ci),
      title = paste0("Model 3: ", prot, " + NIHSS0"),
      filename_base = paste0("ROC_Model3_", prot, "_NIHSS0")
    )
  }
}

# 4) Model 4: Concordant DAP pairs + NIHSS0 
for(combo in names(model4_list)){ 
  m4 <- model4_list[[combo]] 
  if(!is.null(m4$roc)){ 
    plot_roc_stata_clean( 
      roc_obj = m4$roc, 
      auc_val = as.numeric(m4$auc), 
      ci_val = as.numeric(m4$ci), 
      title = paste0("Model 4: ", combo, " + NIHSS0"), 
      filename_base = paste0("ROC_Model4_", combo, "_NIHSS0") 
    ) 
  } 
}


# ------------------- Overlay ROC Plot: NIHSS0 vs Selected Model 4 Pairs -------------------

library(pROC)
library(dplyr)

# -------- Model 1 ROC --------
roc_list <- list(model1$roc)
auc_list <- list(as.numeric(model1$auc))
ci_list  <- list(as.numeric(model1$ci))
labels   <- c("Model 1: NIHSS0")

# -------- Selected Model 4 pairs --------
selected_pairs <- c(
  "GH2.10978-39_ACP2.9237-54",
  "CD2.7100-31_ACP2.9237-54"
)

# Filter only those actually computed in model4_list
valid_pairs <- selected_pairs[selected_pairs %in% names(model4_list)]

# Append ROC objects for selected Model 4 pairs
for(pair in valid_pairs){
  m4 <- model4_list[[pair]]
  if(!is.null(m4$roc)){
    roc_list[[length(roc_list) + 1]] <- m4$roc
    auc_list[[length(auc_list) + 1]] <- as.numeric(m4$auc)
    ci_list[[length(ci_list) + 1]]   <- as.numeric(m4$ci)
    labels <- c(labels, paste0("Model 4: ", pair, " + NIHSS0"))
  }
}

# ------------------- Plot Overlay -------------------
pdf("ROC_Overlay_Model1_vs_Model4_Selected_Pairs.pdf", width = 6, height = 6)

# Styling
colors <- c("blue", "red", "darkgreen", "purple")[1:length(roc_list)]
lty    <- c(2, rep(1, length(roc_list)-1))   # dashed for NIHSS0, solid for biomarker pairs
lwd    <- c(2, rep(3, length(roc_list)-1))   # thicker for protein pairs

# Base plot: Model 1
plot(roc_list[[1]],
     col = colors[1],
     lwd = lwd[1],
     lty = lty[1],
     legacy.axes = TRUE,
     main = "Overlay ROC: NIHSS0 vs NIHSS0 + Selected Model 4 Pairs")
abline(a = 0, b = 1, lty = 2, col = "grey")

# Add curves for selected Model 4 pairs
if(length(roc_list) > 1){
  for(i in 2:length(roc_list)){
    plot(roc_list[[i]],
         col = colors[i],
         lwd = lwd[i],
         lty = lty[i],
         add = TRUE,
         legacy.axes = TRUE)
  }
}

# Legend with AUC and 95% CI
legend_text <- sapply(seq_along(labels), function(i){
  paste0(labels[i],
         ": AUC = ", round(auc_list[[i]], 3),
         " (95% CI: ", round(ci_list[[i]][1], 3),
         "-", round(ci_list[[i]][3], 3), ")")
})

legend("bottomright",
       legend = legend_text,
       col = colors,
       lwd = lwd,
       lty = lty,
       bty = "n",
       cex = 0.55)

dev.off()

##########################################################################################################################################

######################################################################
# Diagnostic Metrics for Model 1 and Model 4 using fixed cutoff (0.3)
######################################################################

library(dplyr)

# ------------------- Bootstrap sensitivity & specificity -------------------
bootstrap_sens_spec <- function(fit, data, cutoff, B = 1000) {
  
  vars <- all.vars(formula(fit))
  tmp  <- data %>% select(all_of(vars)) %>% na.omit()
  n    <- nrow(tmp)
  
  boot <- replicate(B, {
    
    idx <- sample(seq_len(n), replace = TRUE)
    d   <- tmp[idx, ]
    
    preds <- predict(fit, newdata = d, type = "response")
    predc <- ifelse(preds >= cutoff, "LVO", "noLVO")
    
    cm <- table(
      Predicted = factor(predc, levels = c("LVO","noLVO")),
      Actual    = factor(d$Group, levels = c("LVO","noLVO"))
    )
    
    TP <- cm["LVO","LVO"]
    FP <- cm["LVO","noLVO"]
    TN <- cm["noLVO","noLVO"]
    FN <- cm["noLVO","LVO"]
    
    c(
      sens = ifelse(TP + FN > 0, TP / (TP + FN), NA),
      spec = ifelse(TN + FP > 0, TN / (TN + FP), NA)
    )
  })
  
  boot
}

# ------------------- Extract metrics -------------------
extract_sens_spec <- function(fit, data, cutoff, B = 1000) {
  
  vars <- all.vars(formula(fit))
  tmp  <- data %>% select(all_of(vars)) %>% na.omit()
  
  preds <- predict(fit, newdata = tmp, type = "response")
  predc <- ifelse(preds >= cutoff, "LVO", "noLVO")
  
  cm <- table(
    Predicted = factor(predc, levels = c("LVO","noLVO")),
    Actual    = factor(tmp$Group, levels = c("LVO","noLVO"))
  )
  
  TP <- cm["LVO","LVO"]
  FP <- cm["LVO","noLVO"]
  TN <- cm["noLVO","noLVO"]
  FN <- cm["noLVO","LVO"]
  
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  
  boot <- bootstrap_sens_spec(fit, data, cutoff, B)
  
  sens_ci <- quantile(boot["sens", ], c(0.025, 0.975), na.rm = TRUE)
  spec_ci <- quantile(boot["spec", ], c(0.025, 0.975), na.rm = TRUE)
  
  tibble(
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    Sensitivity = sens,
    Sensitivity_CI_lower = sens_ci[1],
    Sensitivity_CI_upper = sens_ci[2],
    Specificity = spec,
    Specificity_CI_lower = spec_ci[1],
    Specificity_CI_upper = spec_ci[2],
    LR_positive = sens / (1 - spec),
    LR_negative = (1 - sens) / spec
  )
}

# ------------------- PPV / NPV -------------------
compute_posttest <- function(sens, spec, pretest = seq(0.05, 0.95, by = 0.05)) {
  
  tibble(
    Pretest_Probability = pretest,
    PPV = (sens * pretest) /
      (sens * pretest + (1 - spec) * (1 - pretest)),
    NPV = (spec * (1 - pretest)) /
      ((1 - sens) * pretest + spec * (1 - pretest))
  )
}

bootstrap_posttest <- function(boot_mat, pretest_grid) {
  
  ppv_mat <- sapply(pretest_grid, function(p) {
    (boot_mat["sens", ] * p) /
      (boot_mat["sens", ] * p + (1 - boot_mat["spec", ]) * (1 - p))
  })
  
  npv_mat <- sapply(pretest_grid, function(p) {
    (boot_mat["spec", ] * (1 - p)) /
      ((1 - boot_mat["sens", ]) * p + boot_mat["spec", ] * (1 - p))
  })
  
  tibble(
    Pretest_Probability = pretest_grid,
    PPV_low  = apply(ppv_mat, 2, quantile, 0.025, na.rm = TRUE),
    PPV_high = apply(ppv_mat, 2, quantile, 0.975, na.rm = TRUE),
    NPV_low  = apply(npv_mat, 2, quantile, 0.025, na.rm = TRUE),
    NPV_high = apply(npv_mat, 2, quantile, 0.975, na.rm = TRUE)
  )
}

# ============================================================
# APPLY TO MODEL 1 + MODEL 4 with fixed manual cutoff
# ============================================================

pretest_models <- list(
  Model1_NIHSS0 = model1
)

if (length(model4_list) > 0) {
  for (nm in names(model4_list)) {
    pretest_models[[paste0("Model4_", nm)]] <- model4_list[[nm]]
  }
}

metrics_list <- list()
curves_list  <- list()
bands_list   <- list()

pretest_grid <- seq(0.05, 0.95, by = 0.05)

for (nm in names(pretest_models)) {
  
  obj <- pretest_models[[nm]]
  fit <- obj$model
  
  if (is.null(fit)) next
  
  # Fixed cutoff
  cutoff <- 0.3
  
  core <- extract_sens_spec(
    fit = fit,
    data = model_data,
    cutoff = cutoff,
    B = 1000
  ) %>%
    mutate(Model = nm,
           Cutoff = cutoff)
  
  curves <- compute_posttest(
    sens = core$Sensitivity,
    spec = core$Specificity,
    pretest = pretest_grid
  ) %>%
    mutate(Model = nm)
  
  boot <- bootstrap_sens_spec(
    fit,
    model_data,
    cutoff = cutoff,
    B = 1000
  )
  
  bands <- bootstrap_posttest(boot, pretest_grid) %>%
    mutate(Model = nm)
  
  metrics_list[[nm]] <- core
  curves_list[[nm]]  <- curves
  bands_list[[nm]]   <- bands
}

metrics_df <- bind_rows(metrics_list)
curves_df  <- bind_rows(curves_list)
bands_df   <- bind_rows(bands_list)

ppv_npv_combined_df <- curves_df %>%
  left_join(bands_df, by = c("Model", "Pretest_Probability"))

# ------------------- SAVE -------------------

write.csv(
  metrics_df,
  "Diagnostic_metrics_Model1_Model4_Prob0.3_YNHH.csv",
  row.names = FALSE
)

write.csv(
  ppv_npv_combined_df,
  "PPV_NPV_Model1_Model4_Prob0.3_YNHH.csv",
  row.names = FALSE
)

##########################################################################################################################################

####################################################################################################
# Overlay PPV & NPV Curves with bootstrap bands and reference lines for Model 1 and select Model 4
####################################################################################################

library(dplyr)

# -------- Select models --------
selected_pairs <- c(
  "GH2.10978-39_ACP2.9237-54",
  "CD2.7100-31_ACP2.9237-54"
)

# Match naming used in your pipeline
model_names <- c(
  "Model1_NIHSS0",
  paste0("Model4_", selected_pairs)
)

valid_models <- model_names[model_names %in% unique(ppv_npv_combined_df$Model)]

# -------- Subset data --------
plot_df <- ppv_npv_combined_df %>%
  filter(Model %in% valid_models)

# -------- Plot settings --------
colors <- c("black", "red", "blue")[seq_along(valid_models)]
lty    <- c(2, rep(1, length(valid_models)-1))
lwd    <- c(2, rep(3, length(valid_models)-1))
vlines <- c(0.25, 0.5)

# ================================
# -------- PPV OVERLAY ----------
# ================================
pdf("PPV_Overlay_Model1_vs_Model4.pdf", width = 6, height = 6)

first <- TRUE
i <- 1

for (m in valid_models) {
  
  sub <- plot_df %>% filter(Model == m)
  
  if (first) {
    plot(sub$Pretest_Probability, sub$PPV,
         type = "l",
         col = colors[i], lwd = lwd[i], lty = lty[i],
         ylim = c(0,1),
         xlab = "Pretest probability",
         ylab = "Post-test probability (positive test)",
         main = "PPV: Model 1 vs Selected Model 4")
    first <- FALSE
  } else {
    lines(sub$Pretest_Probability, sub$PPV,
          col = colors[i], lwd = lwd[i], lty = lty[i])
  }
  
  # Confidence band
  polygon(c(sub$Pretest_Probability,
            rev(sub$Pretest_Probability)),
          c(sub$PPV_low,
            rev(sub$PPV_high)),
          col = adjustcolor(colors[i], alpha.f = 0.15),
          border = NA)
  
  # Annotate at reference lines
  for (xv in vlines) {
    yv <- approx(sub$Pretest_Probability, sub$PPV, xout = xv)$y
    text(xv, yv, labels = round(yv,2),
         pos = 4, col = colors[i], cex = 0.8)
  }
  
  i <- i + 1
}

abline(v = vlines, lty = 2)

legend("bottomright",
       legend = valid_models,
       col = colors,
       lwd = lwd,
       lty = lty,
       bty = "n",
       cex = 0.8)

dev.off()

# ================================
# -------- NPV OVERLAY ----------
# ================================
pdf("NPV_Overlay_Model1_vs_Model4.pdf", width = 6, height = 6)

first <- TRUE
i <- 1

for (m in valid_models) {
  
  sub <- plot_df %>% filter(Model == m)
  
  if (first) {
    plot(sub$Pretest_Probability, sub$NPV,
         type = "l",
         col = colors[i], lwd = lwd[i], lty = lty[i],
         ylim = c(0,1),
         xlab = "Pretest probability",
         ylab = "Post-test probability (negative test)",
         main = "NPV: Model 1 vs Selected Model 4")
    first <- FALSE
  } else {
    lines(sub$Pretest_Probability, sub$NPV,
          col = colors[i], lwd = lwd[i], lty = lty[i])
  }
  
  # Confidence band
  polygon(c(sub$Pretest_Probability,
            rev(sub$Pretest_Probability)),
          c(sub$NPV_low,
            rev(sub$NPV_high)),
          col = adjustcolor(colors[i], alpha.f = 0.15),
          border = NA)
  
  # Annotate
  for (xv in vlines) {
    yv <- approx(sub$Pretest_Probability, sub$NPV, xout = xv)$y
    text(xv, yv, labels = round(yv,2),
         pos = 4, col = colors[i], cex = 0.8)
  }
  
  i <- i + 1
}

abline(v = vlines, lty = 2)

legend("topright",
       legend = valid_models,
       col = colors,
       lwd = lwd,
       lty = lty,
       bty = "n",
       cex = 0.8)

dev.off()

###########################################################################################################################################

###############################################################
# Box plots of significant proteins after adjusting for NIHSS0
###############################################################

library(ggplot2)
library(dplyr)

# ------------------- Identify significant proteins with NIHSS0 -------------------
sig_proteins_nihss <- results %>%
  filter(Covariate == "NIHSS0", p_value < 0.05) %>%
  select(Protein, p_value)

# ------------------- Generate box plots -------------------
for(i in seq_len(nrow(sig_proteins_nihss))){
  prot <- sig_proteins_nihss$Protein[i]
  pval <- sig_proteins_nihss$p_value[i]
  
  # Ensure protein exists in data
  if(!(prot %in% names(df))) next
  
  df_plot <- df %>%
    select(Group, all_of(prot)) %>%
    na.omit()
  
  # Force order: LVO first, then noLVO
  df_plot$Group <- factor(df_plot$Group, levels = c("LVO", "noLVO"))
  
  p <- ggplot(df_plot, aes(x = Group, y = .data[[prot]], fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    labs(title = paste0("Protein: ", prot),
         subtitle = paste0("Bivariate p-value (with NIHSS0): ", signif(pval, 3)),
         x = "Group", y = "Protein level") +
    theme_minimal() +
    theme(
      plot.title = element_text(size=14, face="bold"),
      plot.subtitle = element_text(size=12),
      axis.title = element_text(size=12),
      axis.text = element_text(size=10),
      legend.position = "none"
    )
  
  # Save plots
  ggsave(paste0("Boxplot_", prot, ".pdf"), p, width=5, height=5, device="pdf")
  
  print(p)
}

###########################################################################################################################################

##################################################
# Continuous NRI Calculation (NIHSS, YNHH Cohort)
##################################################

library(dplyr)
library(purrr)
library(nricens)
library(tibble)
library(ggplot2)

# ------------------- Function to compute continuous NRI -------------------
compute_nri_ci <- function(model1, model4, data, nboot = 1000) {
  if ("glm" %in% class(model1)) model1 <- list(model = model1)
  if ("glm" %in% class(model4)) model4 <- list(model = model4)
  
  vars1 <- all.vars(model1$model$formula)
  vars4 <- all.vars(model4$model$formula)
  
  tmp <- data %>% select(all_of(unique(c(vars1, vars4)))) %>% na.omit()
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  p1 <- predict(model1$model, tmp, type = "response")
  p4 <- predict(model4$model, tmp, type = "response")
  
  nri_out <- nribin(
    event = tmp$Group_bin,
    p.std = p1,
    p.new = p4,
    updown = "diff",
    cut = 0,
    niter = nboot
  )
  
  tibble(
    Model4 = paste(vars4[vars4 != "Group"], collapse = " + "),
    NRI_event     = nri_out$nri["NRI+", "Estimate"],
    NRI_event_LCL = nri_out$nri["NRI+", "Lower"],
    NRI_event_UCL = nri_out$nri["NRI+", "Upper"],
    NRI_nonevent     = nri_out$nri["NRI-", "Estimate"],
    NRI_nonevent_LCL = nri_out$nri["NRI-", "Lower"],
    NRI_nonevent_UCL = nri_out$nri["NRI-", "Upper"]
  )
}

# ------------------- Compute NRI for all Model4 variants -------------------
nri_results <- map_df(model4_list, ~compute_nri_ci(model1, .x, model_data))
nri_results

write.csv(nri_results, "NRI_continuous_Model 1_vs_Model4_YNHH.csv", row.names = FALSE)

# ------------------- Function to get scatter plot data -------------------
get_scatter_data <- function(model1, model4, data) {
  vars1 <- all.vars(model1$model$formula)
  vars4 <- all.vars(model4$model$formula)
  
  tmp <- data %>% select(all_of(unique(c(vars1, vars4)))) %>% na.omit()
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  tibble(
    Model4 = paste(vars4[vars4 != "Group"], collapse = " + "),
    ID = 1:nrow(tmp),
    Group = ifelse(tmp$Group_bin == 1, "Case", "Control"),
    Prob_Model1 = predict(model1$model, tmp, type = "response"),
    Prob_Model4 = predict(model4$model, tmp, type = "response")
  )
}

# ------------------- Prepare scatter data -------------------
scatter_data <- map_df(model4_list, ~get_scatter_data(model1, .x, model_data))

# ------------------- Generate and save scatter plots -------------------
unique_models <- unique(scatter_data$Model4)

for (mod in unique_models) {
  plot_data <- scatter_data %>% filter(Model4 == mod)
  
  p <- ggplot(plot_data, aes(x = Prob_Model1, y = Prob_Model4, color = Group)) +
    geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.01, height = 0.01)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Case" = "red", "Control" = "black")) +
    labs(
      title = paste0("Continuous NRI Scatter:\n", mod),
      x = "Predicted probability Model 1",
      y = "Predicted probability Model 4",
      color = "Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10)),
      plot.margin = unit(c(1, 1, 1.5, 1), "cm"),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  # Save as PDF
  pdf_filename <- paste0(
    "Continuous_NRI_Scatter_",
    gsub("[^A-Za-z0-9]+", "_", mod),
    ".pdf"
  )
  
  ggsave(pdf_filename,
         plot = p,
         width = 6,
         height = 6.5,
         device = "pdf")
}

############################################################################################################################################

##############################################
# Integrated Discrimination Improvement (IDI)
# Model 4 vs Model 1 (YNHH cohort)
##############################################

library(dplyr)
library(purrr)
library(tibble)

# ------------------- Function to compute IDI -------------------
compute_idi_ci <- function(model1, model4, data, B = 1000) {
  
  # Ensure consistent structure
  if ("glm" %in% class(model1)) model1 <- list(model = model1)
  if ("glm" %in% class(model4)) model4 <- list(model = model4)
  
  # Extract variables
  vars1 <- all.vars(model1$model$formula)
  vars4 <- all.vars(model4$model$formula)
  
  # Prepare dataset
  tmp <- data %>%
    select(all_of(unique(c(vars1, vars4)))) %>%
    na.omit()
  
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  # Predicted probabilities
  p1 <- predict(model1$model, tmp, type = "response")
  p4 <- predict(model4$model, tmp, type = "response")
  
  # ------------------- IDI calculation -------------------
  mean_case_1 <- mean(p1[tmp$Group_bin == 1])
  mean_ctrl_1 <- mean(p1[tmp$Group_bin == 0])
  
  mean_case_4 <- mean(p4[tmp$Group_bin == 1])
  mean_ctrl_4 <- mean(p4[tmp$Group_bin == 0])
  
  idi <- (mean_case_4 - mean_case_1) - (mean_ctrl_4 - mean_ctrl_1)
  
  # ------------------- Bootstrap -------------------
  n <- nrow(tmp)
  
  boot_idi <- replicate(B, {
    
    idx <- sample(seq_len(n), replace = TRUE)
    d   <- tmp[idx, ]
    
    p1b <- predict(model1$model, d, type = "response")
    p4b <- predict(model4$model, d, type = "response")
    
    yb <- d$Group_bin
    
    mc1 <- mean(p1b[yb == 1])
    m0c1 <- mean(p1b[yb == 0])
    
    mc4 <- mean(p4b[yb == 1])
    m0c4 <- mean(p4b[yb == 0])
    
    (mc4 - mc1) - (m0c4 - m0c1)
  })
  
  # Confidence interval
  ci <- quantile(boot_idi, c(0.025, 0.975), na.rm = TRUE)
  
  # Two-sided bootstrap p-value
  p_val <- 2 * min(
    mean(boot_idi <= 0, na.rm = TRUE),
    mean(boot_idi >= 0, na.rm = TRUE)
  )
  
  # ------------------- Output -------------------
  tibble(
    Model4 = paste(vars4[vars4 != "Group"], collapse = " + "),
    
    Mean_case_Model1 = mean_case_1,
    Mean_control_Model1 = mean_ctrl_1,
    
    Mean_case_Model4 = mean_case_4,
    Mean_control_Model4 = mean_ctrl_4,
    
    IDI = idi,
    IDI_CI_lower = ci[1],
    IDI_CI_upper = ci[2],
    IDI_p_value = p_val
  )
}

# ------------------- Run IDI for all Model 4 -------------------
idi_results <- map_df(model4_list, ~compute_idi_ci(model1, .x, model_data))

# ------------------- View results -------------------
print(idi_results)

# ------------------- Save to CSV -------------------
write.csv(
  idi_results,
  "IDI_Model1_vs_Model4_YNHH.csv",
  row.names = FALSE
)

#############################################################################################################################################

#####################################################################################################################################
# Sub-analysis of Samples within 30 hours of stroke onset from Yale cohort (Sub-YNHH)
#####################################################################################################################################

####################################################################
# Differential abundance analysis within 30 hours (Sub-YNHH cohort)
####################################################################

#Load Packages
library(tidyverse)
library(plyr)
library(dplyr)
library(readxl)
library(writexl)
library(gplots)
library(ggrepel)
library(factoextra)

#Proteomics data file upload
excel_sheets("YNHH_LVO_data.xlsx")
YNHH_LVO_data= excel_sheets("YNHH_LVO_data.xlsx") %>% map(~read_xlsx("YNHH_LVO_data.xlsx",.))
YNHH_LVO_data

#Visualizing the Data
dat = YNHH_LVO_data[[5]]
View(dat)

#Creating a T-test function for multiple experiments
t_test <- function(dt,grp1,grp2){
  # Subset Total Stroke Case group and convert to numeric
  x <- dt[grp1] %>% unlist %>% as.numeric()
  # Subset Healthy Control group and convert to numeric
  y <- dt[grp2] %>% unlist %>% as.numeric()
  # Perform t-test using the mean of x and y
  result <- t.test(x, y)
  # Extract p-values from the results
  p_vals <- tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
} 

#Apply t-test function to data using plyr adply
#.margins = 1, slice by rows, .fun = t_test plus t_test arguments
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:23), grp2 = c(24:39)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:23]), as.numeric(dat[1,24:39]))$p.value

#Bind columns to create transformed data frame
dat_combine = bind_cols(dat, dat_pvals[,40])
View (dat_combine)

#Calculating log-fold change
dat_fc = dat_combine %>% 
  group_by(Protein_ID) %>% 
  dplyr::mutate(mean_LVO_case = mean(c(LVO1,	LVO2,	LVO3,	LVO4,	LVO5,	LVO6,	LVO7,	LVO8,	LVO9,	LVO10,	LVO11,	LVO12, LVO13,	LVO14,	
                                       LVO15,	LVO16,	LVO17,	LVO18,	LVO19,	LVO20,	LVO21,	LVO22)),
                mean_noLVO_case= mean(c(noLVO1,	noLVO2,	noLVO3,	noLVO4,	noLVO5,	noLVO6,	noLVO7,	noLVO8,	noLVO9,	noLVO10,	noLVO11,	
                                        noLVO12, noLVO13,	noLVO14,	noLVO15,	noLVO16)),
                log_fc = mean_LVO_case - mean_noLVO_case,
                log_pval = -1*log10(p_val))
View(dat_fc)

#Save final data with list of final data in csv file
write.csv(dat_fc, "Final_LVO_data_YNHH_30hours.csv")

#Volcano plot of log-fold change on x-axis and log p-value on y-axis
dat_fc %>% ggplot(aes(log_fc,log_pval)) + geom_point()

#Volcano plot at FC 1.2 and p<0.05

VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val))) + geom_point() + theme_minimal()
VP 
#Add vertical lines for Log2 FC and a horizontal line for p-value threshold
VP2= VP + geom_vline(xintercept = c(-0.26, 0.26), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Add a column of NAs
dat_fc$diffexpressed= "NO"
#Set Log2 FC and p-value cut-offs in the new column
dat_fc$diffexpressed[dat_fc$log_fc>0.26 & dat_fc$p_val<0.05] <- "UP"
dat_fc$diffexpressed[dat_fc$log_fc< -0.26 & dat_fc$p_val<0.05] <- "DOWN"
#Re-plot but this time color the points with "diffexpressed"
VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val), col= diffexpressed)) + geom_point() + theme_minimal()
VP
#Add lines as before..
VP2= VP + geom_vline(xintercept = c(-0.26, 0.26), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Change point colors
VP3= VP2 + scale_color_manual(values= c("blue", "black", "red"))
mycolors= c("blue", "red", "black")  
names(mycolors) = c("DOWN", "UP", "NO")
VP3= VP2 + scale_color_manual(values=mycolors)
#Create a new column "proteinlabel" that will contain names of differentially expressed protein IDs
dat_fc$proteinlabel= NA
dat_fc$proteinlabel[dat_fc$diffexpressed != "NO"] <- dat_fc$Protein_ID[dat_fc$diffexpressed != "NO"]
ggplot(data=dat_fc, aes(x= log_fc, y= -log10(p_val), col= diffexpressed, label=proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text()
View(dat_fc)

#Plot the Volcano plot using all layers used so far
X11()
ggplot(data= dat_fc, aes(x=log_fc, y= -log10(p_val), col= diffexpressed, label= proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline (xintercept = c(-0.26, 0.26), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")

#Proteins with significant observations
final_data<-dat_fc %>%
  #Filter for significant observations
  filter(log_pval >= 1.3 & (log_fc >= 0.26 | log_fc <= -0.26)) %>% 
  #Ungroup the data
  ungroup() %>% 
  #Select columns of interest
  select(Protein_ID, LVO1:LVO22, noLVO1:noLVO16, mean_LVO_case, mean_noLVO_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in csv file
write.csv(final_data, "Final_diffproteins_LVO_YNHH_30hours.csv")

########################################################################################################################################

#####################################################################################################
# Bivariable Logistic Regression analysis of Concordant DEPs across GMH and Sub-YNHH (within 30 hrs)
#####################################################################################################

library(readxl)
library(dplyr)
library(purrr)
library(openxlsx)

# ------------ USER INPUTS ------------
input_file  <- "YNHH_LVO_Data.xlsx"
sheet       <- 6

group_col   <- "Group"

meta_vars <- c("Age", "Sex", "DM", "HTN", "CAD", "AFIB", "Smoking", "Prior_Stroke", "CE", "NIHSS0", "GCS")

protein_start_col <- 14
# -------------------------------------

df <- read_excel(input_file, sheet = sheet)

# Ensure binary coding: controls = noLVO (0), cases = LVO (1)
df[[group_col]] <- factor(df[[group_col]], levels = c("noLVO", "LVO"))

protein_cols <- names(df)[protein_start_col:ncol(df)]
protein_cols

# ---------- MODEL FUNCTION ----------
run_model <- function(data, protein, covariate, outcome = "Group") {
  
  tmp <- data %>%
    select(all_of(c(outcome, protein, covariate))) %>%
    na.omit()
  
  if (length(unique(tmp[[outcome]])) < 2)
    return(c(OR=NA, CI_low=NA, CI_high=NA, p=NA))
  
  # Wrap variable names in backticks for safe formula
  formula_text <- paste0("`", outcome, "` ~ `", covariate, "` + `", protein, "`")
  fit <- glm(
    formula = as.formula(formula_text),
    data = tmp,
    family = binomial()
  )
  
  s <- summary(fit)$coefficients
  protein_row <- paste0("`", protein, "`")
  
  if (!(protein_row %in% rownames(s)))
    return(c(OR=NA, CI_low=NA, CI_high=NA, p=NA))
  
  est <- s[protein_row, "Estimate"]
  se  <- s[protein_row, "Std. Error"]
  p   <- s[protein_row, "Pr(>|z|)"]
  
  OR      <- exp(est)
  CI_low  <- exp(est - 1.96 * se)
  CI_high <- exp(est + 1.96 * se)
  
  c(OR=OR, CI_low=CI_low, CI_high=CI_high, p=p)
}

# -------- RUN MODELS ACROSS ALL PROTEINS × COVARIATES --------

results <- map_df(protein_cols, function(prot) {
  map_df(meta_vars, function(meta) {
    out <- run_model(df, protein = prot, covariate = meta)
    
    dplyr::tibble(
      Protein   = prot,
      Covariate = meta,
      OR        = out["OR"],
      CI_low    = out["CI_low"],
      CI_high   = out["CI_high"],
      p_value   = out["p"]
    )
  })
})

# -------- SAVE AS CSV --------
write.csv(results, "Bivariable_logistic_results_YNHH_30hrs.csv", row.names = FALSE)

##########################################################################################################################################

###################################################################
# Prediction models to classify LVO from non-LVO (Sub-YNHH cohort)
###################################################################

library(dplyr)
library(purrr)
library(pROC)
library(tidyr)
library(readxl)

# ------------ USER INPUTS ------------
input_file  <- "YNHH_LVO_Data.xlsx"
sheet       <- 6

group_col   <- "Group"
protein_start_col <- 14

# ------------ READ DATA ------------
df <- read_excel(input_file, sheet = sheet)

# ------------ Extract concordant DAPs ------------
concordant_DAPs <- colnames(df)[protein_start_col:ncol(df)]
print(concordant_DAPs)

# ------------ Subset model data ------------
model_data <- df %>%
  select(all_of(c(group_col, "NIHSS0", concordant_DAPs)))

# Ensure correct factor levels (controls first, cases second)
model_data$Group <- factor(model_data$Group, levels = c("noLVO","LVO"))

# ------------ Helper function for ROC + AUC ------------
run_model_auc <- function(formula_text, data) {
  
  vars <- all.vars(as.formula(formula_text))
  tmp <- data %>% select(all_of(vars)) %>% na.omit()
  
  if(length(unique(tmp$Group)) < 2){
    warning(paste("Skipping model", formula_text, "- only one class present"))
    return(list(model = NULL, roc = NULL, auc = NA, ci = c(NA,NA,NA)))
  }
  
  fit <- glm(as.formula(formula_text), data = tmp, family = binomial)
  preds <- predict(fit, newdata = tmp, type = "response")
  
  roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"), direction = "<")
  auc_val <- auc(roc_obj)
  ci_val  <- ci.auc(roc_obj)
  
  list(model = fit, roc = roc_obj, auc = auc_val, ci = ci_val)
}

# ------------ Helper function for PDF ROC plot only ------------
plot_roc_stata_clean <- function(roc_obj, auc_val, ci_val, title, filename_base){
  
  pdf(paste0(filename_base, ".pdf"), width = 6, height = 6)
  
  plot(roc_obj,
       col = "blue",
       lwd = 2,
       legacy.axes = TRUE,
       main = title,
       print.auc = FALSE)
  
  abline(a = 0, b = 1, lty = 2, col = "grey")
  
  legend("bottomright", 
         legend = paste0(
           "AUC = ", round(auc_val, 3),
           " (95% CI: ", round(ci_val[1], 3), "-", round(ci_val[3], 3), ")"
         ),
         bty = "n")
  
  dev.off()
}

# ------------ MODEL 1: NIHSS0 only ------------
model1 <- run_model_auc("Group ~ NIHSS0", model_data)

# ------------ MODEL 2: individual proteins ------------
model2_list <- map(concordant_DAPs, function(prot) {
  run_model_auc(paste0("Group ~ `", prot, "`"), model_data)
})
names(model2_list) <- concordant_DAPs

# ------------ MODEL 3: protein + NIHSS0 ------------
model3_list <- map(concordant_DAPs, function(prot) {
  run_model_auc(paste0("Group ~ NIHSS0 + `", prot, "`"), model_data)
})
names(model3_list) <- concordant_DAPs

# ------------------- Extract significant concordant DAPs from Model 3 ------------------- 
sig_proteins_m3 <- keep(names(model3_list), function(prot) { 
  fit <- model3_list[[prot]]$model
  if(is.null(fit)) return(FALSE)  # skip if model didn't run
  
  coefs <- summary(fit)$coefficients
  # Match protein term using grep to handle backticks
  prot_row <- grep(prot, rownames(coefs), fixed = TRUE)
  if(length(prot_row) == 0) return(FALSE)  # protein not found
  
  p_val <- coefs[prot_row, "Pr(>|z|)"]
  !is.na(p_val) && p_val < 0.05
})

sig_proteins_m3

# ------------------- 4) Concordant DAP pairs + NIHSS0 (filtered from Model 3) ------------------- 
model4_list <- list()

if(length(sig_proteins_m3) >= 2){ 
  combos <- combn(sig_proteins_m3, 2, simplify = FALSE) 
  
  for(combo in combos){ 
    # wrap protein names in backticks
    combo_backtick <- paste0("`", combo, "`")
    formula4 <- paste("Group ~ NIHSS0 +", paste(combo_backtick, collapse = " + "))
    
    vars <- c("Group", "NIHSS0", combo) 
    tmp <- model_data %>% select(all_of(vars)) %>% na.omit() 
    
    if(length(unique(tmp$Group)) < 2) next 
    
    fit <- glm(as.formula(formula4), data = tmp, family = binomial) 
    coefs <- summary(fit)$coefficients 
    
    # Extract p-values using grep to match backtick names
    protein_rows <- sapply(combo, function(p) grep(p, rownames(coefs), fixed = TRUE))
    protein_p <- coefs[protein_rows, "Pr(>|z|)"]
    
    if(all(!is.na(protein_p)) && all(protein_p < 0.05)){ 
      preds <- predict(fit, newdata = tmp, type = "response") 
      roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO")) 
      
      model4_list[[paste(combo, collapse = "_")]] <- list( 
        model = fit, 
        roc = roc_obj, 
        auc = auc(roc_obj), 
        ci = ci.auc(roc_obj) 
      ) 
    } 
  } 
}

model4_list

# ------------ SAVE AUC RESULTS ------------

# Model 1
auc_results <- tibble(
  Model = "NIHSS0",
  AUC = as.numeric(model1$auc),
  CI_lower = as.numeric(model1$ci[1]),
  CI_upper = as.numeric(model1$ci[3])
)

# Model 2
m2_df <- tibble(
  Model = paste0("Protein ", names(model2_list)),
  AUC = sapply(model2_list, function(x) as.numeric(x$auc)),
  CI_lower = sapply(model2_list, function(x) as.numeric(x$ci[1])),
  CI_upper = sapply(model2_list, function(x) as.numeric(x$ci[3]))
) %>% mutate(across(everything(), as.vector))

auc_results <- bind_rows(auc_results, m2_df)

# Model 3
m3_df <- tibble(
  Model = paste0(names(model3_list), " + NIHSS0"),
  AUC = sapply(model3_list, function(x) as.numeric(x$auc)),
  CI_lower = sapply(model3_list, function(x) as.numeric(x$ci[1])),
  CI_upper = sapply(model3_list, function(x) as.numeric(x$ci[3]))
) %>% mutate(across(everything(), as.vector))

auc_results <- bind_rows(auc_results, m3_df)

# Model 4: two proteins + NIHSS0 
if(length(model4_list) > 0){ 
  m4_df <- tibble( 
    Model = paste0(names(model4_list), " + NIHSS0"), 
    AUC = sapply(model4_list, function(x) if(!is.null(x$auc)) 
      as.numeric(x$auc) else NA_real_), 
    CI_lower = sapply(model4_list, function(x) if(!is.null(x$ci)) 
      as.numeric(x$ci[1]) else NA_real_), 
    CI_upper = sapply(model4_list, function(x) if(!is.null(x$ci)) 
      as.numeric(x$ci[3]) else NA_real_) 
  ) 
  m4_df <- m4_df %>% mutate(across(everything(), ~ as.vector(.))) 
  auc_results <- bind_rows(auc_results, m4_df) 
}

# ------------ DeLong test: Model 1 vs all other models ------------

get_delong_p <- function(roc_ref, roc_comp){
  if(is.null(roc_ref) || is.null(roc_comp)) return(NA_real_)
  
  # Ensure both ROC objects have valid responses
  if(length(roc_ref$response) == 0 || length(roc_comp$response) == 0) return(NA_real_)
  
  # Try DeLong test safely
  out <- tryCatch({
    test <- roc.test(roc_ref, roc_comp, method = "delong")
    as.numeric(test$p.value)
  }, error = function(e) NA_real_)
  
  return(out)
}

# Initialize vector
delong_pvals <- c()

# Model 1 (reference) → p-value = NA
delong_pvals <- c(delong_pvals, NA_real_)

# Model 2
delong_pvals <- c(
  delong_pvals,
  sapply(model2_list, function(x) get_delong_p(model1$roc, x$roc))
)

# Model 3
delong_pvals <- c(
  delong_pvals,
  sapply(model3_list, function(x) get_delong_p(model1$roc, x$roc))
)

# Model 4
if(length(model4_list) > 0){
  delong_pvals <- c(
    delong_pvals,
    sapply(model4_list, function(x) get_delong_p(model1$roc, x$roc))
  )
}

# Add to results
auc_results$DeLong_p_value <- delong_pvals

# ============================================================
# Extract OR, 95% CI, and p-values from all models
# Combine with AUC + DeLong results
# ============================================================

# ============================================================
# Helper function: Extract OR results
# ============================================================
extract_or_results <- function(fit, model_name){
  
  if(is.null(fit)){
    return(NULL)
  }
  
  coefs <- summary(fit)$coefficients
  
  # Remove intercept
  coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
  
  # Calculate OR and 95% CI
  OR  <- exp(coefs[, "Estimate"])
  LCL <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
  UCL <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
  
  # Extract p-values
  pvals <- coefs[, "Pr(>|z|)"]
  
  tibble(
    Model = model_name,
    Variable = rownames(coefs),
    OR = OR,
    CI_lower_OR = LCL,
    CI_upper_OR = UCL,
    P_value = pvals
  )
}

# ============================================================
# Model 1 OR results
# ============================================================
or_results_m1 <- extract_or_results(
  model1$model,
  "NIHSS0"
)

# ============================================================
# Model 2 OR results
# ============================================================
or_results_m2 <- map2(
  model2_list,
  names(model2_list),
  ~ extract_or_results(
    .x$model,
    paste0("Protein ", .y)
  )
)

or_results_m2 <- bind_rows(or_results_m2)

# ============================================================
# Model 3 OR results
# ============================================================
or_results_m3 <- map2(
  model3_list,
  names(model3_list),
  ~ extract_or_results(
    .x$model,
    paste0(.y, " + NIHSS0")
  )
)

or_results_m3 <- bind_rows(or_results_m3)

# ============================================================
# Model 4 OR results
# ============================================================
if(length(model4_list) > 0){
  
  or_results_m4 <- map2(
    model4_list,
    names(model4_list),
    ~ extract_or_results(
      .x$model,
      paste0(.y, " + NIHSS0")
    )
  )
  
  or_results_m4 <- bind_rows(or_results_m4)
  
} else {
  
  or_results_m4 <- tibble()
}

# ============================================================
# Combine all OR results
# ============================================================
all_or_results <- bind_rows(
  or_results_m1,
  or_results_m2,
  or_results_m3,
  or_results_m4
)

# ============================================================
# Merge OR results with AUC results
# ============================================================
final_results <- all_or_results %>%
  left_join(
    auc_results,
    by = "Model"
  )

# ============================================================
# View final results
# ============================================================
View(final_results)

# ============================================================
# Save final combined results
# ============================================================
write.csv(
  final_results,
  "LVO_prediction_models_OR_AUC_DeLong_Sub-YNHH(30hrs).csv",
  row.names = FALSE
)

# ------------ PLOT ROC CURVES ------------

# Model 1
if(!is.null(model1$roc)){
  plot_roc_stata_clean(
    roc_obj = model1$roc,
    auc_val = as.numeric(model1$auc),
    ci_val  = as.numeric(model1$ci),
    title = "Model 1: NIHSS0 only",
    filename_base = "ROC_Model1_NIHSS0"
  )
}

# Model 2
for(prot in names(model2_list)){
  m2 <- model2_list[[prot]]
  if(!is.null(m2$roc)){
    plot_roc_stata_clean(
      roc_obj = m2$roc,
      auc_val = as.numeric(m2$auc),
      ci_val  = as.numeric(m2$ci),
      title = paste0("Model 2: Protein ", prot),
      filename_base = paste0("ROC_Model2_", prot)
    )
  }
}

# Model 3
for(prot in names(model3_list)){
  m3 <- model3_list[[prot]]
  if(!is.null(m3$roc)){
    plot_roc_stata_clean(
      roc_obj = m3$roc,
      auc_val = as.numeric(m3$auc),
      ci_val  = as.numeric(m3$ci),
      title = paste0("Model 3: ", prot, " + NIHSS0"),
      filename_base = paste0("ROC_Model3_", prot, "_NIHSS0")
    )
  }
}

# 4) Model 4: Concordant DAP pairs + NIHSS0 
for(combo in names(model4_list)){ 
  m4 <- model4_list[[combo]] 
  if(!is.null(m4$roc)){ 
    plot_roc_stata_clean( 
      roc_obj = m4$roc, 
      auc_val = as.numeric(m4$auc), 
      ci_val = as.numeric(m4$ci), 
      title = paste0("Model 4: ", combo, " + NIHSS0"), 
      filename_base = paste0("ROC_Model4_", combo, "_NIHSS0") 
    ) 
  } 
}


# ------------------- Overlay ROC Plot: NIHSS0 vs Selected Model 4 Pairs -------------------

library(pROC)
library(dplyr)

# -------- Model 1 ROC --------
roc_list <- list(model1$roc)
auc_list <- list(as.numeric(model1$auc))
ci_list  <- list(as.numeric(model1$ci))
labels   <- c("Model 1: NIHSS0")

# -------- Selected Model 4 pair --------
selected_pairs <- c(
  "GH2.10978-39_ACP2.9237-54"
)

# Filter only those actually computed in model4_list
valid_pairs <- selected_pairs[selected_pairs %in% names(model4_list)]

# Append ROC objects for selected Model 4 pairs
for(pair in valid_pairs){
  m4 <- model4_list[[pair]]
  if(!is.null(m4$roc)){
    roc_list[[length(roc_list) + 1]] <- m4$roc
    auc_list[[length(auc_list) + 1]] <- as.numeric(m4$auc)
    ci_list[[length(ci_list) + 1]]   <- as.numeric(m4$ci)
    labels <- c(labels, paste0("Model 4: ", pair, " + NIHSS0"))
  }
}

# ------------------- Plot Overlay -------------------
pdf("ROC_Overlay_Model1_vs_Model4_Selected_Pair.pdf", width = 6, height = 6)

# Styling
colors <- c("blue", "red", "darkgreen", "purple")[1:length(roc_list)]
lty    <- c(2, rep(1, length(roc_list)-1))   # dashed for NIHSS0, solid for biomarker pairs
lwd    <- c(2, rep(3, length(roc_list)-1))   # thicker for protein pairs

# Base plot: Model 1
plot(roc_list[[1]],
     col = colors[1],
     lwd = lwd[1],
     lty = lty[1],
     legacy.axes = TRUE,
     main = "Overlay ROC: NIHSS0 vs NIHSS0 + Selected Model 4 Pairs")
abline(a = 0, b = 1, lty = 2, col = "grey")

# Add curves for selected Model 4 pairs
if(length(roc_list) > 1){
  for(i in 2:length(roc_list)){
    plot(roc_list[[i]],
         col = colors[i],
         lwd = lwd[i],
         lty = lty[i],
         add = TRUE,
         legacy.axes = TRUE)
  }
}

# Legend with AUC and 95% CI
legend_text <- sapply(seq_along(labels), function(i){
  paste0(labels[i],
         ": AUC = ", round(auc_list[[i]], 3),
         " (95% CI: ", round(ci_list[[i]][1], 3),
         "-", round(ci_list[[i]][3], 3), ")")
})

legend("bottomright",
       legend = legend_text,
       col = colors,
       lwd = lwd,
       lty = lty,
       bty = "n",
       cex = 0.55)

dev.off()

###########################################################################################################################################

##########################################################################
# Box plots of significant proteins after adjusting for NIHSS0 (Sub-YNHH)
##########################################################################

library(ggplot2)
library(dplyr)

# ------------------- Identify significant proteins with NIHSS0 -------------------
sig_proteins_nihss <- results %>%
  filter(Covariate == "NIHSS0", p_value < 0.05) %>%
  select(Protein, p_value)

# ------------------- Generate box plots -------------------
for(i in seq_len(nrow(sig_proteins_nihss))){
  prot <- sig_proteins_nihss$Protein[i]
  pval <- sig_proteins_nihss$p_value[i]
  
  # Ensure protein exists in data
  if(!(prot %in% names(df))) next
  
  df_plot <- df %>%
    select(Group, all_of(prot)) %>%
    na.omit()
  
  # Force order: LVO first, then noLVO
  df_plot$Group <- factor(df_plot$Group, levels = c("LVO", "noLVO"))
  
  p <- ggplot(df_plot, aes(x = Group, y = .data[[prot]], fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    labs(title = paste0("Protein: ", prot),
         subtitle = paste0("Bivariate p-value (with NIHSS0): ", signif(pval, 3)),
         x = "Group", y = "Protein level") +
    theme_minimal() +
    theme(
      plot.title = element_text(size=14, face="bold"),
      plot.subtitle = element_text(size=12),
      axis.title = element_text(size=12),
      axis.text = element_text(size=10),
      legend.position = "none"
    )
  
  # Save plots
  ggsave(paste0("Boxplot_", prot, ".pdf"), p, width=5, height=5, device="pdf")
  
  print(p)
}

##############################################################################################################################################

#####################################################################################################################################
# Sub-analysis of Yale Cohort with Prehospital scale scores (RACE, G-FAST, mG-FAST)
#####################################################################################################################################

##########################################################################################################
# Bivariable Logistic Regression analysis of Concordant DEPs with Prehospital stroke scales (YNHH Cohort)
##########################################################################################################

library(readxl)
library(dplyr)
library(purrr)
library(openxlsx)

# ------------ USER INPUTS ------------
input_file  <- "YNHH_LVO_data.xlsx"
sheet       <- 7

group_col   <- "Group"

meta_vars <- c("RACE", "RACE_cutoff", "G_FAST", "mG_FAST") #Three prehospital stroke scales

protein_start_col <- 18
# -------------------------------------

df <- read_excel(input_file, sheet = sheet)

# Ensure binary coding: controls = noLVO (0), cases = LVO (1)
df[[group_col]] <- factor(df[[group_col]], levels = c("noLVO", "LVO"))

protein_cols <- names(df)[protein_start_col:ncol(df)]
protein_cols

# ---------- MODEL FUNCTION ----------
run_model <- function(data, protein, covariate, outcome = "Group") {
  
  tmp <- data %>%
    select(all_of(c(outcome, protein, covariate))) %>%
    na.omit()
  
  if (length(unique(tmp[[outcome]])) < 2)
    return(c(OR=NA, CI_low=NA, CI_high=NA, p=NA))
  
  # Wrap variable names in backticks for safe formula
  formula_text <- paste0("`", outcome, "` ~ `", covariate, "` + `", protein, "`")
  fit <- glm(
    formula = as.formula(formula_text),
    data = tmp,
    family = binomial()
  )
  
  s <- summary(fit)$coefficients
  protein_row <- paste0("`", protein, "`")
  
  if (!(protein_row %in% rownames(s)))
    return(c(OR=NA, CI_low=NA, CI_high=NA, p=NA))
  
  est <- s[protein_row, "Estimate"]
  se  <- s[protein_row, "Std. Error"]
  p   <- s[protein_row, "Pr(>|z|)"]
  
  OR      <- exp(est)
  CI_low  <- exp(est - 1.96 * se)
  CI_high <- exp(est + 1.96 * se)
  
  c(OR=OR, CI_low=CI_low, CI_high=CI_high, p=p)
}

# -------- RUN MODELS ACROSS ALL PROTEINS × COVARIATES --------

results <- map_df(protein_cols, function(prot) {
  map_df(meta_vars, function(meta) {
    out <- run_model(df, protein = prot, covariate = meta)
    
    dplyr::tibble(
      Protein   = prot,
      Covariate = meta,
      OR        = out["OR"],
      CI_low    = out["CI_low"],
      CI_high   = out["CI_high"],
      p_value   = out["p"]
    )
  })
})

# -------- SAVE AS CSV --------
write.csv(results, "Bivariable_logistic_results_Scales_YNHH.csv", row.names = FALSE)

##############################################################################################################################################

##############################################################################################################################################
#RACE prehospital stroke scale (YNHH)
##############################################################################################################################################

# RACE: Prediction models for LVO (AUC + DeLong comparisons)

library(dplyr)
library(purrr)
library(pROC)
library(tidyr)

# ------------------- Significant proteins -------------------
sig_proteins <- results %>%
  filter(Covariate == "RACE", p_value < 0.05) %>%
  pull(Protein)

# ------------------- Model data -------------------
model_data <- df %>%
  select(all_of(c("Group", "RACE", sig_proteins)))

model_data$Group <- factor(model_data$Group, levels = c("noLVO","LVO"))

# ============================================================
# Helper: ROC + AUC
# ============================================================
run_model_auc <- function(formula_text, data) {
  
  vars <- all.vars(as.formula(formula_text))
  tmp <- data %>% select(all_of(vars)) %>% na.omit()
  
  if(length(unique(tmp$Group)) < 2){
    return(list(model = NULL, roc = NULL, auc = NA, ci = c(NA,NA,NA)))
  }
  
  fit <- glm(as.formula(formula_text), data = tmp, family = binomial)
  preds <- predict(fit, type = "response")
  
  roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"))
  
  list(
    model = fit,
    roc   = roc_obj,
    auc   = auc(roc_obj),
    ci    = ci.auc(roc_obj)
  )
}

# ============================================================
# Model 1 (reference)
# ============================================================
model1 <- run_model_auc("Group ~ RACE", model_data)

# ============================================================
# Model 2 (single protein + RACE)
# ============================================================
model2_list <- map(sig_proteins, function(prot) {
  run_model_auc(paste0("Group ~ RACE + `", prot, "`"), model_data)
})
names(model2_list) <- sig_proteins

# ============================================================
# Model 3 (two proteins + RACE)
# ============================================================
model3_list <- list()

if(length(sig_proteins) >= 2){
  
  combos <- combn(sig_proteins, 2, simplify = FALSE)
  
  for(combo in combos){
    
    formula3 <- paste("Group ~ RACE +", paste0("`", combo, "`", collapse = " + "))
    vars <- c("Group", "RACE", combo)
    
    tmp <- model_data %>% select(all_of(vars)) %>% na.omit()
    
    if(length(unique(tmp$Group)) < 2) next
    
    fit <- glm(as.formula(formula3), data = tmp, family = binomial)
    
    coefs <- summary(fit)$coefficients
    
    prot_rows <- sapply(combo, function(p) grep(p, rownames(coefs), fixed = TRUE))
    protein_p <- coefs[prot_rows, "Pr(>|z|)"]
    
    if(all(!is.na(protein_p)) && all(protein_p < 0.05)){
      
      preds <- predict(fit, type = "response")
      roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"))
      
      model3_list[[paste(combo, collapse = "_")]] <- list(
        model = fit,
        roc   = roc_obj,
        auc   = auc(roc_obj),
        ci    = ci.auc(roc_obj)
      )
    }
  }
}

# ============================================================
# Base AUC table (Model 1)
# ============================================================
auc_results <- tibble(
  Model = "RACE",
  AUC = as.numeric(model1$auc),
  CI_lower = as.numeric(model1$ci[1]),
  CI_upper = as.numeric(model1$ci[3])
)

# ============================================================
# Model 2 table + DeLong vs Model 1
# ============================================================
if(length(model2_list) > 0){
  
  model2_df <- tibble(
    Model = paste0(names(model2_list), " + RACE"),
    AUC = sapply(model2_list, function(x) as.numeric(x$auc)),
    CI_lower = sapply(model2_list, function(x) as.numeric(x$ci[1])),
    CI_upper = sapply(model2_list, function(x) as.numeric(x$ci[3])),
    
    Delong_p_value = sapply(model2_list, function(x) {
      if(is.null(model1$roc) || is.null(x$roc)) return(NA_real_)
      roc.test(model1$roc, x$roc, method = "delong")$p.value
    })
  )
  
  auc_results <- bind_rows(auc_results, model2_df)
}

# ============================================================
# Model 3 table + DeLong vs Model 1
# ============================================================
if(length(model3_list) > 0){
  
  model3_df <- tibble(
    Model = paste0(names(model3_list), " + RACE"),
    AUC = sapply(model3_list, function(x) as.numeric(x$auc)),
    CI_lower = sapply(model3_list, function(x) as.numeric(x$ci[1])),
    CI_upper = sapply(model3_list, function(x) as.numeric(x$ci[3])),
    
    Delong_p_value = sapply(model3_list, function(x) {
      if(is.null(model1$roc) || is.null(x$roc)) return(NA_real_)
      roc.test(model1$roc, x$roc, method = "delong")$p.value
    })
  )
  
  auc_results <- bind_rows(auc_results, model3_df)
}

# ============================================================
# Extract OR, 95% CI, and p-values from all models
# and combine with AUC + DeLong results
# ============================================================

# ============================================================
# Helper function: Extract OR results
# ============================================================
extract_or_results <- function(fit, model_name){
  
  if(is.null(fit)){
    return(NULL)
  }
  
  coefs <- summary(fit)$coefficients
  
  # Remove intercept
  coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
  
  # Calculate OR and 95% CI
  OR  <- exp(coefs[, "Estimate"])
  LCL <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
  UCL <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
  
  # Extract p-values
  pvals <- coefs[, "Pr(>|z|)"]
  
  tibble(
    Model = model_name,
    Variable = rownames(coefs),
    OR = OR,
    CI_lower_OR = LCL,
    CI_upper_OR = UCL,
    P_value = pvals
  )
}

# ============================================================
# Extract OR results: Model 1
# ============================================================
or_results_m1 <- extract_or_results(
  model1$model,
  "RACE"
)

# ============================================================
# Extract OR results: Model 2
# ============================================================
or_results_m2 <- map2(
  model2_list,
  names(model2_list),
  ~ extract_or_results(
    .x$model,
    paste0(.y, " + RACE")
  )
)

or_results_m2 <- bind_rows(or_results_m2)

# ============================================================
# Extract OR results: Model 3
# ============================================================
if(length(model3_list) > 0){
  
  or_results_m3 <- map2(
    model3_list,
    names(model3_list),
    ~ extract_or_results(
      .x$model,
      paste0(.y, " + RACE")
    )
  )
  
  or_results_m3 <- bind_rows(or_results_m3)
  
} else {
  
  or_results_m3 <- tibble()
}

# ============================================================
# Combine all OR results
# ============================================================
all_or_results <- bind_rows(
  or_results_m1,
  or_results_m2,
  or_results_m3
)

# ============================================================
# Merge OR results with AUC results
# ============================================================
final_results <- all_or_results %>%
  left_join(
    auc_results,
    by = "Model"
  )

# ============================================================
# View final combined results
# ============================================================
View(final_results)

# ============================================================
# Save final results
# ============================================================
write.csv(
  final_results,
  "LVO_prediction_models_OR_AUC_Delong_YNHH_RACE.csv",
  row.names = FALSE
)

# ============================================================
# ROC overlay: Model 1 vs selected Model 3 pairs
# ============================================================

library(dplyr)
library(pROC)

# ------------------- Selected Model 3 pairs -------------------
selected_pairs <- c(
  "GH2.10978-39_ACP2.9237-54",
  "CD2.7100-31_ACP2.9237-54"
)

valid_pairs <- selected_pairs[selected_pairs %in% names(model3_list)]

# ------------------- Colors -------------------
cols <- c("blue", "red", "darkgreen")

# ------------------- Open single PDF -------------------
pdf("ROC_Overlay_Model1_vs_Selected pairs_Model3_RACE.pdf",
    width = 6.5, height = 6.5)

# ------------------- Base plot: Model 1 -------------------
plot(model1$roc,
     col = cols[1],
     lwd = 2,
     legacy.axes = TRUE,
     main = "ROC Overlay: Model 1 vs Selected Model 3 Pairs (RACE)")

# ------------------- Add Model 3 curves -------------------
i <- 2

for (pair in valid_pairs) {
  
  m3 <- model3_list[[pair]]
  
  if (is.null(m3$roc)) next
  
  plot(m3$roc,
       col = cols[i],
       lwd = 2,
       add = TRUE)
  
  i <- i + 1
}

# ------------------- Reference line -------------------
abline(a = 0, b = 1, lty = 2, col = "grey")

# ------------------- Legend -------------------
legend_labels <- c(
  paste0("Model 1: RACE (AUC = ", round(model1$auc, 3), ")"),
  paste0(valid_pairs[1], " + RACE (AUC = ", round(model3_list[[valid_pairs[1]]]$auc, 3), ")"),
  paste0(valid_pairs[2], " + RACE (AUC = ", round(model3_list[[valid_pairs[2]]]$auc, 3), ")")
)

legend("bottomright",
       legend = legend_labels,
       col = cols,
       lwd = 2,
       bty = "n")

# ------------------- Close PDF -------------------
dev.off()

##########################################################################################################################################

#############################################################################
# Diagnostic Metrics for Model 1 and Model 3 using fixed cutoff (YNHH, RACE)
#############################################################################

library(dplyr)

# ------------------- Bootstrap sensitivity & specificity -------------------
bootstrap_sens_spec <- function(fit, data, cutoff, B = 1000) {
  
  vars <- all.vars(formula(fit))
  tmp  <- data %>% select(all_of(vars)) %>% na.omit()
  n    <- nrow(tmp)
  
  boot <- replicate(B, {
    
    idx <- sample(seq_len(n), replace = TRUE)
    d   <- tmp[idx, ]
    
    preds <- predict(fit, newdata = d, type = "response")
    predc <- ifelse(preds >= cutoff, "LVO", "noLVO")
    
    cm <- table(
      Predicted = factor(predc, levels = c("LVO","noLVO")),
      Actual    = factor(d$Group, levels = c("LVO","noLVO"))
    )
    
    TP <- cm["LVO","LVO"]
    FP <- cm["LVO","noLVO"]
    TN <- cm["noLVO","noLVO"]
    FN <- cm["noLVO","LVO"]
    
    c(
      sens = ifelse(TP + FN > 0, TP / (TP + FN), NA),
      spec = ifelse(TN + FP > 0, TN / (TN + FP), NA)
    )
  })
  
  boot
}

# ------------------- Extract metrics -------------------
extract_sens_spec <- function(fit, data, cutoff, B = 1000) {
  
  vars <- all.vars(formula(fit))
  tmp  <- data %>% select(all_of(vars)) %>% na.omit()
  
  preds <- predict(fit, newdata = tmp, type = "response")
  predc <- ifelse(preds >= cutoff, "LVO", "noLVO")
  
  cm <- table(
    Predicted = factor(predc, levels = c("LVO","noLVO")),
    Actual    = factor(tmp$Group, levels = c("LVO","noLVO"))
  )
  
  TP <- cm["LVO","LVO"]
  FP <- cm["LVO","noLVO"]
  TN <- cm["noLVO","noLVO"]
  FN <- cm["noLVO","LVO"]
  
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  
  boot <- bootstrap_sens_spec(fit, data, cutoff, B)
  
  sens_ci <- quantile(boot["sens", ], c(0.025, 0.975), na.rm = TRUE)
  spec_ci <- quantile(boot["spec", ], c(0.025, 0.975), na.rm = TRUE)
  
  tibble(
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    Sensitivity = sens,
    Sensitivity_CI_lower = sens_ci[1],
    Sensitivity_CI_upper = sens_ci[2],
    Specificity = spec,
    Specificity_CI_lower = spec_ci[1],
    Specificity_CI_upper = spec_ci[2],
    LR_positive = sens / (1 - spec),
    LR_negative = (1 - sens) / spec
  )
}

# ------------------- PPV / NPV -------------------
compute_posttest <- function(sens, spec, pretest = seq(0.05, 0.95, by = 0.05)) {
  
  tibble(
    Pretest_Probability = pretest,
    PPV = (sens * pretest) /
      (sens * pretest + (1 - spec) * (1 - pretest)),
    NPV = (spec * (1 - pretest)) /
      ((1 - sens) * pretest + spec * (1 - pretest))
  )
}

bootstrap_posttest <- function(boot_mat, pretest_grid) {
  
  ppv_mat <- sapply(pretest_grid, function(p) {
    (boot_mat["sens", ] * p) /
      (boot_mat["sens", ] * p + (1 - boot_mat["spec", ]) * (1 - p))
  })
  
  npv_mat <- sapply(pretest_grid, function(p) {
    (boot_mat["spec", ] * (1 - p)) /
      ((1 - boot_mat["sens", ]) * p + boot_mat["spec", ] * (1 - p))
  })
  
  tibble(
    Pretest_Probability = pretest_grid,
    PPV_low  = apply(ppv_mat, 2, quantile, 0.025, na.rm = TRUE),
    PPV_high = apply(ppv_mat, 2, quantile, 0.975, na.rm = TRUE),
    NPV_low  = apply(npv_mat, 2, quantile, 0.025, na.rm = TRUE),
    NPV_high = apply(npv_mat, 2, quantile, 0.975, na.rm = TRUE)
  )
}

# ============================================================
# APPLY TO MODEL 1 + MODEL 3 with fixed manual cutoff
# ============================================================

pretest_models <- list(
  Model1_RACE = model1
)

if (length(model3_list) > 0) {
  for (nm in names(model3_list)) {
    pretest_models[[paste0("Model3_", nm)]] <- model3_list[[nm]]
  }
}

metrics_list <- list()
curves_list  <- list()
bands_list   <- list()

pretest_grid <- seq(0.05, 0.95, by = 0.05)

for (nm in names(pretest_models)) {
  
  obj <- pretest_models[[nm]]
  fit <- obj$model
  
  if (is.null(fit)) next
  
  # Fixed cutoff
  cutoff <- 0.3
  
  core <- extract_sens_spec(
    fit = fit,
    data = model_data,
    cutoff = cutoff,
    B = 1000
  ) %>%
    mutate(Model = nm,
           Cutoff = cutoff)
  
  curves <- compute_posttest(
    sens = core$Sensitivity,
    spec = core$Specificity,
    pretest = pretest_grid
  ) %>%
    mutate(Model = nm)
  
  boot <- bootstrap_sens_spec(
    fit,
    model_data,
    cutoff = cutoff,
    B = 1000
  )
  
  bands <- bootstrap_posttest(boot, pretest_grid) %>%
    mutate(Model = nm)
  
  metrics_list[[nm]] <- core
  curves_list[[nm]]  <- curves
  bands_list[[nm]]   <- bands
}

metrics_df <- bind_rows(metrics_list)
curves_df  <- bind_rows(curves_list)
bands_df   <- bind_rows(bands_list)

ppv_npv_combined_df <- curves_df %>%
  left_join(bands_df, by = c("Model", "Pretest_Probability"))

# ------------------- SAVE -------------------

write.csv(
  metrics_df,
  "Diagnostic_metrics_Model1_Model3_Prob0.3_RACE_YNHH.csv",
  row.names = FALSE
)

write.csv(
  ppv_npv_combined_df,
  "PPV_NPV_Model1_Model3_Prob0.3_RACE_YNHH.csv",
  row.names = FALSE
)

##########################################################################################################################################

################################################################################################################
# Overlay PPV & NPV Curves with bootstrap bands and reference lines for Model 1 and select Model 3 (YNHH, RACE)
################################################################################################################

library(dplyr)

# -------- Select models --------
selected_pairs <- c(
  "GH2.10978-39_ACP2.9237-54",
  "CD2.7100-31_ACP2.9237-54"
)

# Match naming used in your pipeline
model_names <- c(
  "Model1_RACE",
  paste0("Model3_", selected_pairs)
)

valid_models <- model_names[model_names %in% unique(ppv_npv_combined_df$Model)]

# -------- Subset data --------
plot_df <- ppv_npv_combined_df %>%
  filter(Model %in% valid_models)

# -------- Plot settings --------
colors <- c("black", "red", "blue")[seq_along(valid_models)]
lty    <- c(2, rep(1, length(valid_models)-1))
lwd    <- c(2, rep(3, length(valid_models)-1))
vlines <- c(0.25, 0.5)

# ================================
# -------- PPV OVERLAY ----------
# ================================
pdf("PPV_Overlay_Model1_vs_Model3_RACE.pdf", width = 6, height = 6)

first <- TRUE
i <- 1

for (m in valid_models) {
  
  sub <- plot_df %>% filter(Model == m)
  
  if (first) {
    plot(sub$Pretest_Probability, sub$PPV,
         type = "l",
         col = colors[i], lwd = lwd[i], lty = lty[i],
         ylim = c(0,1),
         xlab = "Pretest probability",
         ylab = "Post-test probability (positive test)",
         main = "PPV: Model 1 vs Selected Model 3")
    first <- FALSE
  } else {
    lines(sub$Pretest_Probability, sub$PPV,
          col = colors[i], lwd = lwd[i], lty = lty[i])
  }
  
  # Confidence band
  polygon(c(sub$Pretest_Probability,
            rev(sub$Pretest_Probability)),
          c(sub$PPV_low,
            rev(sub$PPV_high)),
          col = adjustcolor(colors[i], alpha.f = 0.15),
          border = NA)
  
  # Annotate at reference lines
  for (xv in vlines) {
    yv <- approx(sub$Pretest_Probability, sub$PPV, xout = xv)$y
    text(xv, yv, labels = round(yv,2),
         pos = 4, col = colors[i], cex = 0.8)
  }
  
  i <- i + 1
}

abline(v = vlines, lty = 2)

legend("bottomright",
       legend = valid_models,
       col = colors,
       lwd = lwd,
       lty = lty,
       bty = "n",
       cex = 0.8)

dev.off()

# ================================
# -------- NPV OVERLAY ----------
# ================================
pdf("NPV_Overlay_Model1_vs_Model3_RACE.pdf", width = 6, height = 6)

first <- TRUE
i <- 1

for (m in valid_models) {
  
  sub <- plot_df %>% filter(Model == m)
  
  if (first) {
    plot(sub$Pretest_Probability, sub$NPV,
         type = "l",
         col = colors[i], lwd = lwd[i], lty = lty[i],
         ylim = c(0,1),
         xlab = "Pretest probability",
         ylab = "Post-test probability (negative test)",
         main = "NPV: Model 1 vs Selected Model 3")
    first <- FALSE
  } else {
    lines(sub$Pretest_Probability, sub$NPV,
          col = colors[i], lwd = lwd[i], lty = lty[i])
  }
  
  # Confidence band
  polygon(c(sub$Pretest_Probability,
            rev(sub$Pretest_Probability)),
          c(sub$NPV_low,
            rev(sub$NPV_high)),
          col = adjustcolor(colors[i], alpha.f = 0.15),
          border = NA)
  
  # Annotate
  for (xv in vlines) {
    yv <- approx(sub$Pretest_Probability, sub$NPV, xout = xv)$y
    text(xv, yv, labels = round(yv,2),
         pos = 4, col = colors[i], cex = 0.8)
  }
  
  i <- i + 1
}

abline(v = vlines, lty = 2)

legend("topright",
       legend = valid_models,
       col = colors,
       lwd = lwd,
       lty = lty,
       bty = "n",
       cex = 0.8)

dev.off()

#############################################################################################################################################

####################################
# Continuous NRI Calculation (RACE)
####################################

library(dplyr)
library(purrr)
library(nricens)
library(tibble)
library(ggplot2)

# ------------------- Function to compute continuous NRI -------------------
compute_nri_ci <- function(model1, model3, data, nboot = 1000) {
  if ("glm" %in% class(model1)) model1 <- list(model = model1)
  if ("glm" %in% class(model3)) model3 <- list(model = model3)
  
  vars1 <- all.vars(model1$model$formula)
  vars3 <- all.vars(model3$model$formula)
  
  tmp <- data %>% select(all_of(unique(c(vars1, vars3)))) %>% na.omit()
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  p1 <- predict(model1$model, tmp, type = "response")
  p3 <- predict(model3$model, tmp, type = "response")
  
  nri_out <- nribin(
    event = tmp$Group_bin,
    p.std = p1,
    p.new = p3,
    updown = "diff",
    cut = 0,
    niter = nboot
  )
  
  tibble(
    Model3 = paste(vars3[vars3 != "Group"], collapse = " + "),
    NRI_event     = nri_out$nri["NRI+", "Estimate"],
    NRI_event_LCL = nri_out$nri["NRI+", "Lower"],
    NRI_event_UCL = nri_out$nri["NRI+", "Upper"],
    NRI_nonevent     = nri_out$nri["NRI-", "Estimate"],
    NRI_nonevent_LCL = nri_out$nri["NRI-", "Lower"],
    NRI_nonevent_UCL = nri_out$nri["NRI-", "Upper"]
  )
}

# ------------------- Compute NRI for all Model3 variants -------------------
nri_results <- map_df(model3_list, ~compute_nri_ci(model1, .x, model_data))
nri_results

write.csv(nri_results, "NRI_continuous_RACE_vs_Model3_YNHH.csv", row.names = FALSE)

# ------------------- Function to get scatter plot data -------------------
get_scatter_data <- function(model1, model3, data) {
  vars1 <- all.vars(model1$model$formula)
  vars3 <- all.vars(model3$model$formula)
  
  tmp <- data %>% select(all_of(unique(c(vars1, vars3)))) %>% na.omit()
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  tibble(
    Model3 = paste(vars3[vars3 != "Group"], collapse = " + "),
    ID = 1:nrow(tmp),
    Group = ifelse(tmp$Group_bin == 1, "Case", "Control"),
    Prob_Model1 = predict(model1$model, tmp, type = "response"),
    Prob_Model3 = predict(model3$model, tmp, type = "response")
  )
}

# ------------------- Prepare scatter data -------------------
scatter_data <- map_df(model3_list, ~get_scatter_data(model1, .x, model_data))

# ------------------- Generate and save scatter plots -------------------
unique_models <- unique(scatter_data$Model3)

for (mod in unique_models) {
  plot_data <- scatter_data %>% filter(Model3 == mod)
  
  p <- ggplot(plot_data, aes(x = Prob_Model1, y = Prob_Model3, color = Group)) +
    geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.01, height = 0.01)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Case" = "red", "Control" = "black")) +
    labs(
      title = paste0("Continuous NRI Scatter:\n", mod),
      x = "Predicted probability Model 1",
      y = "Predicted probability Model 3",
      color = "Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10)),
      plot.margin = unit(c(1, 1, 1.5, 1), "cm"),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  # Save as PDF
  pdf_filename <- paste0(
    "Continuous_NRI_Scatter_",
    gsub("[^A-Za-z0-9]+", "_", mod),
    ".pdf"
  )
  
  ggsave(pdf_filename,
         plot = p,
         width = 6,
         height = 6.5,
         device = "pdf")
}

##############################################################################################################################################

##############################################
# Integrated Discrimination Improvement (IDI)
# Model 3 vs Model 1 (YNHH, RACE)
##############################################

library(dplyr)
library(purrr)
library(tibble)

# ------------------- Function to compute IDI -------------------
compute_idi_ci <- function(model1, model3, data, B = 1000) {
  
  # Ensure consistent structure
  if ("glm" %in% class(model1)) model1 <- list(model = model1)
  if ("glm" %in% class(model3)) model3 <- list(model = model3)
  
  # Extract variables
  vars1 <- all.vars(model1$model$formula)
  vars3 <- all.vars(model3$model$formula)
  
  # Prepare dataset
  tmp <- data %>%
    select(all_of(unique(c(vars1, vars3)))) %>%
    na.omit()
  
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  # Predicted probabilities
  p1 <- predict(model1$model, tmp, type = "response")
  p3 <- predict(model3$model, tmp, type = "response")
  
  # ------------------- IDI calculation -------------------
  mean_case_1 <- mean(p1[tmp$Group_bin == 1])
  mean_ctrl_1 <- mean(p1[tmp$Group_bin == 0])
  
  mean_case_3 <- mean(p3[tmp$Group_bin == 1])
  mean_ctrl_3 <- mean(p3[tmp$Group_bin == 0])
  
  idi <- (mean_case_3 - mean_case_1) - (mean_ctrl_3 - mean_ctrl_1)
  
  # ------------------- Bootstrap -------------------
  n <- nrow(tmp)
  
  boot_idi <- replicate(B, {
    
    idx <- sample(seq_len(n), replace = TRUE)
    d   <- tmp[idx, ]
    
    p1b <- predict(model1$model, d, type = "response")
    p3b <- predict(model3$model, d, type = "response")
    
    yb <- d$Group_bin
    
    mc1 <- mean(p1b[yb == 1])
    m0c1 <- mean(p1b[yb == 0])
    
    mc3 <- mean(p3b[yb == 1])
    m0c3 <- mean(p3b[yb == 0])
    
    (mc3 - mc1) - (m0c3 - m0c1)
  })
  
  # Confidence interval
  ci <- quantile(boot_idi, c(0.025, 0.975), na.rm = TRUE)
  
  # Two-sided bootstrap p-value
  p_val <- 2 * min(
    mean(boot_idi <= 0, na.rm = TRUE),
    mean(boot_idi >= 0, na.rm = TRUE)
  )
  
  # ------------------- Output -------------------
  tibble(
    Model3 = paste(vars3[vars3 != "Group"], collapse = " + "),
    
    Mean_case_Model1 = mean_case_1,
    Mean_control_Model1 = mean_ctrl_1,
    
    Mean_case_Model3 = mean_case_3,
    Mean_control_Model3 = mean_ctrl_3,
    
    IDI = idi,
    IDI_CI_lower = ci[1],
    IDI_CI_upper = ci[2],
    IDI_p_value = p_val
  )
}

# ------------------- Run IDI for all Model 3 -------------------
idi_results <- map_df(model3_list, ~compute_idi_ci(model1, .x, model_data))

# ------------------- View results -------------------
print(idi_results)

# ------------------- Save to CSV -------------------
write.csv(
  idi_results,
  "IDI_Model1_vs_Model3_YNHH_RACE.csv",
  row.names = FALSE
)

##############################################################################################################################################

##############################################################################################################################################
#RACE cut off (5 or more)
##############################################################################################################################################

####################################################################
# RACE_cutoff: Prediction models for LVO (AUC + DeLong comparisons)
####################################################################

library(dplyr)
library(purrr)
library(pROC)
library(tidyr)

# ------------------- Significant proteins -------------------
sig_proteins <- results %>%
  filter(Covariate == "RACE_cutoff", p_value < 0.05) %>%
  pull(Protein)

# ------------------- Model data -------------------
model_data <- df %>%
  select(all_of(c("Group", "RACE_cutoff", sig_proteins)))

model_data$Group <- factor(model_data$Group, levels = c("noLVO","LVO"))

# ============================================================
# Helper: ROC + AUC
# ============================================================
run_model_auc <- function(formula_text, data) {
  
  vars <- all.vars(as.formula(formula_text))
  tmp <- data %>% select(all_of(vars)) %>% na.omit()
  
  if(length(unique(tmp$Group)) < 2){
    return(list(model = NULL, roc = NULL, auc = NA, ci = c(NA,NA,NA)))
  }
  
  fit <- glm(as.formula(formula_text), data = tmp, family = binomial)
  preds <- predict(fit, type = "response")
  
  roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"))
  
  list(
    model = fit,
    roc   = roc_obj,
    auc   = auc(roc_obj),
    ci    = ci.auc(roc_obj)
  )
}

# ============================================================
# Model 1 (reference)
# ============================================================
model1 <- run_model_auc("Group ~ RACE_cutoff", model_data)

# ============================================================
# Model 2 (single protein + RACE_cutoff)
# ============================================================
model2_list <- map(sig_proteins, function(prot) {
  run_model_auc(paste0("Group ~ RACE_cutoff + `", prot, "`"), model_data)
})
names(model2_list) <- sig_proteins

# ============================================================
# Model 3 (two proteins + RACE_cutoff)
# ============================================================
model3_list <- list()

if(length(sig_proteins) >= 2){
  
  combos <- combn(sig_proteins, 2, simplify = FALSE)
  
  for(combo in combos){
    
    formula3 <- paste("Group ~ RACE_cutoff +", paste0("`", combo, "`", collapse = " + "))
    vars <- c("Group", "RACE_cutoff", combo)
    
    tmp <- model_data %>% select(all_of(vars)) %>% na.omit()
    
    if(length(unique(tmp$Group)) < 2) next
    
    fit <- glm(as.formula(formula3), data = tmp, family = binomial)
    
    coefs <- summary(fit)$coefficients
    
    prot_rows <- sapply(combo, function(p) grep(p, rownames(coefs), fixed = TRUE))
    protein_p <- coefs[prot_rows, "Pr(>|z|)"]
    
    if(all(!is.na(protein_p)) && all(protein_p < 0.05)){
      
      preds <- predict(fit, type = "response")
      roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"))
      
      model3_list[[paste(combo, collapse = "_")]] <- list(
        model = fit,
        roc   = roc_obj,
        auc   = auc(roc_obj),
        ci    = ci.auc(roc_obj)
      )
    }
  }
}

# ============================================================
# Base AUC table (Model 1)
# ============================================================
auc_results <- tibble(
  Model = "RACE_cutoff",
  AUC = as.numeric(model1$auc),
  CI_lower = as.numeric(model1$ci[1]),
  CI_upper = as.numeric(model1$ci[3])
)

# ============================================================
# Model 2 table + DeLong vs Model 1
# ============================================================
if(length(model2_list) > 0){
  
  model2_df <- tibble(
    Model = paste0(names(model2_list), " + RACE_cutoff"),
    AUC = sapply(model2_list, function(x) as.numeric(x$auc)),
    CI_lower = sapply(model2_list, function(x) as.numeric(x$ci[1])),
    CI_upper = sapply(model2_list, function(x) as.numeric(x$ci[3])),
    
    Delong_p_value = sapply(model2_list, function(x) {
      if(is.null(model1$roc) || is.null(x$roc)) return(NA_real_)
      roc.test(model1$roc, x$roc, method = "delong")$p.value
    })
  )
  
  auc_results <- bind_rows(auc_results, model2_df)
}

# ============================================================
# Model 3 table + DeLong vs Model 1
# ============================================================
if(length(model3_list) > 0){
  
  model3_df <- tibble(
    Model = paste0(names(model3_list), " + RACE_cutoff"),
    AUC = sapply(model3_list, function(x) as.numeric(x$auc)),
    CI_lower = sapply(model3_list, function(x) as.numeric(x$ci[1])),
    CI_upper = sapply(model3_list, function(x) as.numeric(x$ci[3])),
    
    Delong_p_value = sapply(model3_list, function(x) {
      if(is.null(model1$roc) || is.null(x$roc)) return(NA_real_)
      roc.test(model1$roc, x$roc, method = "delong")$p.value
    })
  )
  
  auc_results <- bind_rows(auc_results, model3_df)
}

# ============================================================
# Extract OR, 95% CI, and p-values from all models
# and combine with AUC + DeLong results
# ============================================================

# ============================================================
# Helper function: Extract OR results
# ============================================================
extract_or_results <- function(fit, model_name){
  
  if(is.null(fit)){
    return(NULL)
  }
  
  coefs <- summary(fit)$coefficients
  
  # Remove intercept
  coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
  
  # Calculate OR and 95% CI
  OR  <- exp(coefs[, "Estimate"])
  LCL <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
  UCL <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
  
  # Extract p-values
  pvals <- coefs[, "Pr(>|z|)"]
  
  tibble(
    Model = model_name,
    Variable = rownames(coefs),
    OR = OR,
    CI_lower_OR = LCL,
    CI_upper_OR = UCL,
    P_value = pvals
  )
}

# ============================================================
# Extract OR results: Model 1
# ============================================================
or_results_m1 <- extract_or_results(
  model1$model,
  "RACE_cutoff"
)

# ============================================================
# Extract OR results: Model 2
# ============================================================
or_results_m2 <- map2(
  model2_list,
  names(model2_list),
  ~ extract_or_results(
    .x$model,
    paste0(.y, " + RACE_cutoff")
  )
)

or_results_m2 <- bind_rows(or_results_m2)

# ============================================================
# Extract OR results: Model 3
# ============================================================
if(length(model3_list) > 0){
  
  or_results_m3 <- map2(
    model3_list,
    names(model3_list),
    ~ extract_or_results(
      .x$model,
      paste0(.y, " + RACE_cutoff")
    )
  )
  
  or_results_m3 <- bind_rows(or_results_m3)
  
} else {
  
  or_results_m3 <- tibble()
}

# ============================================================
# Combine all OR results
# ============================================================
all_or_results <- bind_rows(
  or_results_m1,
  or_results_m2,
  or_results_m3
)

# ============================================================
# Merge OR results with AUC results
# ============================================================
final_results <- all_or_results %>%
  left_join(
    auc_results,
    by = "Model"
  )

# ============================================================
# View final combined results
# ============================================================
View(final_results)

# ============================================================
# Save final results
# ============================================================
write.csv(
  final_results,
  "LVO_prediction_models_OR_AUC_Delong_YNHH_RACE_cutoff.csv",
  row.names = FALSE
)

##############################################################################################################################################

##############################################################################################################################################
#G-FAST stroke scale
##############################################################################################################################################

###############################################################
# G-FAST: Prediction models for LVO (AUC + DeLong comparisons)
###############################################################

library(dplyr)
library(purrr)
library(pROC)
library(tidyr)

# ------------------- Significant proteins -------------------
sig_proteins <- results %>%
  filter(Covariate == "G_FAST", p_value < 0.05) %>%
  pull(Protein)

# ------------------- Model data -------------------
model_data <- df %>%
  select(all_of(c("Group", "G_FAST", sig_proteins)))

model_data$Group <- factor(model_data$Group, levels = c("noLVO","LVO"))

# ============================================================
# Helper: ROC + AUC
# ============================================================
run_model_auc <- function(formula_text, data) {
  
  vars <- all.vars(as.formula(formula_text))
  tmp <- data %>% select(all_of(vars)) %>% na.omit()
  
  if(length(unique(tmp$Group)) < 2){
    return(list(model = NULL, roc = NULL, auc = NA, ci = c(NA,NA,NA)))
  }
  
  fit <- glm(as.formula(formula_text), data = tmp, family = binomial)
  preds <- predict(fit, type = "response")
  
  roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"))
  
  list(
    model = fit,
    roc   = roc_obj,
    auc   = auc(roc_obj),
    ci    = ci.auc(roc_obj)
  )
}

# ============================================================
# Model 1 (reference)
# ============================================================
model1 <- run_model_auc("Group ~ G_FAST", model_data)

# ============================================================
# Model 2 (single protein + G-FAST)
# ============================================================
model2_list <- map(sig_proteins, function(prot) {
  run_model_auc(paste0("Group ~ G_FAST + `", prot, "`"), model_data)
})
names(model2_list) <- sig_proteins

# ============================================================
# Model 3 (two proteins + G-FAST)
# ============================================================
model3_list <- list()

if(length(sig_proteins) >= 2){
  
  combos <- combn(sig_proteins, 2, simplify = FALSE)
  
  for(combo in combos){
    
    formula3 <- paste("Group ~ G_FAST +", paste0("`", combo, "`", collapse = " + "))
    vars <- c("Group", "G_FAST", combo)
    
    tmp <- model_data %>% select(all_of(vars)) %>% na.omit()
    
    if(length(unique(tmp$Group)) < 2) next
    
    fit <- glm(as.formula(formula3), data = tmp, family = binomial)
    
    coefs <- summary(fit)$coefficients
    
    prot_rows <- sapply(combo, function(p) grep(p, rownames(coefs), fixed = TRUE))
    protein_p <- coefs[prot_rows, "Pr(>|z|)"]
    
    if(all(!is.na(protein_p)) && all(protein_p < 0.05)){
      
      preds <- predict(fit, type = "response")
      roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"))
      
      model3_list[[paste(combo, collapse = "_")]] <- list(
        model = fit,
        roc   = roc_obj,
        auc   = auc(roc_obj),
        ci    = ci.auc(roc_obj)
      )
    }
  }
}

# ============================================================
# Base AUC table (Model 1)
# ============================================================
auc_results <- tibble(
  Model = "G_FAST",
  AUC = as.numeric(model1$auc),
  CI_lower = as.numeric(model1$ci[1]),
  CI_upper = as.numeric(model1$ci[3])
)

# ============================================================
# Model 2 table + DeLong vs Model 1
# ============================================================
if(length(model2_list) > 0){
  
  model2_df <- tibble(
    Model = paste0(names(model2_list), " + G_FAST"),
    AUC = sapply(model2_list, function(x) as.numeric(x$auc)),
    CI_lower = sapply(model2_list, function(x) as.numeric(x$ci[1])),
    CI_upper = sapply(model2_list, function(x) as.numeric(x$ci[3])),
    
    Delong_p_value = sapply(model2_list, function(x) {
      if(is.null(model1$roc) || is.null(x$roc)) return(NA_real_)
      roc.test(model1$roc, x$roc, method = "delong")$p.value
    })
  )
  
  auc_results <- bind_rows(auc_results, model2_df)
}

# ============================================================
# Model 3 table + DeLong vs Model 1
# ============================================================
if(length(model3_list) > 0){
  
  model3_df <- tibble(
    Model = paste0(names(model3_list), " + G_FAST"),
    AUC = sapply(model3_list, function(x) as.numeric(x$auc)),
    CI_lower = sapply(model3_list, function(x) as.numeric(x$ci[1])),
    CI_upper = sapply(model3_list, function(x) as.numeric(x$ci[3])),
    
    Delong_p_value = sapply(model3_list, function(x) {
      if(is.null(model1$roc) || is.null(x$roc)) return(NA_real_)
      roc.test(model1$roc, x$roc, method = "delong")$p.value
    })
  )
  
  auc_results <- bind_rows(auc_results, model3_df)
}

# ============================================================
# Extract OR, 95% CI, and p-values from all models
# and combine with AUC + DeLong results
# ============================================================

# ============================================================
# Helper function: Extract OR results
# ============================================================
extract_or_results <- function(fit, model_name){
  
  if(is.null(fit)){
    return(NULL)
  }
  
  coefs <- summary(fit)$coefficients
  
  # Remove intercept
  coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
  
  # Calculate OR and 95% CI
  OR  <- exp(coefs[, "Estimate"])
  LCL <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
  UCL <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
  
  # Extract p-values
  pvals <- coefs[, "Pr(>|z|)"]
  
  tibble(
    Model = model_name,
    Variable = rownames(coefs),
    OR = OR,
    CI_lower_OR = LCL,
    CI_upper_OR = UCL,
    P_value = pvals
  )
}

# ============================================================
# Extract OR results: Model 1
# ============================================================
or_results_m1 <- extract_or_results(
  model1$model,
  "G_FAST"
)

# ============================================================
# Extract OR results: Model 2
# ============================================================
or_results_m2 <- map2(
  model2_list,
  names(model2_list),
  ~ extract_or_results(
    .x$model,
    paste0(.y, " + G_FAST")
  )
)

or_results_m2 <- bind_rows(or_results_m2)

# ============================================================
# Extract OR results: Model 3
# ============================================================
if(length(model3_list) > 0){
  
  or_results_m3 <- map2(
    model3_list,
    names(model3_list),
    ~ extract_or_results(
      .x$model,
      paste0(.y, " + G_FAST")
    )
  )
  
  or_results_m3 <- bind_rows(or_results_m3)
  
} else {
  
  or_results_m3 <- tibble()
}

# ============================================================
# Combine all OR results
# ============================================================
all_or_results <- bind_rows(
  or_results_m1,
  or_results_m2,
  or_results_m3
)

# ============================================================
# Merge OR results with AUC results
# ============================================================
final_results <- all_or_results %>%
  left_join(
    auc_results,
    by = "Model"
  )

# ============================================================
# View final combined results
# ============================================================
View(final_results)

# ============================================================
# Save final results
# ============================================================
write.csv(
  final_results,
  "LVO_prediction_models_OR_AUC_Delong_YNHH_G-FAST.csv",
  row.names = FALSE
)

# ============================================================
# ROC overlay
# ============================================================

library(dplyr)
library(pROC)

selected_pairs <- c(
  "GH2.10978-39_ACP2.9237-54",
  "CD2.7100-31_ACP2.9237-54"
)

valid_pairs <- selected_pairs[selected_pairs %in% names(model3_list)]

cols <- c("blue", "red", "darkgreen")

pdf("ROC_Overlay_Model1_vs_Selected_pairs_Model3_G_FAST.pdf",
    width = 6.5, height = 6.5)

plot(model1$roc,
     col = cols[1],
     lwd = 2,
     legacy.axes = TRUE,
     main = "ROC Overlay: Model 1 vs Selected Model 3 Pairs (G_FAST)")

i <- 2

for (pair in valid_pairs) {
  
  m3 <- model3_list[[pair]]
  
  if (is.null(m3$roc)) next
  
  plot(m3$roc,
       col = cols[i],
       lwd = 2,
       add = TRUE)
  
  i <- i + 1
}

abline(a = 0, b = 1, lty = 2, col = "grey")

legend_labels <- c(
  paste0("Model 1: G_FAST (AUC = ", round(model1$auc, 3), ")"),
  paste0(valid_pairs[1], " + G_FAST (AUC = ", round(model3_list[[valid_pairs[1]]]$auc, 3), ")"),
  paste0(valid_pairs[2], " + G_FAST (AUC = ", round(model3_list[[valid_pairs[2]]]$auc, 3), ")")
)

legend("bottomright",
       legend = legend_labels,
       col = cols,
       lwd = 2,
       bty = "n")

dev.off()

##########################################################################################################################################

###############################################################################
# Diagnostic Metrics for Model 1 and Model 3 using fixed cutoff (YNHH, G-FAST)
###############################################################################

library(dplyr)

# ------------------- Bootstrap sensitivity & specificity -------------------
bootstrap_sens_spec <- function(fit, data, cutoff, B = 1000) {
  
  vars <- all.vars(formula(fit))
  tmp  <- data %>% select(all_of(vars)) %>% na.omit()
  n    <- nrow(tmp)
  
  boot <- replicate(B, {
    
    idx <- sample(seq_len(n), replace = TRUE)
    d   <- tmp[idx, ]
    
    preds <- predict(fit, newdata = d, type = "response")
    predc <- ifelse(preds >= cutoff, "LVO", "noLVO")
    
    cm <- table(
      Predicted = factor(predc, levels = c("LVO","noLVO")),
      Actual    = factor(d$Group, levels = c("LVO","noLVO"))
    )
    
    TP <- cm["LVO","LVO"]
    FP <- cm["LVO","noLVO"]
    TN <- cm["noLVO","noLVO"]
    FN <- cm["noLVO","LVO"]
    
    c(
      sens = ifelse(TP + FN > 0, TP / (TP + FN), NA),
      spec = ifelse(TN + FP > 0, TN / (TN + FP), NA)
    )
  })
  
  boot
}

# ------------------- Extract metrics -------------------
extract_sens_spec <- function(fit, data, cutoff, B = 1000) {
  
  vars <- all.vars(formula(fit))
  tmp  <- data %>% select(all_of(vars)) %>% na.omit()
  
  preds <- predict(fit, newdata = tmp, type = "response")
  predc <- ifelse(preds >= cutoff, "LVO", "noLVO")
  
  cm <- table(
    Predicted = factor(predc, levels = c("LVO","noLVO")),
    Actual    = factor(tmp$Group, levels = c("LVO","noLVO"))
  )
  
  TP <- cm["LVO","LVO"]
  FP <- cm["LVO","noLVO"]
  TN <- cm["noLVO","noLVO"]
  FN <- cm["noLVO","LVO"]
  
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  
  boot <- bootstrap_sens_spec(fit, data, cutoff, B)
  
  sens_ci <- quantile(boot["sens", ], c(0.025, 0.975), na.rm = TRUE)
  spec_ci <- quantile(boot["spec", ], c(0.025, 0.975), na.rm = TRUE)
  
  tibble(
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    Sensitivity = sens,
    Sensitivity_CI_lower = sens_ci[1],
    Sensitivity_CI_upper = sens_ci[2],
    Specificity = spec,
    Specificity_CI_lower = spec_ci[1],
    Specificity_CI_upper = spec_ci[2],
    LR_positive = sens / (1 - spec),
    LR_negative = (1 - sens) / spec
  )
}

# ------------------- PPV / NPV -------------------
compute_posttest <- function(sens, spec, pretest = seq(0.05, 0.95, by = 0.05)) {
  
  tibble(
    Pretest_Probability = pretest,
    PPV = (sens * pretest) /
      (sens * pretest + (1 - spec) * (1 - pretest)),
    NPV = (spec * (1 - pretest)) /
      ((1 - sens) * pretest + spec * (1 - pretest))
  )
}

bootstrap_posttest <- function(boot_mat, pretest_grid) {
  
  ppv_mat <- sapply(pretest_grid, function(p) {
    (boot_mat["sens", ] * p) /
      (boot_mat["sens", ] * p + (1 - boot_mat["spec", ]) * (1 - p))
  })
  
  npv_mat <- sapply(pretest_grid, function(p) {
    (boot_mat["spec", ] * (1 - p)) /
      ((1 - boot_mat["sens", ]) * p + boot_mat["spec", ] * (1 - p))
  })
  
  tibble(
    Pretest_Probability = pretest_grid,
    PPV_low  = apply(ppv_mat, 2, quantile, 0.025, na.rm = TRUE),
    PPV_high = apply(ppv_mat, 2, quantile, 0.975, na.rm = TRUE),
    NPV_low  = apply(npv_mat, 2, quantile, 0.025, na.rm = TRUE),
    NPV_high = apply(npv_mat, 2, quantile, 0.975, na.rm = TRUE)
  )
}

# ============================================================
# APPLY TO MODEL 1 + MODEL 3 with fixed manual cutoff
# ============================================================

pretest_models <- list(
  Model1_G_FAST = model1
)

if (length(model3_list) > 0) {
  for (nm in names(model3_list)) {
    pretest_models[[paste0("Model3_", nm)]] <- model3_list[[nm]]
  }
}

metrics_list <- list()
curves_list  <- list()
bands_list   <- list()

pretest_grid <- seq(0.05, 0.95, by = 0.05)

for (nm in names(pretest_models)) {
  
  obj <- pretest_models[[nm]]
  fit <- obj$model
  
  if (is.null(fit)) next
  
  cutoff <- 0.3
  
  core <- extract_sens_spec(
    fit = fit,
    data = model_data,
    cutoff = cutoff,
    B = 1000
  ) %>%
    mutate(Model = nm,
           Cutoff = cutoff)
  
  curves <- compute_posttest(
    sens = core$Sensitivity,
    spec = core$Specificity,
    pretest = pretest_grid
  ) %>%
    mutate(Model = nm)
  
  boot <- bootstrap_sens_spec(
    fit,
    model_data,
    cutoff = cutoff,
    B = 1000
  )
  
  bands <- bootstrap_posttest(boot, pretest_grid) %>%
    mutate(Model = nm)
  
  metrics_list[[nm]] <- core
  curves_list[[nm]]  <- curves
  bands_list[[nm]]   <- bands
}

metrics_df <- bind_rows(metrics_list)
curves_df  <- bind_rows(curves_list)
bands_df   <- bind_rows(bands_list)

ppv_npv_combined_df <- curves_df %>%
  left_join(bands_df, by = c("Model", "Pretest_Probability"))

# ------------------- SAVE -------------------

write.csv(
  metrics_df,
  "Diagnostic_metrics_Model1_Model3_Prob0.3_G_FAST_YNHH.csv",
  row.names = FALSE
)

write.csv(
  ppv_npv_combined_df,
  "PPV_NPV_Model1_Model3_Prob0.3_G_FAST_YNHH.csv",
  row.names = FALSE
)

##########################################################################################################################################

##################################################################################################################
# Overlay PPV & NPV Curves with bootstrap bands and reference lines for Model 1 and select Model 3 (YNHH, G_FAST)
##################################################################################################################

library(dplyr)

# -------- Select models --------
selected_pairs <- c(
  "GH2.10978-39_ACP2.9237-54",
  "CD2.7100-31_ACP2.9237-54"
)

# Match naming used in your pipeline
model_names <- c(
  "Model1_G_FAST",
  paste0("Model3_", selected_pairs)
)

valid_models <- model_names[model_names %in% unique(ppv_npv_combined_df$Model)]

# -------- Subset data --------
plot_df <- ppv_npv_combined_df %>%
  filter(Model %in% valid_models)

# -------- Plot settings --------
colors <- c("black", "red", "blue")[seq_along(valid_models)]
lty    <- c(2, rep(1, length(valid_models)-1))
lwd    <- c(2, rep(3, length(valid_models)-1))
vlines <- c(0.25, 0.5)

# ================================
# -------- PPV OVERLAY ----------
# ================================
pdf("PPV_Overlay_Model1_vs_Model3_G_FAST.pdf", width = 6, height = 6)

first <- TRUE
i <- 1

for (m in valid_models) {
  
  sub <- plot_df %>% filter(Model == m)
  
  if (first) {
    plot(sub$Pretest_Probability, sub$PPV,
         type = "l",
         col = colors[i], lwd = lwd[i], lty = lty[i],
         ylim = c(0,1),
         xlab = "Pretest probability",
         ylab = "Post-test probability (positive test)",
         main = "PPV: Model 1 vs Selected Model 3")
    first <- FALSE
  } else {
    lines(sub$Pretest_Probability, sub$PPV,
          col = colors[i], lwd = lwd[i], lty = lty[i])
  }
  
  polygon(c(sub$Pretest_Probability,
            rev(sub$Pretest_Probability)),
          c(sub$PPV_low,
            rev(sub$PPV_high)),
          col = adjustcolor(colors[i], alpha.f = 0.15),
          border = NA)
  
  for (xv in vlines) {
    yv <- approx(sub$Pretest_Probability, sub$PPV, xout = xv)$y
    text(xv, yv, labels = round(yv,2),
         pos = 4, col = colors[i], cex = 0.8)
  }
  
  i <- i + 1
}

abline(v = vlines, lty = 2)

legend("bottomright",
       legend = valid_models,
       col = colors,
       lwd = lwd,
       lty = lty,
       bty = "n",
       cex = 0.8)

dev.off()

# ================================
# -------- NPV OVERLAY ----------
# ================================
pdf("NPV_Overlay_Model1_vs_Model3_G_FAST.pdf", width = 6, height = 6)

first <- TRUE
i <- 1

for (m in valid_models) {
  
  sub <- plot_df %>% filter(Model == m)
  
  if (first) {
    plot(sub$Pretest_Probability, sub$NPV,
         type = "l",
         col = colors[i], lwd = lwd[i], lty = lty[i],
         ylim = c(0,1),
         xlab = "Pretest probability",
         ylab = "Post-test probability (negative test)",
         main = "NPV: Model 1 vs Selected Model 3")
    first <- FALSE
  } else {
    lines(sub$Pretest_Probability, sub$NPV,
          col = colors[i], lwd = lwd[i], lty = lty[i])
  }
  
  polygon(c(sub$Pretest_Probability,
            rev(sub$Pretest_Probability)),
          c(sub$NPV_low,
            rev(sub$NPV_high)),
          col = adjustcolor(colors[i], alpha.f = 0.15),
          border = NA)
  
  for (xv in vlines) {
    yv <- approx(sub$Pretest_Probability, sub$NPV, xout = xv)$y
    text(xv, yv, labels = round(yv,2),
         pos = 4, col = colors[i], cex = 0.8)
  }
  
  i <- i + 1
}

abline(v = vlines, lty = 2)

legend("topright",
       legend = valid_models,
       col = colors,
       lwd = lwd,
       lty = lty,
       bty = "n",
       cex = 0.8)

dev.off()

#############################################################################################################################################

######################################
# Continuous NRI Calculation (G-FAST)
######################################

library(dplyr)
library(purrr)
library(nricens)
library(tibble)
library(ggplot2)

# ------------------- Function to compute continuous NRI -------------------
compute_nri_ci <- function(model1, model3, data, nboot = 1000) {
  if ("glm" %in% class(model1)) model1 <- list(model = model1)
  if ("glm" %in% class(model3)) model3 <- list(model = model3)
  
  vars1 <- all.vars(model1$model$formula)
  vars3 <- all.vars(model3$model$formula)
  
  tmp <- data %>% select(all_of(unique(c(vars1, vars3)))) %>% na.omit()
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  p1 <- predict(model1$model, tmp, type = "response")
  p3 <- predict(model3$model, tmp, type = "response")
  
  nri_out <- nribin(
    event = tmp$Group_bin,
    p.std = p1,
    p.new = p3,
    updown = "diff",
    cut = 0,
    niter = nboot
  )
  
  tibble(
    Model3 = paste(vars3[vars3 != "Group"], collapse = " + "),
    NRI_event     = nri_out$nri["NRI+", "Estimate"],
    NRI_event_LCL = nri_out$nri["NRI+", "Lower"],
    NRI_event_UCL = nri_out$nri["NRI+", "Upper"],
    NRI_nonevent     = nri_out$nri["NRI-", "Estimate"],
    NRI_nonevent_LCL = nri_out$nri["NRI-", "Lower"],
    NRI_nonevent_UCL = nri_out$nri["NRI-", "Upper"]
  )
}

# ------------------- Compute NRI for all Model3 variants -------------------
nri_results <- map_df(model3_list, ~compute_nri_ci(model1, .x, model_data))
nri_results

write.csv(nri_results, "NRI_continuous_G_FAST_vs_Model3_YNHH.csv", row.names = FALSE)

# ------------------- Function to get scatter plot data -------------------
get_scatter_data <- function(model1, model3, data) {
  vars1 <- all.vars(model1$model$formula)
  vars3 <- all.vars(model3$model$formula)
  
  tmp <- data %>% select(all_of(unique(c(vars1, vars3)))) %>% na.omit()
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  tibble(
    Model3 = paste(vars3[vars3 != "Group"], collapse = " + "),
    ID = 1:nrow(tmp),
    Group = ifelse(tmp$Group_bin == 1, "Case", "Control"),
    Prob_Model1 = predict(model1$model, tmp, type = "response"),
    Prob_Model3 = predict(model3$model, tmp, type = "response")
  )
}

# ------------------- Prepare scatter data -------------------
scatter_data <- map_df(model3_list, ~get_scatter_data(model1, .x, model_data))

# ------------------- Generate and save scatter plots -------------------
unique_models <- unique(scatter_data$Model3)

for (mod in unique_models) {
  plot_data <- scatter_data %>% filter(Model3 == mod)
  
  p <- ggplot(plot_data, aes(x = Prob_Model1, y = Prob_Model3, color = Group)) +
    geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.01, height = 0.01)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Case" = "red", "Control" = "black")) +
    labs(
      title = paste0("Continuous NRI Scatter:\n", mod),
      x = "Predicted probability Model 1",
      y = "Predicted probability Model 3",
      color = "Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10)),
      plot.margin = unit(c(1, 1, 1.5, 1), "cm"),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  # Save as PDF
  pdf_filename <- paste0(
    "Continuous_NRI_Scatter_",
    gsub("[^A-Za-z0-9]+", "_", mod),
    ".pdf"
  )
  
  ggsave(pdf_filename,
         plot = p,
         width = 6,
         height = 6.5,
         device = "pdf")
}

##############################################################################################################################################

##############################################
# Integrated Discrimination Improvement (IDI)
# Model 3 vs Model 1 (YNHH, G-FAST)
##############################################

library(dplyr)
library(purrr)
library(tibble)

# ------------------- Function to compute IDI -------------------
compute_idi_ci <- function(model1, model3, data, B = 1000) {
  
  # Ensure consistent structure
  if ("glm" %in% class(model1)) model1 <- list(model = model1)
  if ("glm" %in% class(model3)) model3 <- list(model = model3)
  
  # Extract variables
  vars1 <- all.vars(model1$model$formula)
  vars3 <- all.vars(model3$model$formula)
  
  # Prepare dataset
  tmp <- data %>%
    select(all_of(unique(c(vars1, vars3)))) %>%
    na.omit()
  
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  # Predicted probabilities
  p1 <- predict(model1$model, tmp, type = "response")
  p3 <- predict(model3$model, tmp, type = "response")
  
  # ------------------- IDI calculation -------------------
  mean_case_1 <- mean(p1[tmp$Group_bin == 1])
  mean_ctrl_1 <- mean(p1[tmp$Group_bin == 0])
  
  mean_case_3 <- mean(p3[tmp$Group_bin == 1])
  mean_ctrl_3 <- mean(p3[tmp$Group_bin == 0])
  
  idi <- (mean_case_3 - mean_case_1) - (mean_ctrl_3 - mean_ctrl_1)
  
  # ------------------- Bootstrap -------------------
  n <- nrow(tmp)
  
  boot_idi <- replicate(B, {
    
    idx <- sample(seq_len(n), replace = TRUE)
    d   <- tmp[idx, ]
    
    p1b <- predict(model1$model, d, type = "response")
    p3b <- predict(model3$model, d, type = "response")
    
    yb <- d$Group_bin
    
    mc1 <- mean(p1b[yb == 1])
    m0c1 <- mean(p1b[yb == 0])
    
    mc3 <- mean(p3b[yb == 1])
    m0c3 <- mean(p3b[yb == 0])
    
    (mc3 - mc1) - (m0c3 - m0c1)
  })
  
  # Confidence interval
  ci <- quantile(boot_idi, c(0.025, 0.975), na.rm = TRUE)
  
  # Two-sided bootstrap p-value
  p_val <- 2 * min(
    mean(boot_idi <= 0, na.rm = TRUE),
    mean(boot_idi >= 0, na.rm = TRUE)
  )
  
  # ------------------- Output -------------------
  tibble(
    Model3 = paste(vars3[vars3 != "Group"], collapse = " + "),
    
    Mean_case_Model1 = mean_case_1,
    Mean_control_Model1 = mean_ctrl_1,
    
    Mean_case_Model3 = mean_case_3,
    Mean_control_Model3 = mean_ctrl_3,
    
    IDI = idi,
    IDI_CI_lower = ci[1],
    IDI_CI_upper = ci[2],
    IDI_p_value = p_val
  )
}

# ------------------- Run IDI for all Model 3 -------------------
idi_results <- map_df(model3_list, ~compute_idi_ci(model1, .x, model_data))

# ------------------- View results -------------------
print(idi_results)

# ------------------- Save to CSV -------------------
write.csv(
  idi_results,
  "IDI_Model1_vs_Model3_YNHH_G_FAST.csv",
  row.names = FALSE
)

##############################################################################################################################################

##############################################################################################################################################
#mG-FAST stroke scale
##############################################################################################################################################

################################################################
# mG-FAST: Prediction models for LVO (AUC + DeLong comparisons)
################################################################

library(dplyr)
library(purrr)
library(pROC)
library(tidyr)

# ------------------- Significant proteins -------------------
sig_proteins <- results %>%
  filter(Covariate == "mG_FAST", p_value < 0.05) %>%
  pull(Protein)

# ------------------- Model data -------------------
model_data <- df %>%
  select(all_of(c("Group", "mG_FAST", sig_proteins)))

model_data$Group <- factor(model_data$Group, levels = c("noLVO","LVO"))

# ============================================================
# Helper: ROC + AUC
# ============================================================
run_model_auc <- function(formula_text, data) {
  
  vars <- all.vars(as.formula(formula_text))
  tmp <- data %>% select(all_of(vars)) %>% na.omit()
  
  if(length(unique(tmp$Group)) < 2){
    return(list(model = NULL, roc = NULL, auc = NA, ci = c(NA,NA,NA)))
  }
  
  fit <- glm(as.formula(formula_text), data = tmp, family = binomial)
  preds <- predict(fit, type = "response")
  
  roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"))
  
  list(
    model = fit,
    roc   = roc_obj,
    auc   = auc(roc_obj),
    ci    = ci.auc(roc_obj)
  )
}

# ============================================================
# Model 1 (reference)
# ============================================================
model1 <- run_model_auc("Group ~ mG_FAST", model_data)

# ============================================================
# Model 2 (single protein + mG_FAST)
# ============================================================
model2_list <- map(sig_proteins, function(prot) {
  run_model_auc(paste0("Group ~ mG_FAST + `", prot, "`"), model_data)
})
names(model2_list) <- sig_proteins

# ============================================================
# Model 3 (two proteins + mG_FAST)
# ============================================================
model3_list <- list()

if(length(sig_proteins) >= 2){
  
  combos <- combn(sig_proteins, 2, simplify = FALSE)
  
  for(combo in combos){
    
    formula3 <- paste("Group ~ mG_FAST +", paste0("`", combo, "`", collapse = " + "))
    vars <- c("Group", "mG_FAST", combo)
    
    tmp <- model_data %>% select(all_of(vars)) %>% na.omit()
    
    if(length(unique(tmp$Group)) < 2) next
    
    fit <- glm(as.formula(formula3), data = tmp, family = binomial)
    
    coefs <- summary(fit)$coefficients
    
    prot_rows <- sapply(combo, function(p) grep(p, rownames(coefs), fixed = TRUE))
    protein_p <- coefs[prot_rows, "Pr(>|z|)"]
    
    if(all(!is.na(protein_p)) && all(protein_p < 0.05)){
      
      preds <- predict(fit, type = "response")
      roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"))
      
      model3_list[[paste(combo, collapse = "_")]] <- list(
        model = fit,
        roc   = roc_obj,
        auc   = auc(roc_obj),
        ci    = ci.auc(roc_obj)
      )
    }
  }
}

# ============================================================
# Base AUC table (Model 1)
# ============================================================
auc_results <- tibble(
  Model = "mG_FAST",
  AUC = as.numeric(model1$auc),
  CI_lower = as.numeric(model1$ci[1]),
  CI_upper = as.numeric(model1$ci[3])
)

# ============================================================
# Model 2 table + DeLong vs Model 1
# ============================================================
if(length(model2_list) > 0){
  
  model2_df <- tibble(
    Model = paste0(names(model2_list), " + mG_FAST"),
    AUC = sapply(model2_list, function(x) as.numeric(x$auc)),
    CI_lower = sapply(model2_list, function(x) as.numeric(x$ci[1])),
    CI_upper = sapply(model2_list, function(x) as.numeric(x$ci[3])),
    
    Delong_p_value = sapply(model2_list, function(x) {
      if(is.null(model1$roc) || is.null(x$roc)) return(NA_real_)
      roc.test(model1$roc, x$roc, method = "delong")$p.value
    })
  )
  
  auc_results <- bind_rows(auc_results, model2_df)
}

# ============================================================
# Model 3 table + DeLong vs Model 1
# ============================================================
if(length(model3_list) > 0){
  
  model3_df <- tibble(
    Model = paste0(names(model3_list), " + mG_FAST"),
    AUC = sapply(model3_list, function(x) as.numeric(x$auc)),
    CI_lower = sapply(model3_list, function(x) as.numeric(x$ci[1])),
    CI_upper = sapply(model3_list, function(x) as.numeric(x$ci[3])),
    
    Delong_p_value = sapply(model3_list, function(x) {
      if(is.null(model1$roc) || is.null(x$roc)) return(NA_real_)
      roc.test(model1$roc, x$roc, method = "delong")$p.value
    })
  )
  
  auc_results <- bind_rows(auc_results, model3_df)
}

# ============================================================
# Extract OR, 95% CI, and p-values from all models
# and combine with AUC + DeLong results
# ============================================================

# ============================================================
# Helper function: Extract OR results
# ============================================================
extract_or_results <- function(fit, model_name){
  
  if(is.null(fit)){
    return(NULL)
  }
  
  coefs <- summary(fit)$coefficients
  
  # Remove intercept
  coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
  
  # Calculate OR and 95% CI
  OR  <- exp(coefs[, "Estimate"])
  LCL <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
  UCL <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
  
  # Extract p-values
  pvals <- coefs[, "Pr(>|z|)"]
  
  tibble(
    Model = model_name,
    Variable = rownames(coefs),
    OR = OR,
    CI_lower_OR = LCL,
    CI_upper_OR = UCL,
    P_value = pvals
  )
}

# ============================================================
# Extract OR results: Model 1
# ============================================================
or_results_m1 <- extract_or_results(
  model1$model,
  "mG_FAST"
)

# ============================================================
# Extract OR results: Model 2
# ============================================================
or_results_m2 <- map2(
  model2_list,
  names(model2_list),
  ~ extract_or_results(
    .x$model,
    paste0(.y, " + mG_FAST")
  )
)

or_results_m2 <- bind_rows(or_results_m2)

# ============================================================
# Extract OR results: Model 3
# ============================================================
if(length(model3_list) > 0){
  
  or_results_m3 <- map2(
    model3_list,
    names(model3_list),
    ~ extract_or_results(
      .x$model,
      paste0(.y, " + mG_FAST")
    )
  )
  
  or_results_m3 <- bind_rows(or_results_m3)
  
} else {
  
  or_results_m3 <- tibble()
}

# ============================================================
# Combine all OR results
# ============================================================
all_or_results <- bind_rows(
  or_results_m1,
  or_results_m2,
  or_results_m3
)

# ============================================================
# Merge OR results with AUC results
# ============================================================
final_results <- all_or_results %>%
  left_join(
    auc_results,
    by = "Model"
  )

# ============================================================
# View final combined results
# ============================================================
View(final_results)

# ============================================================
# Save final results
# ============================================================
write.csv(
  final_results,
  "LVO_prediction_models_OR_AUC_Delong_YNHH_mG-FAST.csv",
  row.names = FALSE
)

# ============================================================
# ROC overlay
# ============================================================

library(dplyr)
library(pROC)

selected_pairs <- c(
  "GH2.10978-39_ACP2.9237-54",
  "CD2.7100-31_ACP2.9237-54"
)

valid_pairs <- selected_pairs[selected_pairs %in% names(model3_list)]

cols <- c("blue", "red", "darkgreen")

pdf("ROC_Overlay_Model1_vs_Selected_pairs_Model3_mG_FAST.pdf",
    width = 6.5, height = 6.5)

plot(model1$roc,
     col = cols[1],
     lwd = 2,
     legacy.axes = TRUE,
     main = "ROC Overlay: Model 1 vs Selected Model 3 Pairs (mG_FAST)")

i <- 2

for (pair in valid_pairs) {
  
  m3 <- model3_list[[pair]]
  
  if (is.null(m3$roc)) next
  
  plot(m3$roc,
       col = cols[i],
       lwd = 2,
       add = TRUE)
  
  i <- i + 1
}

abline(a = 0, b = 1, lty = 2, col = "grey")

legend_labels <- c(
  paste0("Model 1: mG_FAST (AUC = ", round(model1$auc, 3), ")"),
  paste0(valid_pairs[1], " + mG_FAST (AUC = ", round(model3_list[[valid_pairs[1]]]$auc, 3), ")"),
  paste0(valid_pairs[2], " + mG_FAST (AUC = ", round(model3_list[[valid_pairs[2]]]$auc, 3), ")")
)

legend("bottomright",
       legend = legend_labels,
       col = cols,
       lwd = 2,
       bty = "n")

dev.off()

##########################################################################################################################################

################################################################################
# Diagnostic Metrics for Model 1 and Model 3 using fixed cutoff (YNHH, mG_FAST)
################################################################################

library(dplyr)

# ------------------- Bootstrap sensitivity & specificity -------------------
bootstrap_sens_spec <- function(fit, data, cutoff, B = 1000) {
  
  vars <- all.vars(formula(fit))
  tmp  <- data %>% select(all_of(vars)) %>% na.omit()
  n    <- nrow(tmp)
  
  boot <- replicate(B, {
    
    idx <- sample(seq_len(n), replace = TRUE)
    d   <- tmp[idx, ]
    
    preds <- predict(fit, newdata = d, type = "response")
    predc <- ifelse(preds >= cutoff, "LVO", "noLVO")
    
    cm <- table(
      Predicted = factor(predc, levels = c("LVO","noLVO")),
      Actual    = factor(d$Group, levels = c("LVO","noLVO"))
    )
    
    TP <- cm["LVO","LVO"]
    FP <- cm["LVO","noLVO"]
    TN <- cm["noLVO","noLVO"]
    FN <- cm["noLVO","LVO"]
    
    c(
      sens = ifelse(TP + FN > 0, TP / (TP + FN), NA),
      spec = ifelse(TN + FP > 0, TN / (TN + FP), NA)
    )
  })
  
  boot
}

# ------------------- Extract metrics -------------------
extract_sens_spec <- function(fit, data, cutoff, B = 1000) {
  
  vars <- all.vars(formula(fit))
  tmp  <- data %>% select(all_of(vars)) %>% na.omit()
  
  preds <- predict(fit, newdata = tmp, type = "response")
  predc <- ifelse(preds >= cutoff, "LVO", "noLVO")
  
  cm <- table(
    Predicted = factor(predc, levels = c("LVO","noLVO")),
    Actual    = factor(tmp$Group, levels = c("LVO","noLVO"))
  )
  
  TP <- cm["LVO","LVO"]
  FP <- cm["LVO","noLVO"]
  TN <- cm["noLVO","noLVO"]
  FN <- cm["noLVO","LVO"]
  
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  
  boot <- bootstrap_sens_spec(fit, data, cutoff, B)
  
  sens_ci <- quantile(boot["sens", ], c(0.025, 0.975), na.rm = TRUE)
  spec_ci <- quantile(boot["spec", ], c(0.025, 0.975), na.rm = TRUE)
  
  tibble(
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    Sensitivity = sens,
    Sensitivity_CI_lower = sens_ci[1],
    Sensitivity_CI_upper = sens_ci[2],
    Specificity = spec,
    Specificity_CI_lower = spec_ci[1],
    Specificity_CI_upper = spec_ci[2],
    LR_positive = sens / (1 - spec),
    LR_negative = (1 - sens) / spec
  )
}

# ------------------- PPV / NPV -------------------
compute_posttest <- function(sens, spec, pretest = seq(0.05, 0.95, by = 0.05)) {
  
  tibble(
    Pretest_Probability = pretest,
    PPV = (sens * pretest) /
      (sens * pretest + (1 - spec) * (1 - pretest)),
    NPV = (spec * (1 - pretest)) /
      ((1 - sens) * pretest + spec * (1 - pretest))
  )
}

bootstrap_posttest <- function(boot_mat, pretest_grid) {
  
  ppv_mat <- sapply(pretest_grid, function(p) {
    (boot_mat["sens", ] * p) /
      (boot_mat["sens", ] * p + (1 - boot_mat["spec", ]) * (1 - p))
  })
  
  npv_mat <- sapply(pretest_grid, function(p) {
    (boot_mat["spec", ] * (1 - p)) /
      ((1 - boot_mat["sens", ]) * p + boot_mat["spec", ] * (1 - p))
  })
  
  tibble(
    Pretest_Probability = pretest_grid,
    PPV_low  = apply(ppv_mat, 2, quantile, 0.025, na.rm = TRUE),
    PPV_high = apply(ppv_mat, 2, quantile, 0.975, na.rm = TRUE),
    NPV_low  = apply(npv_mat, 2, quantile, 0.025, na.rm = TRUE),
    NPV_high = apply(npv_mat, 2, quantile, 0.975, na.rm = TRUE)
  )
}

# ============================================================
# APPLY TO MODEL 1 + MODEL 3 with fixed manual cutoff
# ============================================================

pretest_models <- list(
  Model1_mG_FAST = model1
)

if (length(model3_list) > 0) {
  for (nm in names(model3_list)) {
    pretest_models[[paste0("Model3_", nm)]] <- model3_list[[nm]]
  }
}

metrics_list <- list()
curves_list  <- list()
bands_list   <- list()

pretest_grid <- seq(0.05, 0.95, by = 0.05)

for (nm in names(pretest_models)) {
  
  obj <- pretest_models[[nm]]
  fit <- obj$model
  
  if (is.null(fit)) next
  
  cutoff <- 0.3
  
  core <- extract_sens_spec(
    fit = fit,
    data = model_data,
    cutoff = cutoff,
    B = 1000
  ) %>%
    mutate(Model = nm,
           Cutoff = cutoff)
  
  curves <- compute_posttest(
    sens = core$Sensitivity,
    spec = core$Specificity,
    pretest = pretest_grid
  ) %>%
    mutate(Model = nm)
  
  boot <- bootstrap_sens_spec(
    fit,
    model_data,
    cutoff = cutoff,
    B = 1000
  )
  
  bands <- bootstrap_posttest(boot, pretest_grid) %>%
    mutate(Model = nm)
  
  metrics_list[[nm]] <- core
  curves_list[[nm]]  <- curves
  bands_list[[nm]]   <- bands
}

metrics_df <- bind_rows(metrics_list)
curves_df  <- bind_rows(curves_list)
bands_df   <- bind_rows(bands_list)

ppv_npv_combined_df <- curves_df %>%
  left_join(bands_df, by = c("Model", "Pretest_Probability"))

# ------------------- SAVE -------------------

write.csv(
  metrics_df,
  "Diagnostic_metrics_Model1_Model3_Prob0.3_mG_FAST_YNHH.csv",
  row.names = FALSE
)

write.csv(
  ppv_npv_combined_df,
  "PPV_NPV_Model1_Model3_Prob0.3_mG_FAST_YNHH.csv",
  row.names = FALSE
)  

##########################################################################################################################################

###################################################################################################################
# Overlay PPV & NPV Curves with bootstrap bands and reference lines for Model 1 and select Model 3 (YNHH, mG_FAST)
###################################################################################################################

library(dplyr)

# -------- Select models --------
selected_pairs <- c(
  "GH2.10978-39_ACP2.9237-54",
  "CD2.7100-31_ACP2.9237-54"
)

# Match naming used in your pipeline
model_names <- c(
  "Model1_mG_FAST",
  paste0("Model3_", selected_pairs)
)

valid_models <- model_names[model_names %in% unique(ppv_npv_combined_df$Model)]

# -------- Subset data --------
plot_df <- ppv_npv_combined_df %>%
  filter(Model %in% valid_models)

# -------- Plot settings --------
colors <- c("black", "red", "blue")[seq_along(valid_models)]
lty    <- c(2, rep(1, length(valid_models)-1))
lwd    <- c(2, rep(3, length(valid_models)-1))
vlines <- c(0.25, 0.5)

# ================================
# -------- PPV OVERLAY ----------
# ================================
pdf("PPV_Overlay_Model1_vs_Model3_mG_FAST.pdf", width = 6, height = 6)

first <- TRUE
i <- 1

for (m in valid_models) {
  
  sub <- plot_df %>% filter(Model == m)
  
  if (first) {
    plot(sub$Pretest_Probability, sub$PPV,
         type = "l",
         col = colors[i], lwd = lwd[i], lty = lty[i],
         ylim = c(0,1),
         xlab = "Pretest probability",
         ylab = "Post-test probability (positive test)",
         main = "PPV: Model 1 vs Selected Model 3")
    first <- FALSE
  } else {
    lines(sub$Pretest_Probability, sub$PPV,
          col = colors[i], lwd = lwd[i], lty = lty[i])
  }
  
  polygon(c(sub$Pretest_Probability,
            rev(sub$Pretest_Probability)),
          c(sub$PPV_low,
            rev(sub$PPV_high)),
          col = adjustcolor(colors[i], alpha.f = 0.15),
          border = NA)
  
  for (xv in vlines) {
    yv <- approx(sub$Pretest_Probability, sub$PPV, xout = xv)$y
    text(xv, yv, labels = round(yv,2),
         pos = 4, col = colors[i], cex = 0.8)
  }
  
  i <- i + 1
}

abline(v = vlines, lty = 2)

legend("bottomright",
       legend = valid_models,
       col = colors,
       lwd = lwd,
       lty = lty,
       bty = "n",
       cex = 0.8)

dev.off()

# ================================
# -------- NPV OVERLAY ----------
# ================================
pdf("NPV_Overlay_Model1_vs_Model3_mG_FAST.pdf", width = 6, height = 6)

first <- TRUE
i <- 1

for (m in valid_models) {
  
  sub <- plot_df %>% filter(Model == m)
  
  if (first) {
    plot(sub$Pretest_Probability, sub$NPV,
         type = "l",
         col = colors[i], lwd = lwd[i], lty = lty[i],
         ylim = c(0,1),
         xlab = "Pretest probability",
         ylab = "Post-test probability (negative test)",
         main = "NPV: Model 1 vs Selected Model 3")
    first <- FALSE
  } else {
    lines(sub$Pretest_Probability, sub$NPV,
          col = colors[i], lwd = lwd[i], lty = lty[i])
  }
  
  polygon(c(sub$Pretest_Probability,
            rev(sub$Pretest_Probability)),
          c(sub$NPV_low,
            rev(sub$NPV_high)),
          col = adjustcolor(colors[i], alpha.f = 0.15),
          border = NA)
  
  for (xv in vlines) {
    yv <- approx(sub$Pretest_Probability, sub$NPV, xout = xv)$y
    text(xv, yv, labels = round(yv,2),
         pos = 4, col = colors[i], cex = 0.8)
  }
  
  i <- i + 1
}

abline(v = vlines, lty = 2)

legend("topright",
       legend = valid_models,
       col = colors,
       lwd = lwd,
       lty = lty,
       bty = "n",
       cex = 0.8)

dev.off()

#############################################################################################################################################

#######################################
# Continuous NRI Calculation (mG-FAST)
#######################################

library(dplyr)
library(purrr)
library(nricens)
library(tibble)
library(ggplot2)

# ------------------- Function to compute continuous NRI -------------------
compute_nri_ci <- function(model1, model3, data, nboot = 1000) {
  if ("glm" %in% class(model1)) model1 <- list(model = model1)
  if ("glm" %in% class(model3)) model3 <- list(model = model3)
  
  vars1 <- all.vars(model1$model$formula)
  vars3 <- all.vars(model3$model$formula)
  
  tmp <- data %>% select(all_of(unique(c(vars1, vars3)))) %>% na.omit()
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  p1 <- predict(model1$model, tmp, type = "response")
  p3 <- predict(model3$model, tmp, type = "response")
  
  nri_out <- nribin(
    event = tmp$Group_bin,
    p.std = p1,
    p.new = p3,
    updown = "diff",
    cut = 0,
    niter = nboot
  )
  
  tibble(
    Model3 = paste(vars3[vars3 != "Group"], collapse = " + "),
    NRI_event     = nri_out$nri["NRI+", "Estimate"],
    NRI_event_LCL = nri_out$nri["NRI+", "Lower"],
    NRI_event_UCL = nri_out$nri["NRI+", "Upper"],
    NRI_nonevent     = nri_out$nri["NRI-", "Estimate"],
    NRI_nonevent_LCL = nri_out$nri["NRI-", "Lower"],
    NRI_nonevent_UCL = nri_out$nri["NRI-", "Upper"]
  )
}

# ------------------- Compute NRI for all Model3 variants -------------------
nri_results <- map_df(model3_list, ~compute_nri_ci(model1, .x, model_data))
nri_results

write.csv(nri_results, "NRI_continuous_mG_FAST_vs_Model3_YNHH.csv", row.names = FALSE)

# ------------------- Function to get scatter plot data -------------------
get_scatter_data <- function(model1, model3, data) {
  vars1 <- all.vars(model1$model$formula)
  vars3 <- all.vars(model3$model$formula)
  
  tmp <- data %>% select(all_of(unique(c(vars1, vars3)))) %>% na.omit()
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  tibble(
    Model3 = paste(vars3[vars3 != "Group"], collapse = " + "),
    ID = 1:nrow(tmp),
    Group = ifelse(tmp$Group_bin == 1, "Case", "Control"),
    Prob_Model1 = predict(model1$model, tmp, type = "response"),
    Prob_Model3 = predict(model3$model, tmp, type = "response")
  )
}

# ------------------- Prepare scatter data -------------------
scatter_data <- map_df(model3_list, ~get_scatter_data(model1, .x, model_data))

# ------------------- Generate and save scatter plots -------------------
unique_models <- unique(scatter_data$Model3)

for (mod in unique_models) {
  plot_data <- scatter_data %>% filter(Model3 == mod)
  
  p <- ggplot(plot_data, aes(x = Prob_Model1, y = Prob_Model3, color = Group)) +
    geom_point(alpha = 0.7, size = 2, position = position_jitter(width = 0.01, height = 0.01)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Case" = "red", "Control" = "black")) +
    labs(
      title = paste0("Continuous NRI Scatter:\n", mod),
      x = "Predicted probability Model 1",
      y = "Predicted probability Model 3",
      color = "Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10)),
      plot.margin = unit(c(1, 1, 1.5, 1), "cm"),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  # Save as PDF
  pdf_filename <- paste0(
    "Continuous_NRI_Scatter_",
    gsub("[^A-Za-z0-9]+", "_", mod),
    ".pdf"
  )
  
  ggsave(pdf_filename,
         plot = p,
         width = 6,
         height = 6.5,
         device = "pdf")
}

##############################################################################################################################################

##############################################
# Integrated Discrimination Improvement (IDI)
# Model 3 vs Model 1 (YNHH, mG-FAST)
##############################################

library(dplyr)
library(purrr)
library(tibble)

# ------------------- Function to compute IDI -------------------
compute_idi_ci <- function(model1, model3, data, B = 1000) {
  
  # Ensure consistent structure
  if ("glm" %in% class(model1)) model1 <- list(model = model1)
  if ("glm" %in% class(model3)) model3 <- list(model = model3)
  
  # Extract variables
  vars1 <- all.vars(model1$model$formula)
  vars3 <- all.vars(model3$model$formula)
  
  # Prepare dataset
  tmp <- data %>%
    select(all_of(unique(c(vars1, vars3)))) %>%
    na.omit()
  
  tmp$Group_bin <- ifelse(tmp$Group == "LVO", 1, 0)
  
  # Predicted probabilities
  p1 <- predict(model1$model, tmp, type = "response")
  p3 <- predict(model3$model, tmp, type = "response")
  
  # ------------------- IDI calculation -------------------
  mean_case_1 <- mean(p1[tmp$Group_bin == 1])
  mean_ctrl_1 <- mean(p1[tmp$Group_bin == 0])
  
  mean_case_3 <- mean(p3[tmp$Group_bin == 1])
  mean_ctrl_3 <- mean(p3[tmp$Group_bin == 0])
  
  idi <- (mean_case_3 - mean_case_1) - (mean_ctrl_3 - mean_ctrl_1)
  
  # ------------------- Bootstrap -------------------
  n <- nrow(tmp)
  
  boot_idi <- replicate(B, {
    
    idx <- sample(seq_len(n), replace = TRUE)
    d   <- tmp[idx, ]
    
    p1b <- predict(model1$model, d, type = "response")
    p3b <- predict(model3$model, d, type = "response")
    
    yb <- d$Group_bin
    
    mc1 <- mean(p1b[yb == 1])
    m0c1 <- mean(p1b[yb == 0])
    
    mc3 <- mean(p3b[yb == 1])
    m0c3 <- mean(p3b[yb == 0])
    
    (mc3 - mc1) - (m0c3 - m0c1)
  })
  
  # Confidence interval
  ci <- quantile(boot_idi, c(0.025, 0.975), na.rm = TRUE)
  
  # Two-sided bootstrap p-value
  p_val <- 2 * min(
    mean(boot_idi <= 0, na.rm = TRUE),
    mean(boot_idi >= 0, na.rm = TRUE)
  )
  
  # ------------------- Output -------------------
  tibble(
    Model3 = paste(vars3[vars3 != "Group"], collapse = " + "),
    
    Mean_case_Model1 = mean_case_1,
    Mean_control_Model1 = mean_ctrl_1,
    
    Mean_case_Model3 = mean_case_3,
    Mean_control_Model3 = mean_ctrl_3,
    
    IDI = idi,
    IDI_CI_lower = ci[1],
    IDI_CI_upper = ci[2],
    IDI_p_value = p_val
  )
}

# ------------------- Run IDI for all Model 3 -------------------
idi_results <- map_df(model3_list, ~compute_idi_ci(model1, .x, model_data))

# ------------------- View results -------------------
print(idi_results)

# ------------------- Save to CSV -------------------
write.csv(
  idi_results,
  "IDI_Model1_vs_Model3_YNHH_mG_FAST.csv",
  row.names = FALSE
)

#############################################################################################################################################

#####################################################################################################################################
# Sub-analysis of Categorical NIHSS score cut off at 10 for YNHH cohort
#####################################################################################################################################

#############################################################################################
# Bivariable Logistic Regression analysis of Concordant DAPs Categorical NIHSS (YNHH cohort)
#############################################################################################

library(readxl)
library(dplyr)
library(purrr)
library(openxlsx)

# ------------ USER INPUTS ------------
input_file  <- "YNHH_LVO_data.xlsx"
sheet       <- 8

group_col   <- "Group"

meta_vars <- c("NIHSS10")

protein_start_col <- 6
# -------------------------------------

df <- read_excel(input_file, sheet = sheet)
df[[group_col]] <- factor(df[[group_col]], levels = c("noLVO", "LVO"))

protein_cols <- names(df)[protein_start_col:ncol(df)]
protein_cols

# ---------- MODEL FUNCTION ----------
run_model <- function(data, protein, covariate, outcome = "Group") {
  
  tmp <- data %>%
    select(all_of(c(outcome, protein, covariate))) %>%
    na.omit()
  
  if (length(unique(tmp[[outcome]])) < 2)
    return(c(OR=NA, CI_low=NA, CI_high=NA, p=NA))
  
  # Wrap variable names in backticks for safe formula
  formula_text <- paste0("`", outcome, "` ~ `", covariate, "` + `", protein, "`")
  fit <- glm(
    formula = as.formula(formula_text),
    data = tmp,
    family = binomial()
  )
  
  s <- summary(fit)$coefficients
  protein_row <- paste0("`", protein, "`")
  
  if (!(protein_row %in% rownames(s)))
    return(c(OR=NA, CI_low=NA, CI_high=NA, p=NA))
  
  est <- s[protein_row, "Estimate"]
  se  <- s[protein_row, "Std. Error"]
  p   <- s[protein_row, "Pr(>|z|)"]
  
  OR      <- exp(est)
  CI_low  <- exp(est - 1.96 * se)
  CI_high <- exp(est + 1.96 * se)
  
  c(OR=OR, CI_low=CI_low, CI_high=CI_high, p=p)
}

# -------- RUN MODELS ACROSS ALL PROTEINS × COVARIATES --------

results <- map_df(protein_cols, function(prot) {
  map_df(meta_vars, function(meta) {
    out <- run_model(df, protein = prot, covariate = meta)
    
    tibble(
      Protein   = prot,
      Covariate = meta,
      OR        = out["OR"],
      CI_low    = out["CI_low"],
      CI_high   = out["CI_high"],
      p_value   = out["p"]
    )
  })
})

# -------- SAVE AS CSV --------
write.csv(results, "Bivariable_logistic_NIHSS_cutoffs_YNHH.csv", row.names = FALSE)

################################################################################################################################################

################################################################
# NIHSS10: Prediction models for LVO (AUC + DeLong comparisons)
################################################################

library(dplyr)
library(purrr)
library(pROC)
library(tidyr)

# ------------------- Significant proteins -------------------
sig_proteins <- results %>%
  filter(Covariate == "NIHSS10", p_value < 0.05) %>%
  pull(Protein)

# ------------------- Model data -------------------
model_data <- df %>%
  select(all_of(c("Group", "NIHSS10", sig_proteins)))

model_data$Group <- factor(model_data$Group, levels = c("noLVO","LVO"))

# ============================================================
# Helper: ROC + AUC
# ============================================================
run_model_auc <- function(formula_text, data) {
  
  vars <- all.vars(as.formula(formula_text))
  tmp <- data %>% select(all_of(vars)) %>% na.omit()
  
  if(length(unique(tmp$Group)) < 2){
    return(list(model = NULL, roc = NULL, auc = NA, ci = c(NA,NA,NA)))
  }
  
  fit <- glm(as.formula(formula_text), data = tmp, family = binomial)
  preds <- predict(fit, type = "response")
  
  roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"))
  
  list(
    model = fit,
    roc   = roc_obj,
    auc   = auc(roc_obj),
    ci    = ci.auc(roc_obj)
  )
}

# ============================================================
# Model 1 (reference)
# ============================================================
model1 <- run_model_auc("Group ~ NIHSS10", model_data)

# ============================================================
# Model 2 (single protein + NIHSS10)
# ============================================================
model2_list <- map(sig_proteins, function(prot) {
  run_model_auc(paste0("Group ~ NIHSS10 + `", prot, "`"), model_data)
})
names(model2_list) <- sig_proteins

# ============================================================
# Model 3 (two proteins + NIHSS10)
# ============================================================
model3_list <- list()

if(length(sig_proteins) >= 2){
  
  combos <- combn(sig_proteins, 2, simplify = FALSE)
  
  for(combo in combos){
    
    formula3 <- paste("Group ~ NIHSS10 +", paste0("`", combo, "`", collapse = " + "))
    vars <- c("Group", "NIHSS10", combo)
    
    tmp <- model_data %>% select(all_of(vars)) %>% na.omit()
    
    if(length(unique(tmp$Group)) < 2) next
    
    fit <- glm(as.formula(formula3), data = tmp, family = binomial)
    
    coefs <- summary(fit)$coefficients
    
    prot_rows <- sapply(combo, function(p) grep(p, rownames(coefs), fixed = TRUE))
    protein_p <- coefs[prot_rows, "Pr(>|z|)"]
    
    if(all(!is.na(protein_p)) && all(protein_p < 0.05)){
      
      preds <- predict(fit, type = "response")
      roc_obj <- roc(tmp$Group, preds, levels = c("noLVO","LVO"))
      
      model3_list[[paste(combo, collapse = "_")]] <- list(
        model = fit,
        roc   = roc_obj,
        auc   = auc(roc_obj),
        ci    = ci.auc(roc_obj)
      )
    }
  }
}

# ============================================================
# Base AUC table (Model 1)
# ============================================================
auc_results <- tibble(
  Model = "NIHSS10",
  AUC = as.numeric(model1$auc),
  CI_lower = as.numeric(model1$ci[1]),
  CI_upper = as.numeric(model1$ci[3])
)

# ============================================================
# Model 2 table + DeLong vs Model 1
# ============================================================
if(length(model2_list) > 0){
  
  model2_df <- tibble(
    Model = paste0(names(model2_list), " + NIHSS10"),
    AUC = sapply(model2_list, function(x) as.numeric(x$auc)),
    CI_lower = sapply(model2_list, function(x) as.numeric(x$ci[1])),
    CI_upper = sapply(model2_list, function(x) as.numeric(x$ci[3])),
    
    Delong_p_value = sapply(model2_list, function(x) {
      if(is.null(model1$roc) || is.null(x$roc)) return(NA_real_)
      roc.test(model1$roc, x$roc, method = "delong")$p.value
    })
  )
  
  auc_results <- bind_rows(auc_results, model2_df)
}

# ============================================================
# Model 3 table + DeLong vs Model 1
# ============================================================
if(length(model3_list) > 0){
  
  model3_df <- tibble(
    Model = paste0(names(model3_list), " + NIHSS10"),
    AUC = sapply(model3_list, function(x) as.numeric(x$auc)),
    CI_lower = sapply(model3_list, function(x) as.numeric(x$ci[1])),
    CI_upper = sapply(model3_list, function(x) as.numeric(x$ci[3])),
    
    Delong_p_value = sapply(model3_list, function(x) {
      if(is.null(model1$roc) || is.null(x$roc)) return(NA_real_)
      roc.test(model1$roc, x$roc, method = "delong")$p.value
    })
  )
  
  auc_results <- bind_rows(auc_results, model3_df)
}

# ============================================================
# Extract OR, 95% CI, and p-values from all models
# and combine with AUC + DeLong results
# ============================================================

# ============================================================
# Helper function: Extract OR results
# ============================================================
extract_or_results <- function(fit, model_name){
  
  if(is.null(fit)){
    return(NULL)
  }
  
  coefs <- summary(fit)$coefficients
  
  # Remove intercept
  coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
  
  # Calculate OR and 95% CI
  OR  <- exp(coefs[, "Estimate"])
  LCL <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
  UCL <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
  
  # Extract p-values
  pvals <- coefs[, "Pr(>|z|)"]
  
  tibble(
    Model = model_name,
    Variable = rownames(coefs),
    OR = OR,
    CI_lower_OR = LCL,
    CI_upper_OR = UCL,
    P_value = pvals
  )
}

# ============================================================
# Extract OR results: Model 1
# ============================================================
or_results_m1 <- extract_or_results(
  model1$model,
  "NIHSS10"
)

# ============================================================
# Extract OR results: Model 2
# ============================================================
or_results_m2 <- map2(
  model2_list,
  names(model2_list),
  ~ extract_or_results(
    .x$model,
    paste0(.y, " + NIHSS10")
  )
)

or_results_m2 <- bind_rows(or_results_m2)

# ============================================================
# Extract OR results: Model 3
# ============================================================
if(length(model3_list) > 0){
  
  or_results_m3 <- map2(
    model3_list,
    names(model3_list),
    ~ extract_or_results(
      .x$model,
      paste0(.y, " + NIHSS10")
    )
  )
  
  or_results_m3 <- bind_rows(or_results_m3)
  
} else {
  
  or_results_m3 <- tibble()
}

# ============================================================
# Combine all OR results
# ============================================================
all_or_results <- bind_rows(
  or_results_m1,
  or_results_m2,
  or_results_m3
)

# ============================================================
# Merge OR results with AUC results
# ============================================================
final_results <- all_or_results %>%
  left_join(
    auc_results,
    by = "Model"
  )

# ============================================================
# View final combined results
# ============================================================
View(final_results)

# ============================================================
# Save final results
# ============================================================
write.csv(
  final_results,
  "LVO_prediction_models_OR_AUC_Delong_YNHH_NIHSS10.csv",
  row.names = FALSE
)

#########################################################################################################################################

##############################################################################################
# Repeat measure samples
##############################################################################################

##########################################################
# Paired t-tests for Before vs After protein measurements
##########################################################

# Load required package
library(dplyr)
library(readxl)

# ------------ USER INPUTS ------------
input_file  <- "YNHH_LVO_Data.xlsx"
sheet       <- 9
# -------------------------------------

data <- read_excel(input_file, sheet = sheet)

# ------------------------------------------------------------
# Define protein names
# ------------------------------------------------------------
proteins <- c("GH2", "CD2", "C1QL2", "ACP2")

# ------------------------------------------------------------
# Initialize results dataframe
# ------------------------------------------------------------
results <- data.frame()

# ------------------------------------------------------------
# Run paired t-tests
# ------------------------------------------------------------
for(prot in proteins){
  
  before_col <- paste0(prot, "_Before")
  after_col  <- paste0(prot, "_After")
  
  # Remove missing pairs
  temp_data <- data %>%
    select(all_of(c(before_col, after_col))) %>%
    na.omit()
  
  # Paired t-test
  ttest <- t.test(
    temp_data[[before_col]],
    temp_data[[after_col]],
    paired = TRUE
  )
  
  # Store results
  results <- rbind(
    results,
    data.frame(
      Protein = prot,
      Mean_Before = mean(temp_data[[before_col]], na.rm = TRUE),
      Mean_After = mean(temp_data[[after_col]], na.rm = TRUE),
      Mean_Difference = mean(temp_data[[after_col]] - temp_data[[before_col]], na.rm = TRUE),
      T_statistic = ttest$statistic,
      DF = ttest$parameter,
      P_value = ttest$p.value,
      CI_Lower = ttest$conf.int[1],
      CI_Upper = ttest$conf.int[2]
    )
  )
}

# ------------------------------------------------------------
# View results
# ------------------------------------------------------------
print(results)

# ------------------------------------------------------------
# Save results to CSV
# ------------------------------------------------------------
write.csv(results,
          "paired_ttest_repeat_samples_YNHH.csv",
          row.names = FALSE)

######################################################
# Optional:
# If to run tests separately for LVO = 1 and LVO = 0
######################################################

results_by_LVO <- data.frame()

for(group in c(0,1)){
  
  subgroup <- data %>% filter(LVO == group)
  
  for(prot in proteins){
    
    before_col <- paste0(prot, "_Before")
    after_col  <- paste0(prot, "_After")
    
    temp_data <- subgroup %>%
      select(all_of(c(before_col, after_col))) %>%
      na.omit()
    
    ttest <- t.test(
      temp_data[[before_col]],
      temp_data[[after_col]],
      paired = TRUE
    )
    
    results_by_LVO <- rbind(
      results_by_LVO,
      data.frame(
        LVO = group,
        Protein = prot,
        Mean_Before = mean(temp_data[[before_col]], na.rm = TRUE),
        Mean_After = mean(temp_data[[after_col]], na.rm = TRUE),
        Mean_Difference = mean(temp_data[[after_col]] - temp_data[[before_col]], na.rm = TRUE),
        T_statistic = ttest$statistic,
        DF = ttest$parameter,
        P_value = ttest$p.value,
        CI_Lower = ttest$conf.int[1],
        CI_Upper = ttest$conf.int[2]
      )
    )
  }
}

# Save subgroup results
write.csv(results_by_LVO,
          "paired_ttest_by_LVO_repeat_samples_YNHH.csv",
          row.names = FALSE)