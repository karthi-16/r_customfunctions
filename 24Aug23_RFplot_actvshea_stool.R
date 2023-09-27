# Env
setwd("/home/minionpak2/Desktop/my_edit/DSTCD_16Sanalysis_filtered/ACTvsHEL/STOOL")
library(data.table)
library(readxl)
library(vegan)
library(plotrix)
library(randomForest)
library(rfUtilities)
library(caret)
library(RColorBrewer)
library(dunn.test)
library(captioner)
library(flextable)
library(officer)
library(phyloseq)
library(tidyr)
library(scales)
library(ggplot2)
library(ROCR)
library(pROC)
load("/home/minionpak2/Desktop/my_edit/DSTCD_16Sanalysis_filtered/ACTvsHEL/STOOL/24Aug23_RFplot_actvshea_stool.RData")

# Metadata for baseline
fread("/home/minionpak2/Desktop/my_edit/DSTCD_16Sanalysis_filtered/ACTvsHEL/STOOL/metah_tab.tsv")->met

fread("/home/minionpak2/Desktop/my_edit/DSTCD_16Sanalysis_filtered/ACTvsHEL/STOOL/Analysis/l6/l6.tsv")->l6
class(l6)

# RF - Level 6----
l6->rfdat_l6
colnames(rfdat_l6)
# Transpose!
# transpose(rfdat_l6, keep.names="#OTU", make.names=1)->rfdat_l6
str(rfdat_l6)
rfdat_l6<- as.data.frame(rfdat_l6)

# Look at and plot non zero values for each OTU
otu_nonzero_counts_l6 <- apply(rfdat_l6, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts_l6, breaks=100, col="red", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")
# R function for removing rare OTU's
remove_rare_l6 <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

otu_table_rare_removed_l6 <- remove_rare_l6(table=rfdat_l6, cutoff_pro=0.05)
dim(otu_table_rare_removed_l6)
# 148 obs 107 variables
class(otu_table_rare_removed_l6)
as.data.frame(otu_table_rare_removed_l6)->otu_table_rare_removed_l6

# Renormalise data
rownames(otu_table_rare_removed_l6)<-otu_table_rare_removed_l6[,1]
rownames(otu_table_rare_removed_l6)

otu_table_rare_removed_l6[,2:107]<-sapply(otu_table_rare_removed_l6[,2:107], as.numeric)

otu_table_rare_removed_norm_l6 <- sweep(otu_table_rare_removed_l6[,-c(1)], 2, colSums(otu_table_rare_removed_l6[,-c(1)]) , '/')*100
# Scale it by z score
otu_table_scaled_l6 <- scale(otu_table_rare_removed_norm_l6, center = TRUE, scale = TRUE)  
otu_table_scaled_dse_l6<-data.frame(t(otu_table_scaled_l6))

# prep data for disease_background
met[-1,]->metar
metar[GROUP==3, disease_cat:="healthy_control"]
metar[GROUP==2, disease_cat:="active_cd"]
as.data.frame(metar)->metar
rownames(metar)<-metar[,1]

otu_table_scaled_dse_l6$dse <- metar[rownames(otu_table_scaled_dse_l6), "disease_cat"]  

as.factor(otu_table_scaled_dse_l6$dse)->otu_table_scaled_dse_l6$dse


boruta_otu_table_l6<-otu_table_rare_removed_l6[,-1]
t(boruta_otu_table_l6)-> boruta_l6


## boruta 
## making a separate data.frame for genus with colnames modified
genus_df <- as.data.frame(boruta_l6)
# Remove hyphens after 'g__' in column names
colnames(genus_df) <- gsub("g__", "g__temp", colnames(genus_df))
colnames(genus_df) <- sub("-", "", colnames(genus_df))
# Rename 'g__temp' back to 'g__'
colnames(genus_df) <- gsub("g__temp", "g__", colnames(genus_df))



# CLR transform
genus <- read.csv('/home/minionpak2/Desktop/my_edit/DSTCD_16Sanalysis_filtered/ACTvsHEL/STOOL/Analysis/L6/L6.tsv', 
                  header=1, check.names = F, row.names = 1, sep = '\t')

genus.c <- t(genus); eps <- 0.5
genus.c <- genus.c * (1 - rowSums(genus.c==0) * eps / rowSums(genus.c))
genus.c[genus.c == 0] <- eps
genus.c <- sweep(genus.c, 1, rowSums(genus.c), '/')
ls <- log(genus.c)
genus.c <- t(ls - rowMeans(ls))
genus.c <- genus.c[, !is.nan(colSums(genus.c))]
genus.c <- as.data.frame(genus.c)

# rf.map <- tibble::column_to_rownames(cloud.results, var = 'Row.names')
rf.ids <- intersect(colnames(genus.c), rownames(metar))
rd.ids <- sort(rf.ids)
genus.clr <- t(genus.c[, rf.ids])
dim(genus.clr)

set.seed(123)
boruta.dse <- Boruta(x=genus.c, y=otu_table_scaled_dse_l6$dse)
length(boruta.dse$finalDecision[boruta.dse$finalDecision == 'Confirmed']) #3
sig.taxa.dse <- lapply(X=names(boruta.dse$finalDecision[boruta.dse$finalDecision == 'Confirmed']),
        FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))
#sig.taxa.dse <- names(boruta.dse$finalDecision[boruta.dse$finalDecision == 'Confirmed'])

sig.taxa.dse
length(sig.taxa.dse)

# Plot the features as heatmap
sig.genera <- genus.c[, colnames(genus.c) %in% sig.taxa.dse]
dys.ids <- as.character(rownames(otu_table_scaled_dse_l6[otu_table_scaled_dse_l6$dse =='active_cd', ]))
healthylike.ids <- as.character(rownames(otu_table_scaled_dse_l6[otu_table_scaled_dse_l6$dse =='healthy_control', ]))

dys_healthy_order <- c(dys.ids, healthylike.ids)

sig.genera.clr <- sig.genera[dys_healthy_order, ]
fullmap2 <- otu_table_scaled_dse_l6[dys_healthy_order, ]

hm.meta <- data.frame(fullmap2[,'dse'], row.names = rownames(fullmap2), col=fullmap2$dse =='healthy_control')
hm.meta$col[hm.meta$col==TRUE] <- '#E69F00'
hm.meta$col[hm.meta$col==FALSE] <- '#009E73'

# Define the color names and corresponding colors
# color_names <- c("Healthy", "active_cd")  # Replace with your actual color names
# colors <- c("#E69F00", "#009E73")  # Replace with your actual color codes
png(paste0("./R objects/boruta.png"),  # create PNG for the heat map
    width = 15*150,                        # 5 x 300 pixels
    height = 6*150,
    res = 300,                              # 300 pixels per inch
    pointsize = 6)                          # smaller font size
#modifying the colanmes for plot
# Use gsub to extract the last character from "g__"
# modifying the featutres ids
features_ids<- colnames(sig.genera.clr)
# Extract the last word based on the presence of "g__" 
last_words_boruta <- sapply(features_ids, function(str) {
  if (grepl("g__", str)) {
    gsub(".*;g__(.)", "\\1", str)
  } else if (grepl("f__", str)) {
    gsub(".+\\b(f__[A-Za-z]+)\\b", "\\1", str)
  } else {
    "No match"  # Handle cases where neither "g__" nor "f__" is found
  }
})
# adding the features ids to the colnames
colnames(sig.genera.clr) <- last_words_boruta

sig<- sig.genera.clr
sig<- data.frame(dse_status= hm.meta$fullmap2....dse.., sig.genera.clr) 

long_data_boruta <- sig%>%
  pivot_longer(cols = -dse_status, 
               names_to = "feature_name", 
               values_to = "feature_value")

##
heatmap_plot_boruta <- ggplot(long_data_boruta, aes(x = dse_status, y = reorder(feature_name, -feature_value), fill = feature_value)) +
  geom_tile(color = "White") +
  scale_fill_gradient(low = "gray", high = "red") +
  labs(title = "Active_CD vs Healthy_control - Boruta Features", x = "Category", y = "Features") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Print the heatmap
png("./R objects/actvsheastool_boruta.png", width = 800, height = 600)
print(heatmap_plot_boruta)
dev.off()

# colnames(sig.genera.clr) <- extracted_names
# heatmap3(
#   x = t(sig.genera.clr),
#   showColDendro = FALSE,
#   showRowDendro = FALSE,
#   ColSideColors = cbind(`Active = Green`=hm.meta$col),
#   cexRow = 1,
#   margins = c(2, 6),
#   Colv = NA)
# # Create a custom legend
# legend("left", legend = color_names, fill = colors )
# dev.off()

# Random forest analysis with 10 fold cross validation (24Aug2023) ----
set.seed(123)
RF_models <- vector("list", 5)
# Loop through 5 runs of cross-validation
for (i in 1:5) {
  RF_models[[i]] <- train(
    x = genus.clr, 
    y = factor(metar$dse_status),
    method = "rf",
    trControl = trainControl(method = "cv", number = 10),
    tuneLength = 18,
    ntree = 1000,
    metric = "Accuracy"
  )
}

# Extract error rates for each run and store in a data frame
error_data <- data.frame()
for (i in 1:5) {
  error_data <- rbind(error_data, data.frame(mtry = RF_models[[i]]$results$mtry,
                                             error = 1 - RF_models[[i]]$results$Accuracy,
                                             run = paste("Run", i)))
}
as.data.table(error_data)->error_data

# Calculate the average error rate for each mtry value across runs
avg_error_data <- aggregate(error ~ mtry, data = error_data, FUN = mean)
as.data.table(avg_error_data)->avg_error_data
error_data[,avg:=avg_error_data[error_data,error,on="mtry"]]
dcast(error_data, mtry+avg ~ run,value.var="error")->kdat
apply(kdat[,c(3:7)], 1, mean)
apply(kdat[,c(3:7)], 1, sd)
min(apply(kdat[,c(3:7)], 1, mean))
## Therefore minimum error + sd is > 0.3200000+ 0.02291739 = 0.3429174, and optimum number of obs is 129

# Plotting the errors for all 5 runs with average line
 png("R objects/error_rf.png")
ggplot() +
  geom_line(data = error_data, aes(x = mtry, y = error, group = run, col = "gray"), alpha = 0.3) +  # Individual error lines
  geom_line(data = avg_error_data, aes(x = mtry, y = error), color = "black", size = 0.5) +  # Average line
  geom_point(data = error_data, aes(x = mtry, y = error, col= "gray")) +  # Individual points
  labs(x = "Number of Variables", y = "Error Rate",
       title = "Error Rate for Different mtry Values in Random Forest (10 Runs)") +
  #scale_color_manual(values = c("blue", "green", "purple", "orange", "pink")) +
  theme_minimal()+
  geom_vline(xintercept=129,col="green")+
  guides( colour= FALSE)

dev.off()

# Extract variable importance scores from each model
importance_scores <- lapply(RF_models, function(model) varImp(model$finalModel))
# Extract variable importance scores for the first model
importance_scores_first_model <- varImp(RF_models[[1]]$finalModel)

# Extract variable importance scores for the 62 variables
# extract the rownames
rownames<- rownames(importance_scores_first_model)
# Convert importance_scores to a data frame
importance_df <- as.data.table(importance_scores_first_model)
importance_df$row_ids <- rownames

# Sort importance scores in descending order and select the top 62 markers
top_129_markers <- head(importance_df[order(-importance_df$Overall), ], 129)
top_129_marker_names <- top_129_markers$row_ids




##### continue with boruta codes for extrcting sampels with this 129 features
# Train Random Forest with Top 62 Markers Only
# Extract indices of the top 62 markers
top <- as.factor(top_62_marker_names)
train_62 <- otu_table_scaled_dse_l6[, c(top, ncol(otu_table_scaled_dse_l6))]

# Sample index for train and test data
set.seed(121)
index_rf <- sample(1:nrow(train_62), size = 0.75 * nrow(train_62))
train <- train_62[index_rf, ]
test <- train_62[-index_rf, ]

# rf_model for train data (train)
train_rf_62 <- randomForest(
  x = train[, 1:(ncol(train) - 1)],
  y = train[, ncol(train)], 
  importance = TRUE
)

# Predicting the probabilities for train data
rf_pro <- predict(train_rf_62, type = "prob")
# predicting the probabilities for test data
rf_test <- predict(train_rf_62, test, type = "prob")

plot(train_rf_62)
# making ROC for the train set 
# Extract predicted probabilities for the positive class
p.train <- rf_pro[, "active_cd"]  # Replace with the actual label

# Convert outcome to a binary factor (0/1)
outcome <- as.factor(ifelse(train$dse == "active_cd", 0, 1))  # Replace with your outcome variable

# Calculate and plot the ROC curve
roc1 <- roc(outcome, p.train,
            ci = TRUE, boot.n = 100, ci.alpha = 0.9, stratified = FALSE,
            plot = TRUE, percent = TRUE, col = 2,
            auc.polygon = FALSE, max.auc.polygon = FALSE, grid = TRUE)

# Customize the plot
plot(roc1)
plot(roc1, col = "sandybrown", add = TRUE)
sens.ci <- ci.se(roc1, specificities = seq(0, 100, 5))
plot(sens.ci, type = "shape", col = "mediumturquoise")
plot(sens.ci, type = "bars",col ="cyan1")
plot(roc1, col = "sandybrown", add = TRUE)
legend("bottomright", c(paste("AUC=", round(roc1$ci[2], 2), "%"),
                        paste("95% CI:", round(roc1$ci[1], 2), "%-", round(roc1$ci[3], 2), "%")))

# making roc for the test test
p.train_t <- rf_test[, "active_cd"]  # Replace with the actual label

# Convert outcome to a binary factor (0/1)
outcome_t <- as.factor(ifelse(test$dse == "active_cd", 0, 1))  # Replace with your outcome variable

# Calculate and plot the ROC curve
roc1_t <- roc(outcome_t, p.train_t,
            ci = TRUE, boot.n = 100, ci.alpha = 0.9, stratified = FALSE,
            plot = TRUE, percent = TRUE, col = 2,
            auc.polygon = FALSE, max.auc.polygon = FALSE, grid = TRUE)

# Customize the plot
sens.ci_t <- ci.se(roc1_t, specificities = seq(0, 100, 5))
plot(roc1_t)
plot(sens.ci_t, type = "shape", col = "mediumturquoise")
plot(sens.ci_t, type = "bars", col="cyan1")
plot(roc1_t, col = "sandybrown", add = TRUE)
legend("bottomright", c(paste("AUC=", round(roc1_t$ci[2], 2), "%"),
                        paste("95% CI:", round(roc1_t$ci[1], 2), "%-", round(roc1_t$ci[3], 2), "%")))

# rf_model for original top 62 variables for rf plot
set.seed(121)
rf_62 <- randomForest(
  x = train_62[, 1:(ncol(train_62) - 1)],
  y = train_62[, ncol(train_62)], 
  importance = TRUE
)

# Plots for rf_ variable importance 
RF_dse_classify_imp_l6 <- as.data.frame(rf_62$importance)
RF_dse_classify_imp_l6$features <- rownames( RF_dse_classify_imp_l6 )
as.data.table(RF_dse_classify_imp_l6)->RF_dse_classify_imp_l6
RF_dse_classify_imp_sorted_l6 <- RF_dse_classify_imp_l6[,.SD[order(-MeanDecreaseAccuracy)]]
barplot(RF_dse_classify_imp_sorted_l6$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")

# Make files
write.csv(RF_dse_classify_imp_sorted_l6, file="Analysis/l6/RF_dse_l6_62.csv")
# edit the feature labels manually and do the clean up.
## clean up random forest plots (in linux) see ACTvsHEA_stool-25jul23.py

# import the csv as features_df


# fread("Analysis/l6/field_g_l6_62")->fie_g_l6_62
# features_df <- read.csv("/home/minionpak2/Desktop/my_edit/DSTCD_16Sanalysis_filtered/ACTvsHEL/STOOL/Analysis/l6/RF_dse_l6_62.csv")
# features_df<- as.data.table(features_df)
# RF_dse_classify_imp_sorted_l6 <- features_df
# RF_dse_classify_imp_sorted_l6 <- RF_dse_classify_imp_sorted_l6[, -1]
RF_dse_classify_imp_sorted_l6[,clean_otus_g:=fie_g_l6_62[RF_dse_classify_imp_sorted_l6, features_g, on = "features"]]
barplot(RF_dse_classify_imp_sorted_l6[1:62,MeanDecreaseAccuracy], names.arg=RF_dse_classify_imp_sorted_l6[1:62,"clean_otus_g"] , ylab="Mean Decrease in Accuracy (Variable Importance)", las=1, ylim=c(0,0.02), main="Classification disease")
as.data.table(RF_dse_classify_imp_sorted_l6)->RF_dse_classify_imp_sorted_l6
RF_dse_classify_imp_sorted_l6[1:62,]->Rfout_l6_62
Rfout_l6_62$new<-apply(Rfout_l6_62[,1:2],1,function(x){grep(max(x),Rfout_l6_62)})
Rfout_l6_62$abund<-apply(Rfout_l6_62[,7],1, function(x){colnames(Rfout_l6_62)[x]})
as.factor(Rfout_l6_62$abund)
Rfout_l6_62$abund<-factor(Rfout_l6_62$abund, labels=c("Active Crohn's", "Healthy Control"))
#Rfout_l6_62 <- na.omit(Rfout_l6_62)

## ---- chunk-13 ----
png("R objects/level-6_rf.png", width =800, height = 800)
par(mar=c(5,18,1,1))
plot(Rfout_l6_62$MeanDecreaseAccuracy, seq_along(Rfout_l6_62$clean_otus_g), yaxt="n", xlab="Decrease in Accuracy/Variable Importance", col=Rfout_l6_62$abund, pch=19, ylab="")
axis(side=2, at=1:62, labels=Rfout_l6_62$clean_otus_g, las=1, cex.axis=0.8)
abline(h=seq(1,62,1), col="gray64", lty=2)
legend(0.015,55.5, levels(Rfout_l6_62$abund), col=1:3, pch=19,cex = 0.9)
title(ylab = "Features at Genus level (selected 62 features from the RF_model)", line = 13, )  # Adjust the line value as needed

dev.off()



## ---- chunk-11 ----
heatmap(as.matrix(Rfout_l6_62[,.(active_cd,healthy_control)]),Rowv=NA,Colv=NA,cexCol=1,labRow=(Rfout_l6_62$clean_otus_g), margins=c(8,15), labCol=c("Active_CD","Healthy_Control"),scale="column",col = colorRampPalette(brewer.pal(8, "Reds"))(60), xlab="Disease_category")
legend(x = "topleft", legend = seq(-0.006,0.0122,0.003), cex = 0.6, fill = colorRampPalette(brewer.pal(8, "Reds"))(7), title=paste("Variable", "\n", " Importance"), bty="n")
## ----
# Assuming you have the 'Rfout_l6_62' data frame



hmr<- data.frame(
  features = Rfout_l6_62$clean_otus_g,
  active_cd = Rfout_l6_62$active_cd,
  healthy_control = Rfout_l6_62$healthy_control
)


data_long <- hmr %>%
  pivot_longer(cols = c(active_cd, healthy_control), names_to = "Category", values_to = "Value")

# Define the breaks and colors for your custom color scale
color_breaks <- c(0.00, 0.005, 0.010, 0.020)
custom_colors <- c("#8bd3c7", "#b2e061", "#fdcce5", "#fd7f6f")  # Replace with your desired colors

# Define the ggplot2 heatmap plot with the custom color scale
heatmap_plot_rf <- ggplot(data_long, aes(x = Category, y = reorder(features, -Value), fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(breaks = color_breaks, colors = custom_colors, guide = "legend") +
  labs(title = "Active_CD vs Healthy_control - RF Features", x = "Category", y = "Features") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the heatmap
png("./R objects/actvsheastool_RF_heatmap.png", width = 800, height = 600)
print(heatmap_plot_rf)
dev.off()

library(patchwork)
p <- heatmap_plot_boruta + heatmap_plot_rf
p


## extracting the top 62 features from rf to plot the heatmap for relative abundances
# Plot the features as heatmap
genera.clr <- genus.clr[, colnames(genus.clr) %in% top_129_marker_names]
# making ids for dys class
dys.ids_rf <- as.character(rownames(metar[metar$dse_status =='ACTIVE',]))
healthylike.ids_rf <- as.character(rownames(metar[metar$dse_status =='HEALTHY',]))
# merging dys class ids
dys_healthy_order_rf <- c(dys.ids_rf, healthylike.ids_rf)

genera.clr_rf <- genera.clr[dys_healthy_order_rf, ]
fullmap_rf <- metar[dys_healthy_order_rf, ]

hm.meta_rf <- data.frame(fullmap_rf[,'dse_status'], row.names = rownames(fullmap_rf), col=fullmap_rf$dse_status=='HEALTHY')
# color coding to dys class
hm.meta_rf$col[hm.meta_rf$col==TRUE] <- '#E69F00'
hm.meta_rf$col[hm.meta_rf$col==FALSE] <- '#009E73'

# modifying the featutres ids
features_ids_rf<- colnames(genera.clr_rf)
# Extract the last word based on the presence of "g__" 
last_words_boruta_rf <- sapply(features_ids_rf, function(str) {
  if (grepl("g__", str)) {
    gsub(".*;g__(.)", "\\1", str)
  } else if (grepl("f__", str)) {
    gsub(".+\\b(f__[A-Za-z]+)\\b", "\\1", str)
  } else {
    "No match"  # Handle cases where neither "g__" nor "f__" is found
  }
})
# adding the features ids to the colnames
colnames(genera.clr_rf) <- last_words_boruta_rf
# making new data.farme
sig_rf<- genera.clr_rf
sig_rf<- data.frame(dse_status= hm.meta_rf$fullmap_rf....dse_status.., genera.clr_rf, sampleid = rownames(sig_rf)) 
# # converting into long data frame for heatmap
# long_data_boruta_rf <- sig_rf%>%
#   pivot_longer(cols = -dys_class, 
#                names_to = "feature_name", 
#                values_to = "feature_value")
# 
## ---- chunk-4 ---
color_breaks <- c(0.00, 2.5, 7.5, 10.5)
custom_colors <- c("#8bd3c7", "#b2e061", "#fdcce5", "#fd7f6f")  # Replace with your desired colors

## ggplot heatmap

## data conversing uisng melt 
long_rf_dse<- melt(sig_rf)

# Convert sampleid to a factor with the same levels as in the original dataframe
long_rf_dse$sampleid <- factor(long_rf_dse$sampleid, levels = unique(sig_rf$sampleid))


## segment
# Create the heatmap plot with selected sample IDs as lines
heatmap_dse <- ggplot(long_rf_dse, aes(x = sampleid, y = reorder(variable, -value), fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(breaks = color_breaks, colors = custom_colors, guide = "legend") +
  labs(title = "ACTIVE vs. HEALTHY - RF Features", x = "", y = "Features") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis labels
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
  ) +
  geom_segment(
    aes(x = sampleid, xend = sampleid, y = -1, yend = -0.5, color = dse_status),
    size = 1
  ) +
  scale_color_manual(values = c("ACTIVE" = "#D5695D", "HEALTHY" = "#088158"))
print(heatmap_dse)

