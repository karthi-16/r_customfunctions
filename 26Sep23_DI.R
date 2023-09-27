# code obtained from saffouri et,al
setwd("/home/minionpak2/Desktop/my_edit/DSTCD_16Sanalysis_filtered/ACTvsHEL/STOOL/")
# packages
library(data.table)
library(readxl)
library(vegan)
library(plotrix)
library(tibble)
library(ggplot2)
library(gridExtra)
library(randomForest)
library(gplots)
library(heatmap3)
library(stringr)
library(Boruta)
library(phyloseq)
library(dplyr)
library(scales)
library(tidyr)
library(dysbiosisR)
library(caret)
load("./26Sep23_DI.RData")
source("./exportrmd_fucntion.R") # to export the chunks to rmd 
#  Dysbiosis index
# using dysbiosis R package.
physeq<- readRDS("/home/minionpak2/Desktop/my_edit/DSTCD_16Sanalysis_filtered/ACTvsHEL/STOOL/physeq.rds")
physeq
# check first five sample sums
sample_sums(physeq)[1:5]
table(sample_data(physeq)$dse_status)
# Dysbiosis Measures
dist.mat <- phyloseq::distance(physeq, "bray")
# get reference samples
ref.samples <- sample_names(subset_samples(physeq, 
                                           dse_status == "HEALTHY"))
## Cloud-based LOcally linear Unbiased Dysbiosis (CLOUD) test 

cloud.results <- cloudStatistic(physeq,
                                dist_mat = dist.mat,
                                reference_samples = ref.samples,
                                ndim=-1,
                                k_num=5)

# Here, the log2Stats can be used as a dysbiosis score.

# order the data 
cloud.results$disease <- factor(cloud.results$dse_status, 
                                levels = c("HEALTHY", "ACTIVE"))
# ROC
roc_4 <- pROC::roc(as.factor(cloud.results$dse_status),
                   cloud.results$log2Stats ,
                   #direction= "<",
                   plot=TRUE,
                   ci = TRUE,
                   auc.polygon=TRUE,
                   max.auc.polygon=TRUE,
                   print.auc=TRUE)

# plot the dysbiosis density plot
p2 <- plotDysbiosis(df=cloud.results,
                    xvar="dse_status",
                    yvar="log2Stats",
                    colors=c(ACTIVE="brown3", HEALTHY="steelblue"),
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis CLOUD Score")+
  theme_bw(base_size = 14)
p2

# define gradiant cols for gradiant plot
volcano <- c("#003f5c", "#58508d","#bc5090","#ff6361", "#ffa600")

# cutoff selection for log2stats based on the saffouri et,al paper and mainly the mean+2*sd in saffouri has changed to mean-2*SD in DST crohn's
mean(cloud.results[cloud.results$dse_status=="HEALTHY",]$log2Stats)-2*sd(cloud.results[cloud.results$dse_status=="HEALTHY",]$log2Stats) #-0.5924791

# Gradiant plot
## ---- chunk-1 ---
p.cloud <- plotDysbiosisGradient(df=cloud.results,
                                 score="log2Stats",
                                 high_line = -0.59,
                                 # low_line = normobiosis_thres,
                                 group_var = "dse_status",
                                 group_colors=c("HEALTHY" = "steelblue", 
                                                "ACTIVE"= "brown3"),
                                 point_size = 2,
                                 bg_colors = rev(volcano),
                                 jitter_width = 0.1) +
  labs(y="Log 2 stats", subtitle = "CLOUD test") +
  # adjust the x and y values to fit to plot
  ggplot2::annotate("text", x = 0.7, y = -0.55,
                    label = "Cut-off -0.59", color="white")
p.cloud
## ---

cutoff.val<- mean(cloud.results[cloud.results$dse_status=="HEALTHY",]$log2Stats)-2*sd(cloud.results[cloud.results$dse_status=="HEALTHY",]$log2Stats)

table(cloud.results$log2Stats > cutoff.val) # no.of true or false values greater than cutoff

# plot the log2stats aganist the disease group
plot(log2(cloud.results$stats) ~ cloud.results$dse_status, 
     xlab = "Disease Status", ylab = "log2(Stats)") 

# boxplot for dse_status
## ---- chunk-2 ---
boxplot<- ggplot(cloud.results, aes(x = dse_status, y = log2Stats, fill = dse_status)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(width = 0.2, height = 0) +  # Jitter points
  
  labs(x = "Disease Status", y = "log2(Stats)") +  # Labels
  theme_classic() +
  theme(axis.text.x = element_text(size = 15),  # Adjust X-axis label size
        axis.text.y = element_text(size = 15),  # Adjust Y-axis label size
        legend.text = element_text(size = 15))   # Adjust legend text size
boxplot
## ---

### CLOUD method avoids the potential bias by using a fresh approach. The log2stats and dys_class are the 
### scores and classifications used in the final analysis.

# Create a new column 'dys_class' based on the 'log2stats' column and the cutoff value
cloud.results$dys_class <- cloud.results$log2Stats 
cloud.results$dys_class <- ifelse(cloud.results$log2Stats < cutoff.val, 'T', 'F')
cloud.results$dys_class[cloud.results$dys_class == 'T'] <- 'dysbiotic'
cloud.results$dys_class[cloud.results$dys_class == 'F'] <- 'healthy-like'

#  Check out the associations for active and healthy
table(cloud.results$dse_status, cloud.results$dys_class)

# # plot log2stats against disease class
boxplot(log2(cloud.results$stats) ~ cloud.results$dys_class,
       xlab = "Disease Class", ylab = "log2(Stats)")

# boxplot for dys_class
## ---- chunk-3 ---
dysplot<- ggplot(cloud.results, aes(x = dys_class, y = log2Stats, fill = dys_class)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(width = 0.2, height = 0) +  # Jitter points
  
  labs(x = "Disease Class", y = "log2(Stats)") +  # Labels
  theme_minimal()+
  theme(axis.text.x = element_text(size = 15),  # Adjust X-axis label size
        axis.text.y = element_text(size = 15),  # Adjust Y-axis label size
        legend.text = element_text(size = 15)) 

dysplot
## ---


####### Run random forest model on the new classifications ######
####### and Boruta feature selection to yield significantly differentiating OTUs

# Read in the genus-level table, generated from the OTU table using QIIME2 
genus <- read.csv('/home/minionpak2/Desktop/my_edit/DSTCD_16Sanalysis_filtered/ACTvsHEL/STOOL/Analysis/L6/L6.tsv', 
                  header=1, check.names = F, row.names = 1, sep = '\t')

# CLR transform
genus.c <- t(genus); eps <- 0.5
genus.c <- genus.c * (1 - rowSums(genus.c==0) * eps / rowSums(genus.c))
genus.c[genus.c == 0] <- eps
genus.c <- sweep(genus.c, 1, rowSums(genus.c), '/')
ls <- log(genus.c)
genus.c <- t(ls - rowMeans(ls))
genus.c <- genus.c[, !is.nan(colSums(genus.c))]
genus.c <- as.data.frame(genus.c)

# rf.map <- tibble::column_to_rownames(cloud.results, var = 'Row.names')
rf.ids <- intersect(colnames(genus.c), rownames(cloud.results))
rd.ids <- sort(rf.ids)
genus.clr <- t(genus.c[, rf.ids])
dim(genus.clr)
# We have done CLR transform for OTU table, no need to use rare oTUs
# # R function for removing rare OTU's
# remove_rare <- function( table , cutoff_pro ) {
#   row2keep <- c()
#   cutoff <- ceiling( cutoff_pro * ncol(table) )  
#   for ( i in 1:nrow(table) ) {
#     row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
#     if ( row_nonzero > cutoff ) {
#       row2keep <- c( row2keep , i)
#     }
#   }
#   return( table [ row2keep , , drop=F ])
# }
# 
# otu_nonzero_counts_L6 <- apply(genus, 1, function(y) sum(length(which(y > 0))))
# hist(otu_nonzero_counts_L6, breaks=100, col="red", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")
# otu_table_rare_removed_L6 <- remove_rare(table=genus, cutoff_pro=0.05)
# dim(otu_table_rare_removed_L6)
# class(otu_table_rare_removed_L6)
# as.data.frame(otu_table_rare_removed_L6)->genus.clr
# genus.clr<- t(genus.clr)


# # Random forest model without cross validation for dys_class
# set.seed(12345)
# rfc.model <-randomForest(x = genus.clr, y=factor(cloud.results$dys_class), ntree = 1000, importance = T, keep.forest = T, mtry = 140)
# rfc.model   # OOB estimate of error rate 10.38%
# # dysbiotic healthy-like class.error
# # dysbiotic            5            7  0.58333333
# # healthy-like         4           90  0.04255319
# 
# rfc.model$importance
# rfc.prob <- predict(rfc.model, x=genus.clr, type = 'prob')
# cloud.results$cloud_prob <- rfc.prob[,1]
# 
# # extracting the importance from rf model
# rf_imp<- rfc.model$importance
# rownames(rf_imp) <- gsub(".*g__([^ ]+)$", "\\1", rownames(rf_imp))
# view(rf_imp) # shows the importance values for dysbiotic and healthy-like



## making a separate data.frame for genus with colnames modified
genus_df <- as.data.frame(genus.clr)
# Remove hyphens after 'g__' in column names
colnames(genus_df) <- gsub("g__", "g__temp", colnames(genus_df))
colnames(genus_df) <- sub("-", "", colnames(genus_df))
# Rename 'g__temp' back to 'g__'
colnames(genus_df) <- gsub("g__temp", "g__", colnames(genus_df))

# Boruta feature selection for the full genus features with dys_class of cloud.results
set.seed(123)
boruta.rfimp <- Boruta(x=genus_df, y=factor(cloud.results$dys_class))
length(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Confirmed']) # 10
sig.taxa <- lapply(X=names(boruta.rfimp$finalDecision[boruta.rfimp$finalDecision == 'Confirmed']),
                   FUN=function (xx) gsub(x=xx, pattern='.',replacement = ';', fixed = T))

sig.taxa
length(sig.taxa)


# Plot the features as heatmap
sig.genera.clr <- genus_df[, colnames(genus_df) %in% sig.taxa]
# making ids for dys class
dys.ids <- as.character(rownames(cloud.results[cloud.results$dys_class=='dysbiotic',]))
healthylike.ids <- as.character(rownames(cloud.results[cloud.results$dys_class=='healthy-like', ]))
# merging dys class ids
dys_healthy_order <- c(dys.ids, healthylike.ids)

sig.genera.clr <- sig.genera.clr[dys_healthy_order, ]
fullmap2 <- cloud.results[dys_healthy_order, ]

hm.meta <- data.frame(fullmap2[,'dys_class'], row.names = rownames(fullmap2), col=fullmap2$dys_class=='healthy-like')
# color coding to dys class
hm.meta$col[hm.meta$col==TRUE] <- '#E69F00'
hm.meta$col[hm.meta$col==FALSE] <- '#009E73'

# modifying the features ids
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
# making new data.farme
sig<- sig.genera.clr
sig<- data.frame(dys_class= hm.meta$fullmap2....dys_class.., sig.genera.clr, sampleid = rownames(sig)) 
# converting into long data frame for heatmap
# long_data_boruta <- sig%>%
#   pivot_longer(cols = -dys_class, 
#                names_to = "feature_name", 
#                values_to = "feature_value")

## data transform with melt using tidyr package 
long_data_boruta<- melt(sig)
# ggplot heatmap
## ---- chunk-4 ---
png("./R objects/actvshea_stool/heatmap_dysbiosis_fea_boruta.png")
heatmap_boruta_fea <- ggplot(long_data_boruta, aes(x = dys_class, y = reorder(feature_name, -feature_value), fill = feature_value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#8bd3c7", high = "#beb9db") +
  labs(title = "Dysbiotic vs. Healthy-Like - Boruta Features", x = "Category", y = "Features") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(heatmap_boruta_fea)
dev.off()
## ---

# # Print the heatmap
# png("./R objects/actvshea_stool/actvsheastool_boruta_dysbiosis_heatmap.png", width = 800, height = 600)
# print(heatmap_plot_boruta)
# dev.off()
# 
# # Create the heatmap with sample _ids
# png(paste0("./R objects/actvshea_stool/boruta_dysclass_persample.png"),  # create PNG for the heat map
#     width = 15*150,                        # 5 x 300 pixels
#     height = 6*150,
#     res = 300,                              # 300 pixels per inch
#     pointsize = 6)                          # smaller font size
# pastel_palette <- colorRampPalette(c("#8bd3c7", "#b2e061", "#7eb0d5", "#fd7f6f"))(150)
# # Generate the heatmap without the built-in legend
# heatmap3(
#   x = t(sig.genera.clr),
#   showColDendro = FALSE,
#   showRowDendro = FALSE,
#   ColSideColors = hm.meta$col,
#   cexRow = 1,
#   margins = c(2, 6),
#   Colv = NA,
#   Rowv = NULL,
#   col = pastel_palette

###running 
# # Create the heatmap without the built-in legend
# pastel_palette <- colorRampPalette(c("#b2e061", "#8bd3c7", "#7eb0d5", "#fd7f6f"))(7)
# png("./R objects/actvshea_stool/heatmap_ps_boruta.png", width = 1000, height = 600)
# heatmap_boruta<- heatmap3(
#   x = t(sig.genera.clr),
#   showColDendro = FALSE,
#   showRowDendro = FALSE,
#   ColSideColors = hm.meta$col,
#   ColSideLabs = "",
#   RowSideLabs = NA,
#   cexRow = 1.5,
#   margins = c(2, 2),
#   Colv = NA,
#   Rowv = NA,
#   col = pastel_palette,
#   labCol = NA,
#   cexCol = 1
#   )
# # Create a custom legend
# darkgreen <- "#006400"
# legend(
#   x = 0.1,
#   y = 0.9,
#   legend = c("Dysbiotic", "Healthy-Like"),  # Legend labels
#   fill = c(" darkgreen", "orange"),  # Corresponding colors
#   title = "Disease Class",  # Optional title for the legend
#   cex = 0.9 # Adjust the size of the legend text
# )
# dev.off()

                   
## ggplot heatmap with sampelids

long_data_boruta$sampleid <- factor(long_data_boruta$sampleid, levels = unique(sig$sampleid))
color_breaks <- c(0.00, 2.5, 7.5, 10.5)
custom_colors <- c("#8bd3c7", "#b2e061", "#fdcce5", "#fd7f6f")  # Replace with your desired colors
## segment
# Create the heatmap plot with selected sample IDs as lines
## ---- chunk-
heatmap <- ggplot(long_data_boruta, aes(x = sampleid, y = reorder(variable, -value), fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(breaks = color_breaks, colors = custom_colors, guide = "legend") +
  labs(title = "Dysbiotic vs. Healthy-Like - Boruta Features", x = "", y = "Features") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis labels
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
  ) +
  geom_segment(
    aes(x = sampleid, xend = sampleid, y = -1, yend = -0.5, color = dys_class),
    size = 1
  ) +
  scale_color_manual(values = c("dysbiotic" = "#D5695D", "healthy-like" = "#088158"))
print(heatmap)
## ---


## Random forest ##
set.seed(123)
# Random forest analysis with 10 fold cross validation ---- codes from rf.script.R
RF_models <- vector("list", 5)
# Loop through 5 runs of cross-validation
for (i in 1:5) {
  RF_models[[i]] <- train(
    x = genus_df, 
    y = cloud.results$dys_class,
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

## Random forest with 10 fold 
# Calculate the average error rate for each mtry value across runs
avg_error_data <- aggregate(error ~ mtry, data = error_data, FUN = mean)
as.data.table(avg_error_data)->avg_error_data
error_data[,avg:=avg_error_data[error_data,error,on="mtry"]]
dcast(error_data, mtry+avg ~ run,value.var="error")->kdat
apply(kdat[,c(3:7)], 1, mean)
# [1] 0.1521818 0.1186667 0.1166667 0.1224848 0.1257879 0.1279697 0.1243030 0.1246667 0.1244848 0.1224848 0.1243030 0.1148485 0.1168485 0.1186667 0.1186667 0.1186667
# [17] 0.1224848 0.1168485
apply(kdat[,c(3:7)], 1, sd)
# [1] 0.008500122 0.004884337 0.004853217 0.006447499 0.004931114 0.010514540 0.008052856 0.008208710 0.008272228 0.006688669 0.010857200 0.010672949 0.011104314
# [14] 0.009008719 0.013117397 0.008521566 0.009446092 0.010725732
min(apply(kdat[,c(3:7)], 1, mean))
#[1] 0.1148485
## Therefore minimum error + sd is 0.1255214, and optimum number of vars is 65

# Plotting the errors for all 5 runs with average line
## ---- chunk-5 ---
ggplot() +
  geom_line(data = error_data, aes(x = mtry, y = error, group = run, col = "gray"), alpha = 0.3) +  # Individual error lines
  geom_line(data = avg_error_data, aes(x = mtry, y = error), color = "black", linewidth = 0.5) +  # Average line
  geom_point(data = error_data, aes(x = mtry, y = error, col= "gray")) +  # Individual points
  labs(x = "Number of Variables", y = "Error Rate",
       title = "Error Rate for Different mtry Values in Random Forest (5 Runs)") +
  #scale_color_manual(values = c("blue", "green", "purple", "orange", "pink")) +
  theme_minimal()+
  geom_vline(xintercept=65,col="green")+
  guides( colour= FALSE)
## ---

# Extract variable importance scores from each model
importance_scores <- lapply(RF_models, function(model) varImp(model$finalModel))
# Extract variable importance scores for the first model
importance_scores_first_model <- varImp(RF_models[[1]]$finalModel)

# Extract variable importance scores for the 34 variables
# extract the rownames
rownames<- rownames(importance_scores_first_model)
# Convert importance_scores to a data frame
importance_df <- as.data.table(importance_scores_first_model)
importance_df$row_ids <- rownames

# Sort importance scores in descending order and select the top 65 markers
top_65_markers <- head(importance_df[order(-importance_df$Overall), ], 65)
top_65_markers_names<- top_65_markers$row_ids


# # Create a subset of the dataframe using column IDs in 'top'
# #subset_df <- genus.clr[, top_34_marker_names]
# 
# train_34 <- genus.clr[, top_34_marker_names]
# colnames(train_34)
# train_34<- as.data.frame(train_34)
# train_34$dys_class<- cloud.results$dys_class
# 
# # rf_model for original top 34 variables 
# set.seed(123)
# rf_34 <- randomForest(
#   x = genus_df,
#   y = factor(cloud.results$dys_class),  
#   importance = TRUE,
#   ntree = 500,
#  mtry = 65
# )
# rf_34
# 
# # Importance values from random forest output
# rf_34_df<- as.data.frame(rf_34$importance)
# # new data.frame for importance
# rfh<- data.frame(
#   features = rownames(rf_34_df),
#   dysbiotic = rf_34_df$dysbiotic,
#   healthy_like = rf_34_df$`healthy-like`)
# # order the values for dysbiotic
# rfhm<- rfh[order(-rfh$dysbiotic), ]
# # copy the features lables
# fea_ids<- rfhm$features
# # Extract the last word based on the presence of "g__" 
# last_words <- sapply(fea_ids, function(str) {
#   if (grepl("g__", str)) {
#     gsub(".*;g__(.)", "\\1", str)
#   } else if (grepl("f__", str)) {
#     gsub(".+\\b(f__[A-Za-z]+)\\b", "\\1", str)
#   } else {
#     "No match"  # Handle cases where neither "g__" nor "f__" is found
#   }
# })
# 
# # import the modified labels to the features.
# rfhm$features<- last_words
# # Reshape the data from wide to long format using tidyr (assuming you have more columns)
# data_long <- rfhm %>%
#   pivot_longer(cols = c(dysbiotic, healthy_like), names_to = "Category", values_to = "Value")
## older codes not used for ggplot heatmap ##
# data_long <- data_long[order(data_long$Value),]
# # Reorder the 'Category' variable based on the 'Value' variable
# # Get the unique values of 'Category' in the desired order
# category_order <- unique(data_long$Category)
# 
# # Set 'Category' as an ordered factor with the same order as 'category_order'
# data_long$Category <- factor(data_long$Category, levels = category_order, ordered = TRUE)

# # Create the heatmap using ggplot2 with 'features' in the same order as in the dataset
# heatmap_plot <- ggplot(data_long, aes(x = Category, y = reorder(features, -Value), fill = Value)) +
#   geom_tile(color = "white") +
#   scale_fill_gradient2(low = "#8bd3c7",mid = "#fdcce5",  high = "#bd7ebe") +
#   labs(title = "Dysbiotic vs. Healthy-Like Heatmap", x = "Category", y = "features") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# # Define the breaks and colors for your custom color scale
# color_breaks <- c(0.00, 0.025, 0.050, 0.070, 0.09)
# custom_colors <- c("#8bd3c7", "#b2e061", "#7eb0d5", "#fdcce5", "#fd7f6f")  # Replace with your desired colors
# 
# # Define the ggplot2 heatmap plot with the custom color scale
# ## ---- cunnk-6 ---
# heatmap_plot_rf <- ggplot(data_long, aes(x = Category, y = reorder(features, -Value), fill = Value)) +
#   geom_tile(color = "white") +
#   scale_fill_gradientn(breaks = color_breaks, colors = custom_colors, guide = "legend") +
#   labs(title = "Dysbiotic vs. Healthy-Like - RF Features", x = "Category", y = "Features") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# print(heatmap_plot_rf)
# ## ---
# 
# # Print the heatmap
# png("./R objects/actvsheastool_RFdysbiosis_heatmap.png", width = 800, height = 600)
# print(heatmap_plot_rf)
# dev.off()
# 
# library(patchwork)
# dy <- heatmap_plot_boruta + heatmap_plot_rf
# ## ---- chunk-7 ---
# dy
# ## ---

### heatmap with sampleids

## extracting the top 65 features from rf to plot the heatmap for relative abundances
# Plot the features as heatmap
genera.clr <- genus_df[, colnames(genus_df) %in% top_65_markers_names]
# making ids for dys class
dys.ids_rf <- as.character(rownames(cloud.results[cloud.results$dys_class=='dysbiotic',]))
healthylike.ids_rf <- as.character(rownames(cloud.results[cloud.results$dys_class=='healthy-like', ]))
# merging dys class ids
dys_healthy_order_rf <- c(dys.ids_rf, healthylike.ids_rf)

genera.clr_rf <- genera.clr[dys_healthy_order_rf, ]
fullmap_rf <- cloud.results[dys_healthy_order_rf, ]

hm.meta_rf <- data.frame(fullmap_rf[,'dys_class'], row.names = rownames(fullmap_rf), col=fullmap_rf$dys_class=='healthy-like')
# color coding to dys class
hm.meta_rf$col[hm.meta_rf$col==TRUE] <- '#E69F00'
hm.meta_rf$col[hm.meta_rf$col==FALSE] <- '#009E73'

# modifying the features ids
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
sig_rf<- data.frame(dys_class= hm.meta_rf$fullmap_rf....dys_class.., genera.clr_rf, sampleid = rownames(sig_rf)) 
# # converting into long data frame for heatmap
# long_data_boruta_rf <- sig_rf%>%
#   pivot_longer(cols = -dys_class, 
#                names_to = "feature_name", 
#                values_to = "feature_value")
# 

## ggplot heatmap

## data conversing uisng melt 
long_rf<- melt(sig_rf)

# Convert sampleid to a factor with the same levels as in the original dataframe
long_rf$sampleid <- factor(long_rf$sampleid, levels = unique(sig_rf$sampleid))


## segment
# Create the heatmap plot with selected sample IDs as lines
## ---- chunk-4 ---
color_breaks <- c(0.00, 2.5, 7.5, 10.5)
custom_colors <- c("#8bd3c7", "#b2e061", "#fdcce5", "#fd7f6f")  # Replace with your desired colors
heatmap_dys <- ggplot(long_rf, aes(x = sampleid, y = reorder(variable, -value), fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(breaks = color_breaks, colors = custom_colors, guide = "legend") +
  labs(title = "Dysbiotic vs. Healthy-Like - RF Features", x = "", y = "Features") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis labels
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
  ) +
  geom_segment(
    aes(x = sampleid, xend = sampleid, y = -1, yend = -0.5, color = dys_class),
    size = 1
  ) +
  scale_color_manual(values = c("dysbiotic" = "#D5695D", "healthy-like" = "#088158"))
print(heatmap_dys)






# png("./R objects/actvshea_stool/heatmap_dysbiosis_fea_rf.png", width = 800, height = 600)
# heatmap_dysbiosis_fea_rf <- ggplot(long_data_boruta_rf, aes(x = reorder(dys_class, -feature_value), y = reorder(feature_name, -feature_value), fill = feature_value)) +
#   geom_tile(color = "white") +
#   scale_fill_gradientn(breaks = color_breaks, colors = custom_colors, guide = "legend") +
#   labs(title = "Dysbiotic vs. Healthy-Like - RF Features", x = "Category", y = "Features") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# print(heatmap_dysbiosis_fea_rf)
# ## ---
# dev.off()
# ## per sample heatmap
# pastel_palette <- colorRampPalette(c("#b2e061", "#8bd3c7", "#7eb0d5", "#fd7f6f"))(7)
# png("./R objects/actvshea_stool/heatmap_dysbiosis_ps_rf.png", width = 1000, height = 600)
# ## ---- chunk--
# heatmap3(
#   x = t(genera.clr_rf[,1:64]),
#   showColDendro = FALSE,
#   showRowDendro = FALSE,
#   ColSideColors = hm.meta_rf$col,
#   ColSideLabs = "",
#   RowSideLabs = NA,
#   cexRow = 0.8,
#   margins = c(2, 2),
#   Colv = NA,
#   Rowv = NA,
#   col = pastel_palette,
#   labCol = NA,
#   cexCol = 1
# )
# # Create a custom legend
# darkgreen <- "#006400"
# legend(
#   x = 0.1,
#   y = 0.9,
#   legend = c("Dysbiotic", "Healthy-Like"),  # Legend labels
#   fill = c(" darkgreen", "orange"),  # Corresponding colors
#   title = "Disease Class",  # Optional title for the legend
#   cex = 0.9 # Adjust the size of the legend text
# )
# ##
# dev.off()
# 
