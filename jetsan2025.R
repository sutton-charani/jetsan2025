# Variables and screen cleaning:
graphics.off(); cat("\014"); rm(list=ls()); options(warn=-1)

library(data.table); library(ggplot2); library(plotly); library(reshape2); library(FactoMineR); library(entropy)
library(e1071); library(moments); library(nonlinearTseries); library(class); library(rpart); library(rpart.plot)
library(randomForest); library(lightgbm)

###############################
mean_impute <- function(df, mean_values=NULL){
  completed_df <- df
  for (j in 1:ncol(completed_df)){
    v <- completed_df[, j]
    if (is.numeric(v)){
      if (sum(is.na(v)) > 0){
        if (is.null(mean_values)){
          v[is.na(v)] <- mean(v, na.rm=T)
          completed_df[, j] <- v
        } else {
          v[is.na(v)] <- mean_values[j]
          completed_df[, j] <- v
        }
      }
    }
  }
  return(completed_df)
}
###############################
compute_descriptors <- function(v, nquant=11, ndfa_sizes=10){
  if (length(unique(v[!is.na(v)])) > 1){
    if (length(v[!is.na(v)]) > 10){
      dfa_model <- dfa(v[!is.na(v)], window.size.range=c(10, round(length(v)/2)), 
                       npoints=ndfa_sizes, do.plot=F)
      alpha <- as.numeric(estimate(dfa_model))
    } else {
      alpha <- NA
    }
    descriptors_vec <- c(mean(v, na.rm=T), var(v, na.rm=T), sd(v, na.rm=T),
                         quantile(v, probs=seq(0, 1, 1/(nquant-1)), na.rm=T),
                         entropy(v[!is.na(v)]), kurtosis(v, na.rm=T), skewness(v, na.rm=T), 
                         alpha)
  } else {
    descriptors_vec <- rep(NA, 39)
  }
  
  names(descriptors_vec) <- c("mean", "var", "sd", paste0('q', seq(1, nquant)), 
                              "entropy", "kurtosis", "skewness", "alpha")
  return(descriptors_vec)
}
###############################

# Data import online
temp <- tempfile(fileext=".csv")
download.file("https://drive.google.com/uc?authuser=0&id=1RP_SSKpmHySDLcH-QAjP2xM2-KZ9SVm1", temp)
df <- read.csv(temp)
# or locally: df <- data.frame(fread('actimetry_jetsan2025_dataset', sep=';'))

df <- subset(df, select=-c(ENMO_rightwrist, ENMO_leftwrist, id_move))
names(df)[1:7] <- c('time', 'acc_right_x', 'acc_right_y', 'acc_right_z', 
                    'acc_left_x', 'acc_left_y', 'acc_left_z')
for(j in 1:7){
  df[, j] <- as.numeric(df[, j])
}
df$acc_eucl_left <- sqrt(df$acc_left_x^2 + df$acc_left_y^2 + df$acc_left_z^2)
df$acc_eucl_right <- sqrt(df$acc_right_x^2 + df$acc_right_y^2 + df$acc_right_z^2)

# PCA
cat('pca, ')
res.pca_test <- PCA(subset(df, select=-c(time, label, id_subject, id)), graph=F)
df$acc_pca <- res.pca_test$ind$coord[, 1]

df$id <- factor(df$id)

start <- Sys.time()
i_id=1
cat(length(levels(df$id)), 'activities to describe ')
descriptors <- c()
activity_labs <- c()
set.seed(1234)
for (i_id in 1:length(levels(df$id))){
  
  # Data loading
  cat(i_id, ' ')
  df_1activity <- df[df$id == unique(df$id)[i_id], ]
  
  # Descriptors computation
  descriptors_right_x_vec <- compute_descriptors(df_1activity$acc_right_x)
  descriptors_right_y_vec <- compute_descriptors(df_1activity$acc_right_y)
  descriptors_right_z_vec <- compute_descriptors(df_1activity$acc_right_z)
  descriptors_right_eucl_vec <- compute_descriptors(df_1activity$acc_eucl_right)
  descriptors_left_x_vec <- compute_descriptors(df_1activity$acc_left_x)
  descriptors_left_y_vec <- compute_descriptors(df_1activity$acc_left_y)
  descriptors_left_z_vec <- compute_descriptors(df_1activity$acc_left_z)
  descriptors_left_eucl_vec <- compute_descriptors(df_1activity$acc_eucl_left)
  descriptors_pca_vec <- compute_descriptors(df_1activity$acc_pca)
  names(descriptors_right_x_vec) <- paste0(names(descriptors_right_x_vec), '_right_x')
  names(descriptors_right_y_vec) <- paste0(names(descriptors_right_y_vec), '_right_y')
  names(descriptors_right_z_vec) <- paste0(names(descriptors_right_z_vec), '_right_z')
  names(descriptors_right_eucl_vec) <- paste0(names(descriptors_right_eucl_vec), '_right_eucl')
  names(descriptors_left_x_vec) <- paste0(names(descriptors_left_x_vec), '_left_x')
  names(descriptors_left_y_vec) <- paste0(names(descriptors_left_y_vec), '_left_y')
  names(descriptors_left_z_vec) <- paste0(names(descriptors_left_z_vec), '_left_z')
  names(descriptors_left_eucl_vec) <- paste0(names(descriptors_left_eucl_vec), '_left_eucl')
  names(descriptors_pca_vec) <- paste0(names(descriptors_pca_vec), '_pca')
  
  descriptors <- rbind(descriptors, 
                       c(descriptors_right_x_vec, descriptors_right_y_vec, descriptors_right_z_vec, 
                         descriptors_right_eucl_vec, descriptors_left_x_vec, descriptors_left_y_vec, 
                         descriptors_left_z_vec, descriptors_left_eucl_vec, descriptors_pca_vec))
  activity_labs <- c(activity_labs, df_1activity$label[1])
}
end <- Sys.time()
end-start

descriptors_df <- data.frame(descriptors, label=activity_labs)

# write.csv(descriptors, 'descriptors.csv', row.names=F)
# descriptors_df  <- data.frame(fread('descriptors.csv'))

# Correlation heatmap
descriptors_right_eucl <- descriptors_df[, grep("right_eucl", names(descriptors_df))]
descriptors_left_eucl <- descriptors_df[, grep("left_eucl", names(descriptors_df))]
cormat_right <- cor(descriptors_right_eucl, method="spearman", use="pairwise.complete.obs")
cormat_left <- cor(descriptors_left_eucl, method="spearman", use="pairwise.complete.obs")
for (j in 1:nrow(cormat_right)){
  colnames(cormat_right)[j] <- strsplit(colnames(cormat_right), '_')[[j]][1]
  rownames(cormat_right)[j] <- strsplit(rownames(cormat_right), '_')[[j]][1]
  colnames(cormat_left)[j] <- strsplit(colnames(cormat_left), '_')[[j]][1]
  rownames(cormat_left)[j] <- strsplit(rownames(cormat_left), '_')[[j]][1]
}
cormat <- round((cormat_right+cormat_left)/2, 2)
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm=TRUE)
ggplot(data=melted_cormat, aes(Var2, Var1, fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, 
                       limit=c(-1,1), space="Lab", name="Spearman\nCorrelation") +
  theme_minimal() + 
  theme() +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label=value), color="black", size=3) +
  theme(
    axis.text.x=element_text(angle=45, vjust=1, size=25, hjust=1),
    axis.text.y=element_text(size=25),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid.major=element_blank(),
    panel.border=element_blank(),
    panel.background=element_blank(),
    axis.ticks=element_blank(),
    legend.justification=c(1, 0),
    legend.position=c(0.5, 0.8),
    legend.direction="horizontal",
    legend.text=element_text(size=15),
    legend.title=element_text(size=20)
  ) +
  guides(fill=guide_colorbar(barwidth=15, barheight=1.5, title.position="top", title.hjust=0.5))

# PCA
cat('pca, ')
res.pca <- PCA(subset(descriptors_df, select=-label), graph=F)
fviz_pca_ind(res.pca, geom="point", habillage=factor(descriptors_df$label), 
             alpha.ind=0.8, pointsize=7, invisible=c("ind.sup", "quali", "var", "quanti.sup"),
             addEllipses=T) + 
  scale_shape_manual(values=rep(20, 10)) +
  scale_color_paletteer_d("pals::glasbey") +
  theme(
    legend.text = element_text(size=13),
    legend.title = element_text(size=0),
    axis.title = element_text(size=18),
    plot.title = element_blank()
  )

# Models comparison pipeline (leave-1-out activity cross valisation)
cat('leave-1-out activity pipeline:', nrow(descriptors_df), 'activities to process')
descriptors_df$label <- factor(descriptors_df$label)
start <- Sys.time()
iactivity <- 1
set.seed(1234)
dummy_accs <- c()
knn_accs <- c(); nb_accs <- c(); tree_accs <- c(); svm_accs <- c(); forest_accs <- c(); lgbm_accs <- c()
knn_accs_single <- c(); nb_accs_single <- c(); tree_accs_single <- c(); svm_accs_single <- c(); forest_accs_single <- c(); lgbm_accs_single <- c()
knn_accs_eucl <- c(); nb_accs_eucl <- c(); tree_accs_eucl <- c(); svm_accs_eucl <- c(); forest_accs_eucl <- c(); lgbm_accs_eucl <- c()
knn_accs_pc1 <- c(); nb_accs_pc1 <- c(); tree_accs_pc1 <- c(); svm_accs_pc1 <- c(); forest_accs_pc1 <- c(); lgbm_accs_pc1 <- c()
for (iactivity in 1:nrow(descriptors_df)){
  cat(iactivity, ' ')
  
  # Train-test data shaping
  df_train <- descriptors_df[-iactivity, ]
  df_test <- descriptors_df[iactivity, ]
  
  # Filtering on descriptors computed on single dimensions (x, y, z), on euclidean norm or on pc1
  df_train_single <- df_train[, c(names(df_train)[c(grep('_x', names(df_train)), 
                                                    grep('_y', names(df_train)), 
                                                    grep('_z', names(df_train)))], 'label')]
  df_test_single <- df_test[, c(names(df_test)[c(grep('_x', names(df_test)), 
                                                 grep('_y', names(df_test)), 
                                                 grep('_z', names(df_test)))], 'label')]
  df_train_eucl <- df_train[, c(names(df_train)[grep('_eucl', names(df_train))], 'label')]
  df_test_eucl <- df_test[, c(names(df_test)[grep('_eucl', names(df_test))], 'label')]
  df_train_pc1 <- df_train[, c(names(df_train)[grep('_pca', names(df_train))], 'label')]
  df_test_pc1 <- df_test[, c(names(df_test)[grep('_pca', names(df_test))], 'label')]
  
  # PCA
  res.pca_train <- PCA(subset(df_train, select=-c(label)), ncp=10, graph=F) 
  res.pca_train_single <- PCA(subset(df_train_single, select=-c(label)), ncp=10, graph=F) 
  res.pca_train_eucl <- PCA(subset(df_train_eucl, select=-c(label)), ncp=10, graph=F) 
  res.pca_train_pc1 <- PCA(subset(df_train_pc1, select=-c(label)), ncp=10, graph=F) 
  df_train_pca <- data.frame(res.pca_train$ind$coord, label=df_train$label)
  df_train_single_pca <- data.frame(res.pca_train_single$ind$coord, label=df_train$label)
  df_train_eucl_pca <- data.frame(res.pca_train_eucl$ind$coord, label=df_train$label)
  df_train_pc1_pca <- data.frame(res.pca_train_pc1$ind$coord, label=df_train$label)
  df_test_pca <- data.frame(predict(res.pca_train, df_test)$coord, label=df_test$label)
  df_test_single_pca <- data.frame(predict(res.pca_train_single, df_test_single)$coord, label=df_test$label)
  df_test_eucl_pca <- data.frame(predict(res.pca_train_eucl, df_test_eucl)$coord, label=df_test$label)
  df_test_pc1_pca <- data.frame(predict(res.pca_train_pc1, df_test_pc1)$coord, label=df_test$label)
  
  # Training
  naive_bayes_model <- naiveBayes(label ~ ., df_train_pca)
  tree <- rpart(label ~ ., df_train)
  svm_model <- svm(label ~ ., df_train_pca)
  forest <- randomForest(label ~ ., df_train)
  lgbm_model <- lightgbm(data=as.matrix(subset(df_train, select=-label)), 
                         label=df_train$label, verbose=-1)
  naive_bayes_model_single <- naiveBayes(label ~ ., df_train_single_pca)
  tree_single <- rpart(label ~ ., df_train_single)
  svm_model_single <- svm(label ~ ., df_train_single_pca)
  forest_single <- randomForest(label ~ ., df_train_single)
  lgbm_model_single <- lightgbm(data=as.matrix(subset(df_train_single, select=-label)), 
                                label=df_train_single$label, verbose=-1)
  naive_bayes_model_eucl <- naiveBayes(label ~ ., df_train_eucl_pca)
  tree_eucl <- rpart(label ~ ., df_train_eucl)
  svm_model_eucl <- svm(label ~ ., df_train_eucl_pca)
  forest_eucl <- randomForest(label ~ ., df_train_eucl)
  lgbm_model_eucl <- lightgbm(data=as.matrix(subset(df_train_eucl, select=-label)), 
                              label=df_train_eucl$label, verbose=-1)
  naive_bayes_model_pc1 <- naiveBayes(label ~ ., df_train_pc1_pca)
  tree_pc1 <- rpart(label ~ ., df_train_pc1)
  svm_model_pc1 <- svm(label ~ ., df_train_pc1_pca)
  forest_pc1 <- randomForest(label ~ ., df_train_pc1)
  lgbm_model_pc1 <- lightgbm(data=as.matrix(subset(df_train_pc1, select=-label)), 
                             label=df_train_pc1$label, verbose=-1)
  
  # Predictions
  dummy_preds <- rep("climb_the_stairs", nrow(df_test))
  knn_preds <- as.character(knn(subset(df_train_pca, select=-label), 
                                subset(df_test_pca, select=-label), 
                                df_train_pca$label))
  nb_preds <- as.character(predict(naive_bayes_model, df_test_pca))
  tree_preds <- as.character(predict(tree, df_test, type="class"))
  svm_preds <- as.character(predict(svm_model, df_test_pca))
  forest_preds <- as.character(predict(forest, df_test))
  lgbm_preds <- as.character(predict(lgbm_model, as.matrix(subset(df_test, select=-label)), type="class"))
  
  knn_single_preds <-as.character(knn(subset(df_train_single_pca, select=-label), 
                                      subset(df_test_single_pca, select=-label), 
                                      df_train_single_pca$label))
  nb_single_preds <- as.character(predict(naive_bayes_model_single, df_test_single_pca))
  tree_single_preds <- as.character(predict(tree_single, df_test_single, type="class"))
  svm_single_preds <- as.character(predict(svm_model_single, df_test_single_pca))
  forest_single_preds <- as.character(predict(forest_single, df_test_single))
  lgbm_single_preds <- as.character(predict(lgbm_model_single, 
                                            as.matrix(subset(df_test_single, select=-label)), type="class"))
  
  knn_eucl_preds <- as.character(knn(subset(df_train_eucl_pca, select=-label), 
                                     subset(df_test_eucl_pca, select=-label), 
                                     df_train_eucl_pca$label))
  nb_eucl_preds <- as.character(predict(naive_bayes_model_eucl, df_test_eucl_pca))
  tree_eucl_preds <- as.character(predict(tree_eucl, df_test_eucl, type="class"))
  svm_eucl_preds <- as.character(predict(svm_model_eucl, df_test_eucl_pca))
  forest_eucl_preds <- as.character(predict(forest_eucl, df_test_eucl))
  lgbm_eucl_preds <- as.character(predict(lgbm_model_eucl, 
                                          as.matrix(subset(df_test_eucl, select=-label)), type="class"))
  
  knn_preds_pc1 <- as.character(knn(subset(df_train_pc1_pca, select=-label), 
                                    subset(df_test_pca, select=-label), 
                                    df_train_pc1_pca$label))
  nb_preds_pc1 <- as.character(predict(naive_bayes_model_pc1, df_test_pc1_pca))
  tree_preds_pc1 <- as.character(predict(tree_pc1, df_test_pc1, type="class"))
  svm_preds_pc1 <- as.character(predict(svm_model_pc1, df_test_pc1_pca))
  forest_preds_pc1 <- as.character(predict(forest_pc1, df_test_pc1))
  lgbm_preds_pc1 <- as.character(predict(lgbm_model_pc1, 
                                         as.matrix(subset(df_test_pc1, select=-label)), type="class"))
  
  # Evaluation
  dummy_accs <- c(dummy_accs, mean(dummy_preds == df_test$label))
  knn_accs <- c(knn_accs, mean(knn_preds == df_test$label))
  nb_accs <- c(nb_accs, mean(nb_preds == df_test$label))
  tree_accs <- c(tree_accs, mean(tree_preds == df_test$label))
  svm_accs <- c(svm_accs, mean(svm_preds == df_test$label))
  forest_accs <- c(forest_accs, mean(forest_preds == df_test$label))
  lgbm_accs <- c(lgbm_accs, mean(lgbm_preds == df_test$label))
  
  knn_accs_single <- c(knn_accs_single, mean(knn_single_preds == df_test$label))
  nb_accs_single <- c(nb_accs_single, mean(nb_single_preds == df_test$label))
  tree_accs_single <- c(tree_accs_single, mean(tree_single_preds == df_test$label))
  svm_accs_single <- c(svm_accs_single, mean(svm_single_preds == df_test$label))
  forest_accs_single <- c(forest_accs_single, mean(forest_single_preds == df_test$label))
  lgbm_accs_single <- c(lgbm_accs_single, mean(lgbm_single_preds == df_test$label))
  
  knn_accs_eucl <- c(knn_accs_eucl, mean(knn_eucl_preds == df_test$label))
  nb_accs_eucl <- c(nb_accs_eucl, mean(nb_eucl_preds == df_test$label))
  tree_accs_eucl <- c(tree_accs_eucl, mean(tree_eucl_preds == df_test$label))
  svm_accs_eucl <- c(svm_accs_eucl, mean(svm_eucl_preds == df_test$label))
  forest_accs_eucl <- c(forest_accs_eucl, mean(forest_eucl_preds == df_test$label))
  lgbm_accs_eucl <- c(lgbm_accs_eucl, mean(lgbm_eucl_preds == df_test$label))
  
  knn_accs_pc1 <- c(knn_accs_pc1, mean(knn_preds_pc1 == df_test$label))
  nb_accs_pc1 <- c(nb_accs_pc1, mean(nb_preds_pc1 == df_test$label))
  tree_accs_pc1 <- c(tree_accs_pc1, mean(tree_preds_pc1 == df_test$label))
  svm_accs_pc1 <- c(svm_accs_pc1, mean(svm_preds_pc1 == df_test$label))
  forest_accs_pc1 <- c(forest_accs_pc1, mean(forest_preds_pc1 == df_test$label))
  lgbm_accs_pc1 <- c(lgbm_accs_pc1, mean(lgbm_preds_pc1 == df_test$label))
}
end <- Sys.time()
end-start
results <- data.frame(dummy=dummy_accs, 
                      knn=knn_accs, knn_single=knn_accs_single,  knn_eucl=knn_accs_eucl, knn_pc1=knn_accs_pc1,
                      naive_bayes=nb_accs, naive_bayes_single=nb_accs_single, naive_bayes_eucl=nb_accs_eucl, naive_bayes_pc1=nb_accs_pc1,
                      tree=tree_accs, tree_single=tree_accs_single, tree_eucl=tree_accs_eucl, tree_pc1=tree_accs_pc1,
                      svm=svm_accs,  svm_single=svm_accs_single, svm_eucl=svm_accs_eucl, svm_pc1=svm_accs_pc1,
                      forest=forest_accs, forest_single=forest_accs_single, forest_eucl=forest_accs_eucl, forest_pc1=forest_accs_pc1,
                      lightGBM=lgbm_accs, lightGBM_single=lgbm_accs_single, lightGBM_eucl=lgbm_accs_eucl, lightGBM_pc1=lgbm_accs_pc1)

#write.csv(results, './data/results_aliceDB.csv', row.names=F)

results <- data.frame(fread('./data/results_aliceDB.csv'))

round(colMeans(results), 2)

melted_results <- melt(colMeans(results))
melted_results
melted_results$model <- row.names(melted_results)
row.names(melted_results) <- NULL

# Columns renaming
v <- c()
for (i in 1:nrow(melted_results)){
  if (grepl("single", melted_results$model[i])){
    v <- c(v, "single")
    melted_results$model[i] <- substr(melted_results$model[i], 1, nchar(melted_results$model[i])-7)
  } else {
    if (grepl("eucl", melted_results$model[i])){
      v <- c(v, "eucl")
      melted_results$model[i] <- substr(melted_results$model[i], 1, nchar(melted_results$model[i])-5)
    } else {
      if (grepl("pc1", melted_results$model[i])){
        v <- c(v, "pc1")
        melted_results$model[i] <- substr(melted_results$model[i], 1, nchar(melted_results$model[i])-4)
      } else {
        v <- c(v, "all")
      }
    }
  }
}
melted_results$target_dimension <- v
melted_results$target_dimension <- factor(melted_results$target_dimension, levels=c("eucl", "pc1", "single", "all"))
melted_results$model <- factor(melted_results$model, levels=c("dummy", "knn", "naive_bayes", "tree", "svm", "forest", "lightGBM"))

ggplot(melted_results, aes(x=model)) +
  geom_bar(aes(y=value, fill=target_dimension), stat="identity", position="dodge") +
  geom_hline(yintercept=max(melted_results$value), lty="dashed") +
  xlab("") + ylab("accuracy") +
  theme_bw() +
  theme(
    axis.text.x=element_text(size=23),
    axis.title.y=element_text(size=27),
    axis.text.y = element_text(size=18),
    legend.position="top",
    legend.text = element_text(size=20),
    legend.title = element_text(size=23)
  ) + scale_fill_manual(values=wes_palette("GrandBudapest2", n = 4)) +
  scale_y_continuous(breaks=c(seq(0,75,0.25), round(max(melted_results$value), 2)), limits=c(0,1)) +
  labs(fill = "computed time series")

descriptors_df$label <- factor(descriptors_df$label)

n_forests <- 1000
set.seed(1234)
imps <- c()
for (i in 1:n_forests){
  cat(i, ' ')
  global_forest <- randomForest(label ~ ., descriptors_df, importance=T)
  imp <- rowMeans(global_forest$importance[, 11:12])
  imps <- rbind(imps, imp)
}
imp <- colMeans(imps)
imp <- imp / sum(imp)
imp <- imp[order(imp, decreasing=T)]
imp_df <- data.frame(names=factor(names(imp), levels=names(imp)), imp)
row.names(imp_df) <- NULL

ggplot(imp_df[1:15, ], aes(names, imp)) + 
  geom_col(fill='cadetblue3') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=32),
        axis.text.y = element_text(size=32),
        axis.title  = element_blank(),
        plot.title  = element_text(hjust=0.5))
head(imp, 30)

global_tree <- rpart(label ~ ., descriptors_df)
cols <- list("chartreuse", "khaki1", "burlywood1", "lightcoral", "hotpink",
             "darkolivegreen3", "cornflowerblue", "lightcyan2", "cornsilk", "cyan")
prp(global_tree, extra=2, varlen=0, faclen=0, box.palette=cols, 
    branch.type=5, split.border.col=2, split.col="red", cex=1.8, 
    legend.x=NA, compress=T, space=1, gap=0) 








