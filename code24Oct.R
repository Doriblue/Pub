library(virtualspecies)
library(rnaturalearth)
library(tidyverse)
library(fuzzySim)
library(mgcv)
library(ranger)
library(PresenceAbsence)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggstatsplot)
library(gghalves)
library(biomod2)
library(readr)
library(vegan)
library(FSA)
library(raster)
library(RStoolbox)
library(caTools)
library(car)
library(quantmod)
library(MASS)
library(corrplot)
library(rstatix)
library(groupdata2)
library(ROCR)
library(raster)
library(ranger)
library(dplyr)
library(sp)
library(virtualspecies)
library(corrplot)
library(ggsignif)

# Set up environment layers
filenames <- list.files(path = "C:/Users/fpt/Desktop/Publication/envi_vari/", pattern = "wc2.1_10m_bio_", full.names = TRUE)
raster_list <- lapply(filenames, stack)
rasters <- stack(raster_list)
Asia <- ne_countries(scale="medium", type="map_units", returnclass="sf", country=c('Myanmar','Thailand','Laos','Cambodia','Vietnam'))
rasters <- mask(rasters, Asia)
names(rasters) <- paste0("bio_", 1:nlayers(rasters))
ex <- extent(90, 112, 4, 30)
rasters <- crop(rasters, ex)

# Multicollinearity check
df <- as.data.frame(rasters, xy = TRUE) %>% na.omit()
model <- lm(x ~ ., data = df)

df_x <- df[,3:21]                                       # Independent variables 
var <- cor(df_x)                                         # Independent variables correlation matrix 
var_inv <- ginv(var)                                     # Inverse correlation matrix 
colnames(var_inv) <- colnames(df_x)                      # Rename row and column names
rownames(var_inv) <- colnames(df_x)
corrplot(var_inv, method = 'number', is.corr = F)        # Visualize multicollinearity

# Remove highly correlated layers
layers_to_remove <- c('bio_3', 'bio_4', 'bio_5', 'bio_6', 'bio_13', 'bio_14')
for (layer in layers_to_remove) {
  rasters <- dropLayer(rasters, layer)
}
print(rasters)
#myRandNum <- sample(1:nlayers(rasters), size = 5, replace = FALSE)
myRandNum <- c(9,12,3,13,7)

no.species<- 5
# ranger_sp function definition
# Modified ranger_sp function remains unchanged
ranger_sp <- function(x, TrainValue = 0.8, TestValue = 0.2) {
  library(groupdata2)
  library(ROCR)
  set.seed(999)
  inds <- partition(x, p = c(train = TrainValue, test = TestValue))
  train <- as.data.frame(inds[1])
  test <- as.data.frame(inds[2])
  validation <- test
  training <- train
  Model <- ranger(training$pres ~ ., data = training, importance = 'impurity')
  return(list(Model = Model, TestData = test))
}

# Modified Random.SDM function to handle flexible sample prevalence
Random.SDM <- function(alpha = -0.05, value = NULL, Predictors = NULL, SamplePoints = NULL, species_prevalence = NULL, seed = 999) {
  set.seed(seed)
  # Select the same set of 5 raster layers from Predictors for all iterations
  Outputs.l <- list()
  
  for (i in 1:no.species) {  # assuming no.species is predefined
    set.seed(seed)
    print(paste("Generating random species for iteration", i))
    
    # Use the same set of predictors for each iteration
    random.sp <- virtualspecies::generateRandomSp(Predictors[[myRandNum]],
                                                  convert.to.PA = FALSE,
                                                  species.type = "additive",
                                                  realistic.sp = TRUE,
                                                  plot = FALSE)
    set.seed(seed)
    new.pres <- convertToPA(random.sp, 
                            beta = "random",
                            alpha = alpha, 
                            plot = FALSE, 
                            species.prevalence = species_prevalence)
    
    print("Sampling occurrences")
    available_cells <- sum(!is.na(values(Predictors[[myRandNum[1]]])))  # count available cells in the first predictor
    print(paste("Available cells:", available_cells))
    
    if (SamplePoints > available_cells) {
      stop("SamplePoints exceeds the number of available cells")
    }
    
    presence.points <- sampleOccurrences(new.pres,
                                         n = SamplePoints, 
                                         type = "presence-absence",
                                         sample.prevalence = value,
                                         detection.probability = 1,
                                         correct.by.suitability = FALSE,
                                         plot = FALSE)
    
    PresAbs <- presence.points$sample.points[, c("x", "y", "Observed")]
    coordinates(PresAbs) <- ~x + y 
    crs(PresAbs) <- crs(Predictors) 
    values <- raster::extract(Predictors[[myRandNum]], PresAbs, df = TRUE)
    values <- values[,-1]
    
    modSpecies <- data.frame(pres = PresAbs@data[,1], values[1:ncol(values)])
    
    preds <- as.data.frame(Predictors[[myRandNum]]) %>% drop_na()
    predsXY <- as.data.frame(Predictors[[myRandNum]], xy = TRUE) %>% drop_na()
    
    # Use ranger_sp function to train the model
    Model <- ranger_sp(modSpecies, TrainValue = 0.8, TestValue = 0.2)
    Pred <- predict(Model$Model, data = preds, predict.all = FALSE, num.trees = Model$Model$num.trees)
    
    prediction <- data.frame(predsXY[,1:2], Pred$predictions)
    Pred <- rasterFromXYZ(prediction, crs = "+proj=longlat +datum=WGS84 +no_defs")
    
    Outputs <- list("Mods" = Model, "Raster" = Pred, "ModelDatabase" = modSpecies, "Occurrences" = PresAbs, "TestData"= Model$TestData)
    Outputs.l[[i]] <- Outputs 
  }
  
  return(Outputs.l)
}

# Sample prevalence and species prevalence lists
sample_prevalence_list <- c(0.05, 0.1, 0.2, 0.5)
species_prevalence_list <- c(0.5, 0.03)
non_overlapping_time_period_list <- c("t1", "t2", "t3", "t4")
overlapping_time_period_list <- c("T1", "T2", "T3", "T4")
replicate_list <- 1:5  # Assuming 50 replicates
repetitions<-5
# Initialize a list to store the down-sampling sizes for non-overlapping periods
downsampling_sizes <- list()

# Function to run models and store down-sampling sizes
run_model_and_store_sizes <- function(sample_prevalence, species_prevalence, sample_size) {
  # Run the model with the provided sample size
  result <- Random.SDM(
    alpha = -0.05,
    value = sample_prevalence,
    species_prevalence = species_prevalence,
    Predictors = rasters,
    SamplePoints = sample_size
  )
  
  return(result)
}

# Initialize the Results list
Results <- list()

# Loop through each species prevalence
for (species_prevalence in species_prevalence_list) {
  species_results <- list()  # Initialize results for this species
  
  # Loop through each sample prevalence
  for (sample_prevalence in sample_prevalence_list) {
    sample_results <- list()  # Initialize results for this sample prevalence
    # Loop through each replicate
    for (replicate in 1:repetitions) {
      # Use the appropriate sample size for each time period
      # Adjust accordingly if needed for other time periods
        sample_size <- 1000
    
      
      # Run the model and store the result
      result <- run_model_and_store_sizes(sample_prevalence, species_prevalence, sample_size)
      
      # Store the result for this replicate
    }
    
    # Add the sample results to the species results
    species_results[[paste("SamplePrevalence", sample_prevalence, sep = "_")]] <- result
  }
  
  # Store results for this species prevalence
  Results[[paste("SpeciesPrevalence", species_prevalence, sep = "_")]] <- species_results
}


# Function to accumulate occurrences for the same replicate across sample prevalence levels
accumulate_occurrences <- function(results, species_prevalence, sample_prevalence, replicate) {
  accumulated_occurrences <- NULL
  accumulated_coords <- NULL
  total_sample_size <- 0
  
  # Check if the specific species prevalence exists in results
  species_key <- paste0("SpeciesPrevalence_", species_prevalence)
  sample_key <- paste0("SamplePrevalence_", sample_prevalence)
  
  # Ensure the species and sample prevalence exist in the results
  if (!species_key %in% names(results)) {
    stop(paste("SpeciesPrevalence", species_prevalence, "not found in results"))
  }
  if (!sample_key %in% names(results[[species_key]])) {
    stop(paste("SamplePrevalence", sample_prevalence, "not found in results for SpeciesPrevalence", species_prevalence))
  }
  
  # Loop through the results to accumulate presence points and coordinates
  for (i in seq_along(results[[species_key]][[sample_key]])) {
    # Extract the 'ModelDatabase' for the current replicate
    model_data <- results[[species_key]][[sample_key]][[replicate]][["ModelDatabase"]]
    
    # Check if 'ModelDatabase' exists and is not NULL
    if (!is.null(model_data)) {
      # Filter for presence points (pres == 1)
      pres_occurrences <- model_data[model_data[["pres"]] == 1, ]
      
      # Accumulate presence occurrences
      if (is.null(accumulated_occurrences)) {
        accumulated_occurrences <- pres_occurrences
      } else {
        accumulated_occurrences <- rbind(accumulated_occurrences, pres_occurrences)
      }
      
      # Now extract the coordinates from the 'Occurrences' object
      coords <- results[[species_key]][[sample_key]][[replicate]][["Occurrences"]]@coords
      
      # Accumulate only the coordinates for presence points
      pres_coords <- coords[model_data[["pres"]] == 1, ]
      
      if (is.null(accumulated_coords)) {
        accumulated_coords <- pres_coords
      } else {
        accumulated_coords <- rbind(accumulated_coords, pres_coords)
      }
      
      # Calculate the total sample size (assuming a base size of 1000)
      total_sample_size <- 1000
      
      # Randomly downsample from accumulated occurrences and coordinates if necessary
      if (nrow(accumulated_occurrences) > total_sample_size) {
        # Create an index to sample occurrences and coordinates in sync
        sampled_indices <- sample(nrow(accumulated_occurrences), total_sample_size)
        
        accumulated_occurrences <- accumulated_occurrences[sampled_indices, ]
        accumulated_coords <- accumulated_coords[sampled_indices, ]
      }
    } else {
      warning(paste("ModelDatabase not found for SpeciesPrevalence", species_prevalence, "SamplePrevalence", sample_prevalence, "Replicate", replicate))
    }
  }
  
  return(list(occurrences = accumulated_occurrences, coords = accumulated_coords))
}


# Initialize the OverlappingResults list
OverlappingResults <- list()

# Loop through each species prevalence
for (species_prevalence in species_prevalence_list) {
  species_overlapping_results <- list()  # Initialize results for this species
  
  # Loop through each sample prevalence
  for (sample_prevalence_index in 1:length(sample_prevalence_list)) {
    sample_prevalence <- sample_prevalence_list[sample_prevalence_index]
    overlapping_sample_results <- list()  # Initialize results for this sample prevalence
    
    # Calculate cumulative sample prevalence for overlapping time periods
    cumulative_sample_prevalence <- cumsum(sample_prevalence_list[1:sample_prevalence_index])
    
    # Loop through each replicate
    for (replicate in 1:repetitions) {
      # Initialize accumulated occurrences for overlapping time periods
      accumulated_occurrences <- NULL
      
      # Use the cumulative sample prevalence for the current time period
      sample_size <- round(1000 * cumulative_sample_prevalence[sample_prevalence_index])  # Adjust sample size based on cumulative prevalence
      
      # Accumulate occurrences from previous non-overlapping periods
      accumulated_occurrences <- accumulate_occurrences(
        Results, species_prevalence, sample_prevalence, replicate
      )
      
      # Select occurrences for the overlapping period from the accumulated set
      selected_occurrences <- accumulated_occurrences$coords
      presence_count<- nrow(selected_occurrences)
      
      # Calculate number of absences needed
      absences_needed <- 1000 - sample_size
      
      # Prepare to store results
      Outputs <- list()
      
      # Check if absences are needed
        # Generate random absences using BIOMOD_FormatingData
        response_var <- rep(1, times = sample_size)
        pres_coords <- selected_occurrences[sample(1:presence_count, sample_size), ] #select randomly the number of presences within the presence pool that was accumulated
        abs_coords <- matrix(NA, ncol = 2, nrow = absences_needed)
        colnames(abs_coords)<- colnames(pres_coords)
        all_coords <- rbind(pres_coords, abs_coords)
        resp_name <- paste0(
          "SpeciesPrevalence_", species_prevalence, 
          "_SamplePrevalence_", sample_prevalence, 
          "_Replicate_", replicate
        )
        
        # Extract the independent variable names from the model
        independent_variable_names <- Results[[paste0("SpeciesPrevalence_", species_prevalence)]][[paste0("SamplePrevalence_", sample_prevalence)]][[replicate]][["Mods"]][["Model"]][["forest"]][["independent.variable.names"]]
        
        # Create the expl.var variable as a string containing the species prevalence, sample prevalence, and replicate
        expl_var <- rasters[[independent_variable_names]]
        expl_var <- stack(expl_var)
        
        # Now use expl_var in your BIOMOD_FormatingData call
        formatted_data <- BIOMOD_FormatingData(
          resp.name = resp_name,              # Dynamic response name
          resp.var = response_var,            # The presence/absence data
          expl.var = expl_var,                # Raster layers as predictors
          resp.xy = pres_coords,              # Coordinates of presence and placeholder absence points
          PA.nb.rep = 1,                      # One absence dataset
          PA.nb.absences = absences_needed,   # Number of absence points needed
          PA.strategy = "random"              # Random sampling strategy for pseudo-absences
        )
        
        # Extract absences and combine with presences
        final_occurrences <- cbind(formatted_data@data.species, formatted_data@data.env.var)
        colnames(final_occurrences)[1] <- "pres"
        final_occurrences$pres[is.na(final_occurrences$pres)] <- 0
        
        # Run the SDM model using ranger_sp
        modSpecies <- final_occurrences  # Assuming final_occurrences contains presence/absence data with predictors
        modSpecies$pres <- as.factor(modSpecies$pres)
        Model <- ranger_sp(modSpecies, TrainValue = 0.8, TestValue = 0.2)
        
        # Use the model to make predictions
        preds <- as.data.frame(expl_var) %>% drop_na()  # Exclude the response column by name
        predsXY <- as.data.frame(expl_var, xy = TRUE) %>% drop_na()  # Assuming the first two columns are coordinates
        Pred <- predict(Model$Model, data = preds, predict.all = FALSE, num.trees = Model$Model$num.trees)
        
        # Create a raster from the predictions
        prediction <- data.frame(predsXY[, 1:2], Pred$predictions)
        prediction <- prediction %>% distinct(x, y, .keep_all = TRUE)
        grid <- expand.grid(x = seq(min(prediction$x), max(prediction$x), by = 0.1),
                            y = seq(min(prediction$y), max(prediction$y), by = 0.1))
        
        # Join the grid with your prediction data to find missing values
        full_prediction <- grid %>%
          left_join(prediction, by = c("x", "y"))
        
        # Create raster from the complete prediction data
        PredRaster <- rasterFromXYZ(full_prediction[, c("x", "y", "Pred.predictions")], 
                                    crs = "+proj=longlat +datum=WGS84 +no_defs")
        
        # Store the outputs for this replicate
        Outputs <- list(
          "Mods" = Model,
          "Raster" = PredRaster,
          "ModelDatabase" = modSpecies,
          "Occurrences" = final_occurrences,
          "TestData" = Model$TestData
        )
      
      # Store the result for this replicate
        overlapping_sample_results[[replicate]]<- Outputs
        }
    
    # Add the sample results to the species results
    species_overlapping_results[[paste("SamplePrevalence", sample_prevalence, sep = "_")]] <- overlapping_sample_results
  }
  
  # Store results for this species prevalence
  OverlappingResults[[paste("SpeciesPrevalence", species_prevalence, sep = "_")]] <- species_overlapping_results
}





# The overlapping_results now has the same structure as the non-overlapping Results list
predictions_list <- list()

# Loop through each species prevalence
for (species_prevalence in species_prevalence_list) {
  
  # Loop through each sample prevalence
  for (i in seq_along(sample_prevalence_list)) {
    sample_prevalence <- sample_prevalence_list[i]
    
    # Calculate cumulative sample prevalence up to the current sample prevalence
    cumulative_sample_prevalence <- sum(sample_prevalence_list[1:i])
    
    # Loop through each replicate for Non-overlapping temporal setting
    for (replicate in replicate_list) {
      # Access the result for Non-overlapping setting
      result_NO <- Results[[paste0("SpeciesPrevalence_", species_prevalence)]][[paste0("SamplePrevalence_", sample_prevalence)]][[replicate]]
      
      # Check if result_NO is not NULL and contains predictions
      if (!is.null(result_NO) && "Mods" %in% names(result_NO) && "Model" %in% names(result_NO[["Mods"]]) && "predictions" %in% names(result_NO[["Mods"]][["Model"]])) {
        predictions_no <- result_NO[["Mods"]][["Model"]][["predictions"]]
        
        # Ensure predictions_no is a valid vector
        if (length(predictions_no) > 0) {
          # Append the predictions to the list for Non-overlapping setting
          predictions_list[[length(predictions_list) + 1]] <- data.frame(
            Prediction = predictions_no,
            TemporalSetting = "Non-overlapping",
            SamplePrevalence = sample_prevalence,
            SpeciesPrevalence = species_prevalence
          )
        } else {
          cat("Empty predictions for Non-overlapping. Species:", species_prevalence, "Sample Prevalence:", sample_prevalence, "Replicate:", replicate, "\n")
        }
      } else {
        cat("No predictions found for Non-overlapping. Species:", species_prevalence, "Sample Prevalence:", sample_prevalence, "Replicate:", replicate, "\n")
      }
      
      # Access the result for Overlapping setting
      result_O <- OverlappingResults[[paste0("SpeciesPrevalence_", species_prevalence)]][[paste0("SamplePrevalence_", sample_prevalence)]][[replicate]]
      
      # Check if result_O is not NULL and contains predictions
      if (!is.null(result_O) && "Mods" %in% names(result_O) && "Model" %in% names(result_O[["Mods"]]) && "predictions" %in% names(result_O[["Mods"]][["Model"]])) {
        predictions_o <- result_O[["Mods"]][["Model"]][["predictions"]]
        
        # Ensure predictions_o is a valid vector
        if (length(predictions_o) > 0) {
          # Append the predictions to the list for Overlapping setting
          predictions_list[[length(predictions_list) + 1]] <- data.frame(
            Prediction = predictions_o,
            TemporalSetting = "Overlapping",
            SamplePrevalence = cumulative_sample_prevalence,
            SpeciesPrevalence = species_prevalence
          )
        } else {
          cat("Empty predictions for Overlapping. Species:", species_prevalence, "Sample Prevalence:", sample_prevalence, "Replicate:", replicate, "\n")
        }
      } else {
        cat("No predictions found for Overlapping. Species:", species_prevalence, "Sample Prevalence:", sample_prevalence, "Replicate:", replicate, "\n")
      }
    }
  }
}

# Convert predictions_list to a data frame
if (length(predictions_list) > 0) {
  predictions_df <- do.call(rbind, predictions_list)
} else {
  predictions_df <- data.frame()  # Create an empty data frame if no predictions found
}

# Check the final predictions_df
print(predictions_df)

as.factor(predictions_df$SamplePrevalence)

write.csv(predictions_df, file = "predictions_dataframe.csv")


#####Evaluation metrics#####
# Initialize an empty list to store the evaluation results
Mer_Results <- list()

# Merge non-overlapping and overlapping results
for (species_idx in seq_along(species_prevalence_list)) {
  species_prevalence <- species_prevalence_list[species_idx]
  
  # Create a sublist for Non-overlapping and Overlapping settings
  temporal_results <- list(
    "Non-overlapping" = list(),
    "Overlapping" = list()
  )
  
  # Loop through each sample prevalence
  for (sample_idx in seq_along(sample_prevalence_list)) {
    sample_prevalence <- sample_prevalence_list[sample_idx]
    
    # Initialize lists for Non-overlapping and Overlapping results by sample prevalence
    temporal_results[["Non-overlapping"]][[paste("SamplePrevalence_", sample_prevalence, sep = "")]] <- list()
    temporal_results[["Overlapping"]][[paste("SamplePrevalence_", sample_prevalence, sep = "")]] <- list()
    
    # Loop through each replicate
    for (replicate in replicate_list) {
      
      # Access Non-overlapping result and store in the list
      result_NO <- Results[[paste0("SpeciesPrevalence_", species_prevalence)]][[paste0("SamplePrevalence_", sample_prevalence)]][[replicate]]
      
      if (!is.null(result_NO)) {
        temporal_results[["Non-overlapping"]][[paste0("SamplePrevalence_", sample_prevalence)]][[paste0("Replicate_", replicate)]] <- result_NO
      }
      
      # Access Overlapping result and store in the list
      result_O <- OverlappingResults[[paste0("SpeciesPrevalence_", species_prevalence)]][[paste0("SamplePrevalence_", sample_prevalence)]][[replicate]]
      
      if (!is.null(result_O)) {
        temporal_results[["Overlapping"]][[paste0("SamplePrevalence_", sample_prevalence)]][[paste0("Replicate_", replicate)]] <- result_O
      }
    }
  }
  
  # Add the species-specific temporal results to Mer_Results
  Mer_Results[[paste("SpeciesPrevalence_", species_prevalence, sep = "")]] <- temporal_results
}



library(ROCR)
library(ecospat)

evaluate_model <- function(model, test_data) {
  # Extract the predictors (environmental variables) and response (P/A) from test_data
  test_predictors <- test_data[, -which(names(test_data) == "pres")]
  test_response <- test_data$pres
  test_response <- as.numeric(as.character(test_response))
  
  
  # Predict using the test data
  pred <- predict(model, data = test_predictors, predict.all = F,num.trees=model$num.trees)  # Ensure predict.all is FALSE for standard predictions
  pred$predictions <- as.numeric(as.character(pred$predictions))
  
  # Check if predictions are available and of correct length
  if (is.null(pred$predictions) || length(pred$predictions) != length(test_response)) {
    stop("Mismatch in length between predictions and test responses or predictions are NULL.")
  }
  
  # Calculate AUC using ROCR
  pred_obj <- prediction(pred$predictions, test_response)
  auc_obj <- performance(pred_obj, measure = "auc")
  auc_value <- auc_obj@y.values[[1]]
  auc_value<- round(auc_value,4)
  
  # Calculate TSS
  df.prediction<- pred$predictions
  # data_thres<- data.frame(ID=1:nrow(test_data),test_data[,1], df.prediction)
  # threshold <- PresenceAbsence::optimal.thresholds(DATA= data_thres, opt.methods = 'MaxSens+Spec')
  # threshold<- threshold[,2]
  # predicted_binary <- ifelse(pred$predictions >= threshold, 1, 0)
  cmx<- table(Predicted = df.prediction, Actual = test_data[,1])
  sensitivity <- cmx[2, 2] / sum(cmx[,2 ])  # True Positive / (True Positive + False Negative)
  specificity <- cmx[1, 1] / sum(cmx[,1 ])  # True Negative / (True Negative + False Positive)
  
  # Calculate TSS
  tss_value <- sensitivity + specificity - 1
  tss_value <- round(tss_value, 4)
  
  
  return(list(AUC = auc_value, TSS = tss_value))
  
}
Evaluation_Results <- list()

# Loop through each species prevalence
for (species in names(Mer_Results)) {
  species_results <- Mer_Results[[species]]

  # Initialize a list to store results for this species
  species_eval <- list()

  # Loop through temporal settings (e.g., different time periods)
  for (temporal_setting in names(species_results)) {
    temporal_results <- species_results[[temporal_setting]]

    # Initialize a list to store results for this temporal setting
    temporal_eval <- list()

    # Loop through sample prevalence
    for (sample in names(temporal_results)) {
      sample_prevalence <- temporal_results[[sample]]

      # Initialize a list to store results for this sample prevalence
      sample_eval <- list()

      # Loop through each replicate
      for (replicate in 1:repetitions) {
        result <- sample_prevalence[[replicate]]

        # Access the model and test data
        model <- result$Mods$Model  # Access the model
        test_data <- result$Mods$TestData  # Access the test data

        # Ensure that Predictions are not NULL
      
          # Calculate evaluation metrics
          eval_metrics <- evaluate_model(model, test_data)

          # Store the evaluation results for this replicate
          sample_eval[[replicate]] <- eval_metrics
       
      }

      # Store the sample evaluation results for this sample prevalence
      temporal_eval[[sample]] <- sample_eval
    }

    # Store the temporal evaluation results for this temporal setting
    species_eval[[temporal_setting]] <- temporal_eval
  }

  # Store the species evaluation results
  Evaluation_Results[[species]] <- species_eval
}


# Initialize an empty list to store results
extracted_results <- list()

# Loop through species in Evaluation_Results
for (species in names(Evaluation_Results)) {
  species_eval <- Evaluation_Results[[species]]
  
  # Loop through temporal settings (e.g., different time periods)
  for (temporal_setting in names(species_eval)) {
    temporal_eval <- species_eval[[temporal_setting]]
    
    # Loop through sample prevalence
    for (sample in names(temporal_eval)) {
      sample_eval <- temporal_eval[[sample]]
      
      # Loop through each replicate
      for (replicate in 1:length(sample_eval)) {
        eval_metrics <- sample_eval[[replicate]]
        
        # Extract the AUC and TSS values
        auc_value <- eval_metrics$AUC
        tss_value <- eval_metrics$TSS
        
        # Create a data frame row with the relevant information
        extracted_results <- rbind(extracted_results, 
                                   data.frame(
                                     Species = species,
                                     Temporal_Setting = temporal_setting,
                                     Sample_Prevalence = sample,
                                     Replicate = replicate,
                                     AUC = auc_value,
                                     TSS = tss_value,
                                     stringsAsFactors = FALSE
                                   ))
      }
    }
  }
}

# Convert the list to a data frame
evaluation_df <- as.data.frame(extracted_results)

# View the resulting data frame
head(evaluation_df)
evaluation_df<- unique(evaluation_df[,-4])
evaluation_df<- evaluation_df[-c(5,9),]

# Apply the function to create the TimePeriod column
evaluation_df$TimePeriod <- rep(c("t1", "t2", "t3", "t4", "T1","T2", "T3", "T4" ), times=2)
evaluation_df$SamplePrevalence<- rep(c("0.05", "0.1", "0.2", "0.5", "0.05","0.15", "0.35", "0.85" ), times=2)
evaluation_df$SpeciesPrevalence <- rep(c("common", "Rare"), times=8)
# Apply the function to create the new SamplePrevalence column
evaluation_df$TimePeriod <- as.character(evaluation_df$TimePeriod)
evaluation_df$TimePeriod <- as.factor(evaluation_df$TimePeriod)
write.csv(evaluation_df, "evaluation_df.csv", row.names=F)

# Print the dataframe
# Load ggplot2
library(ggplot2)
evaluation_df<- read.csv("evaluation_df.csv")
evaluation_df$TimePeriod <- as.factor(evaluation_df$TimePeriod)
evaluation_df$SamplePrevalence<- as.factor(evaluation_df$SamplePrevalence)
evaluation_df$SpeciesPrevalence<- as.factor(evaluation_df$SpeciesPrevalence)
evaluation_df$Temporal_Setting<- as.factor(evaluation_df$Temporal_Setting)
str(evaluation_df)
# Plot AUC
AUC<-ggplot(evaluation_df, aes(x = Temporal_Setting, y = AUC, shape = SpeciesPrevalence, color=TimePeriod)) +
  geom_point(size=3) +
  labs(title = "AUC",
       x = "Temporal Setting",
       y = "AUC") +
  theme_minimal()+
  theme(
    legend.background = element_blank(),
    legend.position = 'bottom',
    legend.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.text.x = element_text(size = 10, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.text.y = element_text(size = 10, face = 'bold'),
    axis.ticks.y = element_blank(),
    text = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12, angle = 0),
    legend.key.size = unit(1.0, 'cm')
  )

# Plot TSS
TSS<-ggplot(evaluation_df, aes(x = Temporal_Setting, y = TSS, shape = SpeciesPrevalence, color=TimePeriod)) +
  geom_point(size=3) +
  labs(title = "TSS",
       x = "Temporal Setting",
       y = "TSS") +
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(size = 12, face = 'bold'),
        axis.text.y = element_text(size = 10, face = 'bold'),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12),
        strip.text = element_text(size = 12)        )
com<- AUC+TSS
ggsave(plot=com, filename = "Performance.png", width = 12, height = 5, dpi=800)





##################EXTRACT RASTERS##############
# Initialize a list to store the extracted rasters
# Initialize the main list to store extracted rasters
Extracted_Rasters <- list()

# Loop through each species in Mer_Results
for (species_key in names(Mer_Results)) {
  species_results <- Mer_Results[[species_key]]
  
  # Initialize a list to store results for this species
  species_rasters <- list()
  
  # Loop through temporal settings (e.g., different time periods)
  for (temporal_setting_key in names(species_results)) {
    temporal_results <- species_results[[temporal_setting_key]]
    
    # Initialize a list to store results for this temporal setting
    temporal_setting <- list()
    
    # Loop through sample prevalence
    for (sample_key in names(temporal_results)) {
      sample_prevalence <- temporal_results[[sample_key]]
      
      # Initialize a list to store results for this sample prevalence
      sample_prevalence_rasters <- list()
      
      # Loop through each replicate
      for (replicate_key in names(sample_prevalence)) {
        result <- sample_prevalence[[replicate_key]]
        
        # Extract the raster
        raster_prediction <- result[["Raster"]]
        
        # Create a unique key for the raster
        raster_key <- paste(species_key, temporal_setting_key, sample_key, replicate_key, sep = "_")
        
        # Store the raster in the sample prevalence list
        sample_prevalence_rasters[[raster_key]] <- raster_prediction
      }
      
      # Store the sample prevalence results for this sample key
      temporal_setting[[sample_key]] <- sample_prevalence_rasters
    }
    
    # Store the temporal evaluation results for this temporal setting
    species_rasters[[temporal_setting_key]] <- temporal_setting
  }
  
  # Store the species-specific rasters in the main Extracted_Rasters list
  Extracted_Rasters[[species_key]] <- species_rasters
}


# Initialize a list to store the CV results
library(raster)

# Initialize a list to store CV results
library(raster)

# Initialize a list to store CV results

CV_Results <- list()

# Loop through each species in Extracted_Rasters
for (species_key in names(Extracted_Rasters)) {
  species_rasters <- Extracted_Rasters[[species_key]]
  
  # Initialize a list to store results for this species
  species_cv <- list()
  
  # Loop through each temporal setting
  for (temporal_setting_key in names(species_rasters)) {
    temporal_rasters <- species_rasters[[temporal_setting_key]]
    
    # Extract raster layers from the temporal_rasters object
    if (is.list(temporal_rasters)) {
      # Flatten the nested structure to get all raster layers
      raster_layers <- unlist(lapply(temporal_rasters, function(x) {
        if (is.list(x)) {
          # If it's a list, return its raster layers
          return(unlist(x[sapply(x, inherits, "RasterLayer")], recursive = FALSE))
        } else if (inherits(x, "RasterLayer")) {
          # If it's already a raster layer, return it
          return(list(x))
        }
        return(NULL)  # Return NULL for non-raster elements
      }), recursive = FALSE)
      
      # Check if we have valid raster layers
      if (length(raster_layers) > 0) {
        # Create a raster stack
        raster_stack <- stack(raster_layers)
        
        # Calculate the CV (coefficient of variation)
        cv_raster <- cv(raster_stack)
        
        # Store the CV results for this temporal setting
        species_cv[[temporal_setting_key]] <- cv_raster
      } else {
        warning(paste("No valid raster layers found for", species_key, "and", temporal_setting_key))
      }
    } else {
      warning(paste("Temporal rasters for", species_key, "and", temporal_setting_key, "are not a list."))
    }
  }
  
  # Store the species-specific CV results in the main CV_Results list
  CV_Results[[species_key]] <- species_cv
}





# Initialize a list to store the difference in CV results
library(raster)

DiffCV_Results <- list()
# Loop through each species prevalence in CV_Results
for (species_key in names(CV_Results)) {
  species_rasters <- CV_Results[[species_key]]
  
  # Check if both raster layers exist
  if (all(c("Overlapping", "Non-overlapping") %in% names(species_rasters))) {
    # Resample one raster to match the other
    raster_non_overlap_resampled <- resample(species_rasters[["Non-overlapping"]],
                                             species_rasters[["Overlapping"]],
                                             method = "bilinear")  # Use "bilinear" or "ngb" depending on your data
    
    # Now perform the subtraction
    diff_cv <- species_rasters[["Overlapping"]] - raster_non_overlap_resampled
    
    # Convert the raster to a data frame and clean up
    df_diff_cv <- as.data.frame(diff_cv, xy = TRUE) %>% drop_na()
    colnames(df_diff_cv) <- c("x", "y", "values")
    
    # Store the data frame in the list
    DiffCV_Results[[paste(species_key, "_DiffCV", sep = "")]] <- df_diff_cv
  } else {
    warning(paste("One of the rasters does not exist for species:", species_key))
  }
}

DiffCV_Results$SpeciesPrevalence_0.5_DiffCV$Species<- "Common species"
DiffCV_Results$SpeciesPrevalence_0.03_DiffCV$Species<- "Rare species"
# Combine all the individual data frames into one
df_diffTot <- do.call(rbind, DiffCV_Results)

# Create the half-violin plot
diffcv_plot <- ggplot(df_diffTot, aes(Species, values)) +
  geom_half_violin(alpha = 0.6, side = "l", fill = "#0072B2") +
  geom_half_boxplot(nudge = 0.05, outlier.color = NA, side = "r", fill = "#0072B2") +
  geom_hline(aes(yintercept = 0), color = "gray70", size = 0.6) +
  labs(x = "Species", y = "Difference value between CVs (over - nonover)") +
  theme_light() +
  theme(
    legend.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    legend.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.text.x = element_text(size = 10, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.text.y = element_text(size = 10, face = 'bold'),
    axis.ticks.y = element_blank(),
    text = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12, angle = 0),
    legend.key.size = unit(1.0, 'cm')
  )

# Save the plot
ggsave(plot = diffcv_plot,
       filename = "diff_cv.jpeg",
       width = 12,
       height = 8,
       dpi = 600)




########################Statistics#################
predictions_df<- read.csv("C:/Users/fpt/Desktop/Publication/Publication/predictions_dataframe.csv")

kw_timeperiod<- kruskal.test(predictions_df$Prediction, predictions_df$TimePeriod )
kw_temporalsetting<- kruskal.test(predictions_df$Prediction, predictions_df$TemporalSetting )

dunn_timeperiod<- dunnTest(predictions_df$Prediction~predictions_df$TimePeriod)
dunn_temporal <- dunnTest(predictions_df$Prediction, predictions_df$TemporalSetting, method = "bonferroni")

# Load necessary libraries
library(ggstatsplot)

# Assuming predictions_df is your dataframe
ggbetweenstats(
  data = predictions_df,
  x = TemporalSetting,               # Grouping variable
  y = Prediction,                    # Response variable
  title = "Predictions by Temporal Setting", # Title of the plot
  xlab = "Temporal Setting",         # X-axis label
  ylab = "Prediction",               # Y-axis label
  plot.type = "boxplot",             # Type of plot
  pairwise.display = "p-value",      # Show pairwise comparisons p-values
  var.equal = TRUE                   # Assuming equal variances
)

ggbetweenstats(
  data = predictions_df,
  x = TimePeriod,               # Grouping variable
  y = Prediction,                     # Response variable
  title = "Predictions by Time period", # Title of the plot
  xlab = "TimePeriod",          # X-axis label
  ylab = "Prediction",                # Y-axis label
  plot.type = "boxplot",              # Type of plot
  pairwise.display = "p-value",       # Show pairwise comparisons p-values
  ggtheme = theme_minimal(),          # Theme for the plot
  var.equal = TRUE                    # Assuming equal variances
)
# Assuming 'TemporalSetting' has values such as "Overlapping" and "Non-overlapping"
# First, filter the dataset for the desired temporal setting (e.g., "Non-overlapping")

filtered_data <- predictions_df %>%
  filter(TemporalSetting == "Non-overlapping")  # Adjust to filter by the desired TemporalSetting

# Create the plot comparing predictions across time periods within the filtered temporal setting
ggbetweenstats(
  data = filtered_data,                     # Use the filtered dataset
  x = TimePeriod,                           # Grouping variable (TimePeriod)
  y = Prediction,                           # Response variable
  title = "Predictions by Time Period (Non-overlapping)",  # Adjust the title
  xlab = "Time Period",                     # X-axis label
  ylab = "Prediction",                      # Y-axis label
  plot.type = "boxplot",                    # Type of plot
  pairwise.display = "p-value",             # Show pairwise comparisons p-values
  ggtheme = theme_minimal(),                # Theme for the plot
  var.equal = TRUE                          # Assuming equal variances
)
# Assuming 'TemporalSetting' has values such as "Overlapping" and "Non-overlapping"
# First, filter the dataset for the desired temporal setting (e.g., "Non-overlapping")

filtered_data <- predictions_df %>%
  filter(TemporalSetting == "Overlapping")  # Adjust to filter by the desired TemporalSetting

# Create the plot comparing predictions across time periods within the filtered temporal setting
ggbetweenstats(
  data = filtered_data,                     # Use the filtered dataset
  x = TimePeriod,                           # Grouping variable (TimePeriod)
  y = Prediction,                           # Response variable
  title = "Predictions by Time Period (Overlapping)",  # Adjust the title
  xlab = "Time Period",                     # X-axis label
  ylab = "Prediction",                      # Y-axis label
  plot.type = "boxplot",                    # Type of plot
  pairwise.display = "p-value",             # Show pairwise comparisons p-values
  ggtheme = theme_minimal(),                # Theme for the plot
  var.equal = TRUE                          # Assuming equal variances
)

