# Load necessary libraries
library(fastTopics)
library(ggplot2)
library(dplyr)
library(purrr)

# --- 1. Define the k values and model types to assess ---
# These should match the values used in your sbatch array job.
k_values <- c(2,3,4,6, 8, 10, 12, 14,16,18,20,22,26)
model_types <- c("poisson_nmf", "multinom_topic_model")
out.dir <- "/gpfs/scratchfs01/site/u/lucast3/AstrocyteAtlas/AnalysisScripts/Mouse/TopicModeling"

# --- 2. Load the data and create a results data frame ---
# Create an empty list to store all the results.
all_results <- list()
#counts <- GetAssayData(data,"RNA")
for (k in k_values) {
  filename <- file.path(out.dir,"results",paste0("fasttopics_k",k,"_results.rds"))
  
  if (file.exists(filename)) {
    message(paste("Loading file:", filename))
    fits <- readRDS(filename)
    
    for (model_type in model_types) {
      if (model_type %in% names(fits)) {
        fit <- fits[[model_type]]
        num_cells <- dim(fit$F)[1]
        num_genes <-dim(fit$L)[1]
        
        # Access log-likelihood and BIC directly from the fit object.
        loglik <- -fit$loss
        k_params <- (num_cells +num_genes) * k
        
        # Calculate the total number of data points.
        N <- as.numeric(num_cells) * as.numeric(num_genes)
        
        # Calculate BIC.
        bic <- -2 * loglik + k_params * log(N)
        # Store results in a list.
        result_row <- data.frame(
          k = k,
          model = model_type,
          log_likelihood = loglik,
          BIC = bic
        )
        
        all_results[[length(all_results) + 1]] <- result_row
      }
    }
  } else {
    warning(paste("File not found:", filename))
  }
}

# Combine the list of data frames into a single data frame.
results_df <- bind_rows(all_results)

# --- 4. Plot the results ---
# Plot 1: Log-likelihood vs. k
loglik_plot <- ggplot(results_df, aes(x = k, y = log_likelihood, color = model)) +
  geom_point(size = 3) +
  geom_line() +
  labs(
    title = "Model Log-likelihood vs. Number of Topics (k)",
    x = "Number of Topics (k)",
    y = "Log-likelihood",
    color = "Model Type"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = k_values)
print(loglik_plot)
ggsave(file.path(out.dir,"log_likelihood_vs_k.png"), plot = loglik_plot, width = 8, height = 6)

# Plot 2: BIC vs. k
bic_plot <- ggplot(results_df, aes(x = k, y = BIC, color = model)) +
  geom_point(size = 3) +
  geom_line() +
  labs(
    title = "Bayesian Information Criterion (BIC) vs. Number of Topics (k)",
    x = "Number of Topics (k)",
    y = "BIC",
    color = "Model Type"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = k_values)
print(bic_plot)
ggsave(file.path(out.dir,"bic_vs_k.png"), plot = bic_plot, width = 8, height = 6)


message("Assessment complete. Plots saved as PNG files.")