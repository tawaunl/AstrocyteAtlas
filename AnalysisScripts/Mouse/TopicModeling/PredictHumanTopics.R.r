# Load the fastTopics package
library(fastTopics)

# Assuming `fit` is the result of running fastTopics on the subset counts table
# and `counts_rest` is the counts table for the rest of the dataset.

# Predict topics for the rest of the dataset
predicted_topics <- predict(fit, counts_rest)

# View the predicted topic proportions
print(predicted_topics$L)

# Save the predicted topic proportions to a file if needed
write.csv(predicted_topics$L, "predicted_topics.csv", row.names = FALSE)