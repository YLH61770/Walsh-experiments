data <- read.csv("data.csv", header = TRUE)

# 3:1 split ratio (75% training set, 25% test set)
train_ratio <- 0.75
n <- nrow(data)

# Generate random indices
train_indices <- sample(1:n, size = round(train_ratio * n))

# Split the data
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]

# Save the split datasets
write.csv(train_data, "data_train_split.csv", row.names = FALSE)
write.csv(test_data, "data_test_split.csv", row.names = FALSE)

# Check splitting results
cat("Original data size:", n, "\n")
cat("Training set size:", nrow(train_data), "\n")  
cat("Test set size:", nrow(test_data), "\n")
cat("Split ratio:", round(nrow(train_data)/n, 2), ":", round(nrow(test_data)/n, 2), "\n")