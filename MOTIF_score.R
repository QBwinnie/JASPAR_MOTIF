
# Assuming your file is named 'motif_results.csv'
motif_data <- read.csv('~/you_input.csv')

# Group by TF and calculate summary statistics
tf_summary <- aggregate(relScore ~ TF, data = motif_data, 
                        FUN = function(x) c(sum = sum(x),
                                            count = length(x),
                                            mean = mean(x),
                                            min = min(x),
                                            max = max(x)))

# Convert the results to a more readable format
tf_summary <- data.frame(
  TF = tf_summary$TF,
  sum_relScore = tf_summary$relScore[,1],
  count = tf_summary$relScore[,2],
  mean_relScore = tf_summary$relScore[,3],
  min_relScore = tf_summary$relScore[,4],
  max_relScore = tf_summary$relScore[,5]
)

# Round numeric columns to 3 decimal places
tf_summary[,2:6] <- round(tf_summary[,2:6], 3)

# Sort by total relScore in descending order
tf_summary <- tf_summary[order(-tf_summary$sum_relScore),]

# Print the summary
print(tf_summary)

# Print detailed information for a specific TF (e.g., MZF1)
print("\nDetailed sites for MZF1:")
mzf1_sites <- subset(motif_data8, TF == "MZF1(var.2)", 
                     select=c("start", "end", "relScore", "siteSeqs"))
print(mzf1_sites)

#Save summary to csv
write.csv (tf_summary,'your_output.csv')
