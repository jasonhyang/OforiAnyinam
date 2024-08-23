library(qsmooth)

# Load Counts Table
mtb_counts <- read.csv("\\Mtb_counts.csv", header=TRUE)

# Group by condition
sample_labels <- c("WT1", "WT2", "WT3", "WT_BDQ1", "WT_BDQ2", "WT_BDQ3", 
                   "katG1", "katG2", "katG3", "katG_BDQ1”, "katG_BDQ2", "katG_BDQ3")

group_labels <- c(“WT", "WT", "WT", "WT_BDQ", "WT_BDQ", "WT_BDQ",
                  "katG", "katG", "katG", "katG_BDQ", "katG_BDQ", "katG_BDQ")

# Create a data frame with sample labels and their corresponding groups
sample_df <- data.frame(
  Sample = sample_labels,
  Group = group_labels
)
print(sample_df)

# Quantile Normalize
qs <- qsmooth(mtb_counts[,-1], group_labels)
dat <- log2(qsmoothData(qs) + 1)
rownames(dat) <- mtb_counts[,1]

write.csv(dat, file = "sample_groups.csv", row.names = FALSE)
