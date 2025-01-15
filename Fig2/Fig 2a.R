library(readr)
library(ggplot2)
library(gridExtra)

# Import data
data <- read.table("data_cf.txt", header = TRUE, sep = "\t")

# Determine the column name of the y variable
y_columns <- colnames(data)[2:9]
# Determine the column name of the x variable
x_columns <- colnames(data)[1:1]

# Create a blank graphic list
plots_list <- list()

# Loop for each y variable
for (y_column in y_columns) {
  for (x_column in x_columns) {
    # Use the lm() function to build a polynomial regression model
    formula_string <- paste(y_column, "~", x_column, "+ I(", x_column, "^2)")
    regression_model <- lm(formula_string , data = data)
    ## Extract the p value after adjusting R2
    summary_model <- summary(regression_model)
    f_statistic <- summary_model$fstatistic
    p_value_all <- pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE)
    p_value_label_all <- ifelse(p_value_all < 0.001, "***", ifelse(p_value_all < 0.01, "**", ifelse(p_value_all < 0.05, "*", "")))
    
    # Draw scatter plots, regression lines and curve fitting
    plot <- ggplot(data, aes_string(x = x_column, y = y_column)) +
      geom_smooth(method = "lm", formula = y ~ poly(x, degree = 2), color = "#652884", fill = ifelse(p_value_all < 0.05, "#F5BE8F", "#C1E0DB"), linetype = ifelse(p_value_all < 0.05, "solid", "dashed"), alpha = 0.5) +  # 添加多项式回归线及置信区间，并设置阴影透明度
      geom_point(color = "#652884", alpha = 0.5) +  # Change the scatter color to green and set transparency
      labs(x = "Distance (m)", y = paste(y_column," CO2 flux (μmol/m²/s)")) +  # Add axis label
      ylim(0.5, 5.5) +  # Set the ordinate range to start at 0
      theme_minimal() +  
      theme(panel.grid = element_blank(),  
            axis.title.x = element_text(size = 12, family = "serif", face = "bold"),  
            axis.title.y = element_text(size = 12, family = "serif", face = "bold"),  
            axis.text.x = element_text(size = 10, family = "serif"),  
            axis.text.y = element_text(size = 10, family = "serif"),  
            axis.line = element_line(color = "black"),  
            axis.ticks = element_line(color = "black")) +  
      geom_text(x = min(data[[x_column]]), y = 0.5, 
                label = paste("R² =", round(summary(regression_model)$adj.r.squared, 3), p_value_label_all),  # Add R² and significance tags
                hjust = 0, vjust = 0,  
                size = 4, family = "serif", face = "bold")  
    
    # Add the graph to the graph list
    plot_name <- paste(y_column, "_vs_", x_column, ".png", sep = "")
    plots_list[[plot_name]] <- plot
  }
}

# Combine all graphics into one graph using grid.arrange
combined_plot <- grid.arrange(grobs = plots_list, ncol = 3)

# Save the merged graphic
ggsave("CF.pdf", plot = combined_plot, device = "pdf", width = 9, height = 9, units = "in")

ggsave("CF.png", plot = combined_plot, device = "png", width = 9, height = 9, units = "in")

# Output the results of summary(regression_model) to a txt file. Only the last summary_model of the loop can be output here, so it needs to be calculated separately without using the for loop.
summary_txt <- capture.output(summary_model)
writeLines(summary_txt, "summary_regression_model.txt")
