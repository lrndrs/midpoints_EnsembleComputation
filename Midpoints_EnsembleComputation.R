# ------------------------------------------------------------------
# Compute Palaeo midpoints with a Error-in-Variables approach, 
# Compute Overall Uncertainty by Drawing from MCMC ensemble
# 1) d18O
# 2) d13C_init
# 3) temporal relationship between indicators assuming a fully correlated age model

## This script is related to the publication
## Endres, L., C. Pérez-Mejı́as, R. Ivanovic, et al. “Interplay of North Atlantic Freshening and Deep Convection during the Last Deglaciation Constrained by Iberian Speleothems.” EGUsphere 2025 (2025): 1–37. https://doi.org/10.5194/egusphere-2025-3911.
## 
## License:
## This code is released under the MIT License.
## See the LICENSE file in this repository for details.
##
## If you use this code, please cite the associated publication.
# ------------------------------------------------------------------

# --------------------
# d18O
# --------------------
# Data and Library Loading, Depth definition
# --------------------

### Library loading
library(ggplot2)
library(plotly)
library(Bchron)
library(dplyr)
library(knitr)
library(Cairo)
library(openxlsx)

### Data Loading (i block for testing percentile)
GLA.pub <- readRDS("GLApub_withd13Cinit.rds") # Data / Author Age Model
a1 <- readRDS("Glas_Age2023_model_2024-04-22.rds") # Age Model Ensemble

### Assigned Measurement Uncertainty
sigma_y = 0.1 # permill (d18O, d13C)

### Find Depths for which the midpoint in between should be computed
### These have been picked manually using plotly, and for the adjacent computation have just also been changed manually
ggplotly(ggplot() + geom_line(GLA.pub,mapping=aes(x=depth_sample,y=d18O_measurement)) + 
           geom_point(GLA.pub,mapping=aes(x=depth_sample,y=d18O_measurement)))

### Depths (mm)
mid_d18O_depth_0 <- c(37.20,39.20)
mid_d18O_depth_1 <- c(20.20,20.95)
mid_d18O_depth_2 <- c(15.35,15.50) 
mid_d18O_depth_3 <- c(14.45,14.75) 
mid_d18O_depth_4 <- c(12.20,12.80) 
mid_d18O_depth_5 <- c(11.10,11.45) 

mid_d18O_depths <- list(
  cp0 = mid_d18O_depth_0, 
  cp1 = mid_d18O_depth_1,  # ~17.9k
  cp2 = mid_d18O_depth_2,  # ~16.1k
  cp3 = mid_d18O_depth_3,  # ~15.3k
  cp4 = mid_d18O_depth_4,   # ~14.7k
  cp5 = mid_d18O_depth_5   # ~14.7k
)

### Assigned d18O values
mid_d18O_value_0 <- sapply(mid_d18O_depth_0, function(d) {
  GLA.pub$d18O_measurement[
    which.min(abs(GLA.pub$depth_sample - d))
  ]
})
mid_d18O_value_1 <- sapply(mid_d18O_depth_1, function(d) {
  GLA.pub$d18O_measurement[
    which.min(abs(GLA.pub$depth_sample - d))
  ]
})
mid_d18O_value_2 <- sapply(mid_d18O_depth_2, function(d) {
  GLA.pub$d18O_measurement[
    which.min(abs(GLA.pub$depth_sample - d))
  ]
})
mid_d18O_value_3 <- sapply(mid_d18O_depth_3, function(d) {
  GLA.pub$d18O_measurement[
    which.min(abs(GLA.pub$depth_sample - d))
  ]
})
mid_d18O_value_4 <- sapply(mid_d18O_depth_4, function(d) {
  GLA.pub$d18O_measurement[
    which.min(abs(GLA.pub$depth_sample - d))
  ]
})
mid_d18O_value_5 <- sapply(mid_d18O_depth_5, function(d) {
  GLA.pub$d18O_measurement[
    which.min(abs(GLA.pub$depth_sample - d))
  ]
})


mid_d18O_values <- list(
  cp0 = mid_d18O_value_0,  # ~23.8k
  cp1 = mid_d18O_value_1,  # ~17.9k
  cp2 = mid_d18O_value_2,  # ~16.1k
  cp3 = mid_d18O_value_3,  # ~15.3k
  cp4 = mid_d18O_value_4,   # ~14.7k
  cp5 = mid_d18O_value_5   # ~14.7k
)

# --------------------
# Define the Summarising Dataframe for output
# --------------------
summary_df <- data.frame(
  event_med = character(),
  event        = character(),
  event_gra    = character(),
  
  mid_age_author = numeric(),
  mid_d18O_author = numeric(),
  start_age_author = numeric(),
  end_age_author = numeric(),
  duration_author = numeric(),
  slope_author = numeric(),
  
  mid_age_lwr  = numeric(),
  mid_age_med  = numeric(),
  mid_age_upr  = numeric(),
  
  mid_ka_uncert_neg = numeric(),
  mid_ka_med = numeric(),
  mid_ka_uncert_neg = numeric(),
  
  mid_d18O_lwr = numeric(),
  mid_d18O_med = numeric(),
  mid_d18O_upr = numeric(),
  
  start_age_lwr = numeric(),
  start_age_med = numeric(),
  start_age_upr = numeric(),
  
  end_age_lwr = numeric(),
  end_age_med = numeric(),
  end_age_upr = numeric(),
  
  duration_lwr = numeric(),
  duration_med = numeric(),
  duration_upr = numeric(),
  
  slope_lwr = numeric(),
  slope_med = numeric(),
  slope_upr = numeric(),
  
  depth_lwr = numeric(),
  depth_upr = numeric(),
  
  name        = character(),
  
  stringsAsFactors = FALSE
)



# --------------------
# Prediction: Predict ages for respective depth from the age model ensemble
# --------------------

#Load age model data
mid_d18O_ages_0 <- predict(a1,newPositions=mid_d18O_depth_0)
mid_d18O_ages_1 <- predict(a1,newPositions=mid_d18O_depth_1)
mid_d18O_ages_2 <- predict(a1,newPositions=mid_d18O_depth_2)
mid_d18O_ages_3 <- predict(a1,newPositions=mid_d18O_depth_3)
mid_d18O_ages_4 <- predict(a1,newPositions=mid_d18O_depth_4)
mid_d18O_ages_5 <- predict(a1,newPositions=mid_d18O_depth_5)

mid_d18O_ages <- list(
  cp0 = mid_d18O_ages_0,  # ~17.9k
  cp1 = mid_d18O_ages_1,  # ~17.9k
  cp2 = mid_d18O_ages_2,  # ~16.1k
  cp3 = mid_d18O_ages_3,  # ~15.3k
  cp4 = mid_d18O_ages_4,   # ~14.7k
  cp5 = mid_d18O_ages_5   # ~14.7k
)


# -------------------------------------
# LOOP OVER MCMC ENSEMBLE & 4 d18O Midpoints
# -------------------------------------

for (name in names(mid_d18O_depths)) {
  # Storage matrix:
  # columns: slope, intercept, midpoint_x, midpoint_y
  n <- length(mid_d18O_ages[[name]][,1])
  results <- matrix(NA, nrow = n, ncol = 8)
  colnames(results) <- c("beta", "alpha", "mid_x", "mid_y","x1","x2","y1","y2")

  for (i in 1:n) {
    
    # Draw ages from MCMC
    x1 <- mid_d18O_ages[[name]][i,1]
    x2 <- mid_d18O_ages[[name]][i,2]
    
    y1_obs <- mid_d18O_values[[name]][1]
    y2_obs <- mid_d18O_values[[name]][2]
    
    # Add measurement noise #Do i need this?? how can I frame this (maybe it is good indeed)
    y1 <- rnorm(1, y1_obs, sigma_y)
    y2 <- rnorm(1, y2_obs, sigma_y)
    
    x <- c(x1, x2)
    y <- c(y1, y2)
    
    beta  <- (y2 - y1) / (x2 - x1) # slope
    alpha <- y1 - beta * x1 # intercept
    
    mid_x <- (x1 + x2) / 2
    mid_y <- alpha + beta * mid_x
    
    results[i, ] <- c(beta, alpha, mid_x, mid_y,x1,x2,y1,y2)

  }
  
  # Compute also the values for the high percentile age model
  y1 <- y1_obs
  y2 <- y2_obs
  
  x <- sapply(mid_d18O_depths[[name]], function(d) {
    GLA.pub$interp_age[
      which.min(abs(GLA.pub$depth_sample - d))
    ]
  }) 
  x1 <- x[1]
  x2 <- x[2]
  
  beta  <- (y2 - y1) / (x2 - x1) # slope
  alpha <- y1 - beta * x1 # intercept
  
  mid_x <- (x1 + x2) / 2
  mid_y <- alpha + beta * mid_x
  
  author_results <- c(beta, alpha, mid_x, mid_y)
  author_dx <- x2-x1
  list_result[[name]] <- results
  
  # Prepare for Export
  # Quantiles
  mid_xx <- quantile(results[,3], probs = c(0.025, 0.5,0.975), na.rm = TRUE)
  mid_yy <- quantile(results[,4], probs = c(0.025, 0.5,0.975), na.rm = TRUE)
  nomen_m <- sprintf("%.1fk_S", mid_xx[2] / 1000)
  nomen <- sprintf("%.1fk_S", author_results[3] / 1000)
  nomen_gra <- sprintf("%.1fk", author_results[3] / 1000)
  
  end_x <- quantile(mid_d18O_ages[[name]][,1], probs = c(0.025, 0.5,0.975), na.rm = TRUE)
  start_x <- quantile(mid_d18O_ages[[name]][,2], probs = c(0.025, 0.5,0.975), na.rm = TRUE)
  dx <- quantile(mid_d18O_ages[[name]][,2]-mid_d18O_ages[[name]][,1], probs = c(0.025, 0.5,0.975), na.rm = TRUE)
  slope <- quantile(-results[,1], probs = c(0.025, 0.5,0.975), na.rm = TRUE)  

  # --- Create row ---
  new_row <- data.frame(
    event_med = nomen_m,
    event = nomen,
    event_gra = nomen_gra,
    
    mid_age_author = author_results[3],
    mid_d18O_author = author_results[4],
    start_age_author = x2,
    start_d18O_author = y2_obs,
    end_age_author = x1,
    end_d18O_author = y1_obs,
    duration_author = author_dx,
    slope_author = -author_results[1],
    
    mid_age_lwr = mid_xx[1],
    mid_age_med = mid_xx[2],
    mid_age_upr = mid_xx[3],
    
    mid_ka_uncert_neg = (author_results[3]-mid_xx[1])/1000,
    mid_ka_med = mid_xx[2]/1000,
    mid_ka_uncert_pos = (mid_xx[3]-author_results[3])/1000,
    
    mid_d18O_lwr = mid_yy[1],
    mid_d18O_med = mid_yy[2],
    mid_d18O_upr = mid_yy[3],
    
    start_age_lwr = start_x[1],
    start_age_med = start_x[2],
    start_age_upr = start_x[3],
    
    end_age_lwr = end_x[1],
    end_age_med = end_x[2],
    end_age_upr = end_x[3],
    
    duration_lwr = dx[1],
    duration_med = dx[2],
    duration_upr = dx[3],
    
    slope_lwr = slope[1],
    slope_med = slope[2],
    slope_upr = slope[3],
    
    depth_lwr = mid_d18O_depths[[name]][1],
    depth_upr = mid_d18O_depths[[name]][2],
    
    name = name
    
  )
  
  summary_df <- rbind(summary_df, new_row)
  
}


# -------------------------------------
# Formatting, and Preparing the Output
# -------------------------------------

# Formatting
summary_df$mid_age_lwr <- summary_df$mid_age_lwr / 1000
summary_df$mid_age_med <- summary_df$mid_age_med / 1000
summary_df$mid_age_upr <- summary_df$mid_age_upr / 1000
summary_df$mid_age_author <- summary_df$mid_age_author / 1000


summary_df$start_age_lwr <- summary_df$start_age_lwr / 1000
summary_df$start_age_med <- summary_df$start_age_med / 1000
summary_df$start_age_author <- summary_df$start_age_author / 1000
summary_df$start_age_upr <- summary_df$start_age_upr / 1000

summary_df$end_age_lwr <- summary_df$end_age_lwr / 1000
summary_df$end_age_med <- summary_df$end_age_med / 1000
summary_df$end_age_author <- summary_df$end_age_author / 1000
summary_df$end_age_upr <- summary_df$end_age_upr / 1000

format_ci <- function(author, lwr, upr, digits = 2) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"),
    author, lwr, upr
  )
}

formatted_table <- data.frame(
  Event = summary_df$event,
  
  Midpoint_ka = format_ci(summary_df$mid_age_author,
                          summary_df$mid_age_lwr,
                          summary_df$mid_age_upr,
                          digits = 2),
  
  Duration_yr = format_ci(summary_df$duration_author,
                          summary_df$duration_lwr,
                          summary_df$duration_upr,
                          digits = 0),
  
  Slope_permil_per_yr = format_ci(summary_df$slope_author,
                                  summary_df$slope_lwr,
                                  summary_df$slope_upr,
                                  digits = 4),
  Slope_unit = 'd18O',
  Depth = paste(summary_df$depth_lwr,'-',summary_df$depth_upr),
  
  stringsAsFactors = FALSE
)

### Export as Excel and Latex table
# Excel
# Create a new workbookx
wb <- createWorkbook()

#add d18O data
addWorksheet(wb, "GLA d18O")  # Add a sheet
writeData(wb, sheet = "GLA d18O", summary_df)  # Write dataframe to the sheet<

# Save the workbook to an Excel file (dont save for now)
saveWorkbook(wb, "Transition_Analysis_v4.xlsx", overwrite = TRUE)


# Latex Table Output
# Rename columns for LaTeX output
latex_table <- formatted_table %>%
  rename(
    `Midpoint (ka)` = Midpoint_ka,
    `Duration (yr)` = Duration_yr,
    `Slope ($\\permil$ yr$^{-1}$)` = Slope_permil_per_yr,
    `Slope units` = Slope_unit,
    `Sample Depth`= Depth
  )

# Generate LaTeX code
latex_code <- kable(
  latex_table,
  format = "latex",
  booktabs = TRUE,
  align = "lccc",
  caption = "Posterior estimates of transition characteristics. Values are median (2.5--97.5\\% credible interval).",
  escape = FALSE
)

cat(latex_code)

## Generate Markdown
md_table <- kable(
  latex_table,
  format = "markdown",
  align = "lccc"
)

cat(md_table, sep = "\n")


# -------------------------------------
# Further Diagnostic Plot (In Paper Supplementary)
# -------------------------------------

list_result <- list_result[summary_df$name]

cairo_pdf("figure_d18O.pdf", width = 7*1.5, height = 5*1.5, family = "Helvetica")
par(mfrow = c(2,3),
    mar = c(3.2, 4.2, 1.8, 0.8),
    mgp = c(2.1, 0.6, 0),
    cex = 0.85,
    tcl = -0.3
    )
#summary_df$name <- c("cp0","cp1","cp2","cp3","cp4")

j = 0

for (name in names(list_result)) {
  j = j+1
  res <- list_result[[name]]
  
  # Extract columns
  beta  <- res[,1]
  alpha <- res[,2]
  mid_x <- res[,3]
  x1    <- res[,5]
  x2    <- res[,6]
  y1    <- res[,7]
  y2    <- res[,8]
  
  # Median midpoint
  mid_med <- median(mid_x, na.rm = TRUE)
  
  # Define plotting window (±600 years)
  xlim_range <- c(mid_med - 600, mid_med + 600)
  
  # Determine y-range from all posterior endpoints
  ylim_range <- range(c(y1, y2), na.rm = TRUE)
  
  # Create empty plot
  plot(NA,
       xlim = xlim_range,
       ylim = ylim_range,
       xlab = "Age (yr BP)",
       ylab = expression(delta^{18}*O~"\u2030"),
       main = summary_df[j,"event"])
  
  
  
  # Draw posterior lines
  for (i in 1:nrow(res)) {
    
    lines(c(x1[i], x2[i]),
          c(y1[i], y2[i]),
          col = rgb(0,0,1,0.05))  # transparent blue
  }
  
  # Highlight median line
  med_index <- which.min(abs(mid_x - mid_med))
  
  # Add original data
  lines(GLA.pub$interp_age, GLA.pub$d18O_measurement,
        col = "grey")
  points(GLA.pub$interp_age, GLA.pub$d18O_measurement,
         pch = 16, cex = 0.4, col = "black")
  
  df_plot <- summary_df[j,]
  lines(c(df_plot$start_age_author*1000, df_plot$end_age_author*1000),
        c(df_plot$start_d18O_author,df_plot$end_d18O_author),
        col = "yellow",
        lwd = 2)

  
  lines(c(x1[med_index], x2[med_index]),
        c(y1[med_index], y2[med_index]),
        col = "red",
        lwd = 2,lty="dashed")
  
  points(
    df_plot$mid_age_author, df_plot$mid_d18O_author,
         pch = 2, cex = 1, col = "yellow")
  
  arrows(
    x0 = df_plot$mid_age_lwr*1000,
    y0 = df_plot$mid_d18O_author,
    x1 = df_plot$mid_age_upr*1000,
    y1 = df_plot$mid_d18O_author,
    angle = 90,
    code = 3,
    length = 0.05,
    col = "yellow",
    lwd = 1.5
  )
  
  
}
dev.off()


# d13C_init
# --------------------
# Data and Library Loading, Depth definition
# --------------------
### Depths for which the midpoint in between should be computed
### These have been picked manually using plotly

ggplotly(ggplot() + geom_line(GLA.pub,mapping=aes(x=depth_sample,y=d13C_init)) + 
           geom_point(GLA.pub,mapping=aes(x=depth_sample,y=d13C_init)))
mid_d13C_depth_1 <- c(10.7,11.50) # should be ca 17.9k
mid_d13C_depth_2 <- c(17.75,18.30) # should be ca 16.1k
mid_d13C_depth_3 <- c(38.2,39.6) # should be ca 16.1k



mid_d13C_depths <- list(
  cp1 = mid_d13C_depth_1,  
  cp2 = mid_d13C_depth_2,  
  cp3 = mid_d13C_depth_3
)

### d13C values
mid_d13C_value_1 <- sapply(mid_d13C_depth_1, function(d) {
  GLA.pub$d13C_init[
    which.min(abs(GLA.pub$depth_sample - d))
  ]
})
mid_d13C_value_2 <- sapply(mid_d13C_depth_2, function(d) {
  GLA.pub$d13C_init[
    which.min(abs(GLA.pub$depth_sample - d))
  ]
})

mid_d13C_value_3 <- sapply(mid_d13C_depth_3, function(d) {
  GLA.pub$d13C_init[
    which.min(abs(GLA.pub$depth_sample - d))
  ]
})


mid_d13C_values <- list(
  cp1 = mid_d13C_value_1,  # ~14.7k
  cp2 = mid_d13C_value_2,  # ~17.1k
  cp3 = mid_d13C_value_3  # ~15.3k
  #cp4 = mid_d13C_value_4   # ~14.7k
)




# Midpoint 1
# Predict ages from age model ensemble
#Load age model data
mid_d13C_ages_1 <- predict(a1,newPositions=mid_d13C_depth_1)
mid_d13C_ages_2 <- predict(a1,newPositions=mid_d13C_depth_2)
mid_d13C_ages_3 <- predict(a1,newPositions=mid_d13C_depth_3)


mid_d13C_ages <- list(
  cp1 = mid_d13C_ages_1,  # ~17.9k
  cp2 = mid_d13C_ages_2,  # ~16.1k
  cp3 = mid_d13C_ages_3  # ~15.3k
  #cp4 = mid_d13C_ages_4   # ~14.7k
)

# --------------------
# Define the Summarising Dataframe for output
# --------------------
summary_df_d13C <- data.frame(
  event        = character(),
  event_gra    = character(),
  
  mid_age_author = numeric(),
  mid_d13C_author = numeric(),
  start_age_author = numeric(),
  end_age_author = numeric(),
  duration_author = numeric(),
  slope_author = numeric(),
  
  mid_age_lwr  = numeric(),
  mid_age_med  = numeric(),
  mid_age_upr  = numeric(),
  
  mid_ka_uncert_neg = numeric(),
  mid_ka_med = numeric(),
  mid_ka_uncert_neg = numeric(),
  
  mid_d13C_lwr = numeric(),
  mid_d13C_med = numeric(),
  mid_d13C_upr = numeric(),
  
  start_age_lwr = numeric(),
  start_age_med = numeric(),
  start_age_upr = numeric(),
  
  end_age_lwr = numeric(),
  end_age_med = numeric(),
  end_age_upr = numeric(),
  
  duration_lwr = numeric(),
  duration_med = numeric(),
  duration_upr = numeric(),
  
  slope_lwr = numeric(),
  slope_med = numeric(),
  slope_upr = numeric(),
  
  stringsAsFactors = FALSE
)




# -------------------------------------
# LOOP OVER MCMC ENSEMBLE & 4 d13C Breakpoints
# -------------------------------------

for (name in names(mid_d13C_depths)) {
  # Storage matrix:
  # columns: slope, intercept, midpoint_x, midpoint_y
  n <- length(mid_d13C_ages[[name]][,1])
  results <- matrix(NA, nrow = n, ncol = 8)
  colnames(results) <- c("beta", "alpha", "mid_x", "mid_y","x1","x2","y1","y2")
  
  for (i in 1:n) {
    
    # Draw ages from MCMC
    x1 <- mid_d13C_ages[[name]][i,1]
    x2 <- mid_d13C_ages[[name]][i,2]
    
    y1_obs <- mid_d13C_values[[name]][1]
    y2_obs <- mid_d13C_values[[name]][2]
    
    # Add measurement noise #Do i need this?? how can I frame this (maybe it is good indeed)
    y1 <- rnorm(1, y1_obs, sigma_y)
    y2 <- rnorm(1, y2_obs, sigma_y)
    
    x <- c(x1, x2)
    y <- c(y1, y2)
    
    beta  <- (y2 - y1) / (x2 - x1) # slope
    alpha <- y1 - beta * x1 # intercept
    
    mid_x <- (x1 + x2) / 2
    mid_y <- alpha + beta * mid_x
    
    results[i, ] <- c(beta, alpha, mid_x, mid_y,x1,x2,y1,y2)
    
  }
  
  # Compute also the values for the author age model
  y1 <- y1_obs
  y2 <- y2_obs
  
  x <- sapply(mid_d13C_depths[[name]], function(d) {
    GLA.pub$interp_age[
      which.min(abs(GLA.pub$depth_sample - d))
    ]
  }) 
  x1 <- x[1]
  x2 <- x[2]
  
  beta  <- (y2 - y1) / (x2 - x1) # slope
  alpha <- y1 - beta * x1 # intercept
  
  mid_x <- (x1 + x2) / 2
  mid_y <- alpha + beta * mid_x
  
  author_results <- c(beta, alpha, mid_x, mid_y)
  author_dx <- x2-x1
  list_result[[name]] <- results
  
  # Prepare for Export
  # Quantiles
  mid_xx <- quantile(results[,3], probs = c(0.025, 0.5,0.975), na.rm = TRUE)
  mid_yy <- quantile(results[,4], probs = c(0.025, 0.5,0.975), na.rm = TRUE)
  nomen <- sprintf("%.1fk_T", author_results[3] / 1000)
  nomen_gra <- sprintf("%.1fk", author_results[3] / 1000)
  
  end_x <- quantile(mid_d13C_ages[[name]][,1], probs = c(0.025, 0.5,0.975), na.rm = TRUE)
  start_x <- quantile(mid_d13C_ages[[name]][,2], probs = c(0.025, 0.5,0.975), na.rm = TRUE)
  dx <- quantile(mid_d13C_ages[[name]][,2]-mid_d13C_ages[[name]][,1], probs = c(0.025, 0.5,0.975), na.rm = TRUE)
  slope <- quantile(-results[,1], probs = c(0.025, 0.5,0.975), na.rm = TRUE)  
  
  # --- Create row ---
  new_row <- data.frame(
    event = nomen,
    event_gra = nomen_gra,
    
    mid_age_author = author_results[3],
    mid_d13C_author = author_results[4],
    start_age_author = x2,
    start_d13C_author = y2_obs,
    end_age_author = x1,
    end_d13C_author = y1_obs,
    duration_author = author_dx,
    slope_author = -author_results[1],
    
    mid_age_lwr = mid_xx[1],
    mid_age_med = mid_xx[2],
    mid_age_upr = mid_xx[3],
    
    mid_ka_uncert_neg = (author_results[3]-mid_xx[1])/1000,
    mid_ka_med = mid_xx[2]/1000,
    mid_ka_uncert_pos = (mid_xx[3]-author_results[3])/1000,
    
    mid_d13C_lwr = mid_yy[1],
    mid_d13C_med = mid_yy[2],
    mid_d13C_upr = mid_yy[3],
    
    start_age_lwr = start_x[1],
    start_age_med = start_x[2],
    start_age_upr = start_x[3],
    
    end_age_lwr = end_x[1],
    end_age_med = end_x[2],
    end_age_upr = end_x[3],
    
    duration_lwr = dx[1],
    duration_med = dx[2],
    duration_upr = dx[3],
    
    slope_lwr = slope[1],
    slope_med = slope[2],
    slope_upr = slope[3]
  )
  
  summary_df_d13C <- rbind(summary_df_d13C, new_row)
  
}

# -------------------------------------
# Formatting, and Preparing the Output
# -------------------------------------
# Formatting
summary_df_d13C$mid_age_lwr <- summary_df_d13C$mid_age_lwr / 1000
summary_df_d13C$mid_age_med <- summary_df_d13C$mid_age_med / 1000
summary_df_d13C$mid_age_upr <- summary_df_d13C$mid_age_upr / 1000
summary_df_d13C$mid_age_author <- summary_df_d13C$mid_age_author / 1000


summary_df_d13C$start_age_lwr <- summary_df_d13C$start_age_lwr / 1000
summary_df_d13C$start_age_med <- summary_df_d13C$start_age_med / 1000
summary_df_d13C$start_age_author <- summary_df_d13C$start_age_author / 1000
summary_df_d13C$start_age_upr <- summary_df_d13C$start_age_upr / 1000

summary_df_d13C$end_age_lwr <- summary_df_d13C$end_age_lwr / 1000
summary_df_d13C$end_age_med <- summary_df_d13C$end_age_med / 1000
summary_df_d13C$end_age_author <- summary_df_d13C$end_age_author / 1000
summary_df_d13C$end_age_upr <- summary_df_d13C$end_age_upr / 1000

format_ci <- function(author, lwr, upr, digits = 2) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"),
    author, lwr, upr
  )
}

formatted_table <- data.frame(
  Event = summary_df_d13C$event,
  
  Midpoint_ka = format_ci(summary_df_d13C$mid_age_author,
                          summary_df_d13C$mid_age_lwr,
                          summary_df_d13C$mid_age_upr,
                          digits = 2),
  
  Duration_yr = format_ci(summary_df_d13C$duration_author,
                          summary_df_d13C$duration_lwr,
                          summary_df_d13C$duration_upr,
                          digits = 0),
  
  Slope_permil_per_yr = format_ci(summary_df_d13C$slope_author,
                                  summary_df_d13C$slope_lwr,
                                  summary_df_d13C$slope_upr,
                                  digits = 4),
  Slope_unit = 'd13C',
  
  stringsAsFactors = FALSE
)

### Export as excel and latex table
# Excel
# Create a new workbookx
wb <- createWorkbook()

#add d13C data
addWorksheet(wb, "GLA d13C")  # Add a sheet
writeData(wb, sheet = "GLA d13C", summary_df_d13C)  # Write dataframe to the sheet<

# Save the workbook to an Excel file
saveWorkbook(wb, "Transition_Analysis_v4_d13C.xlsx", overwrite = TRUE)


library(dplyr)

# Rename columns for LaTeX output
latex_table <- formatted_table %>%
  rename(
    `Midpoint (ka)` = Midpoint_ka,
    `Duration (yr)` = Duration_yr,
    `Slope ($\\permil$ yr$^{-1}$)` = Slope_permil_per_yr,
    `Slope units` = Slope_unit
  )

# Generate LaTeX code
library(knitr)

latex_code <- kable(
  latex_table,
  format = "latex",
  booktabs = TRUE,
  align = "lccc",
  caption = "Posterior estimates of transition characteristics. Values are median (2.5--97.5\\% credible interval).",
  escape = FALSE
)

cat(latex_code)


# -------------------------------------
# Further Diagnostic Plot (In Paper Supplementary)
# -------------------------------------

cairo_pdf("figure_d13C.pdf", width = 7*1.5, height = 2.5*1.5, family = "Helvetica")
par(mfrow = c(1,3),
    mar = c(3.2, 4.2, 1.8, 0.8),
    mgp = c(2.1, 0.6, 0),
    cex = 0.85,
    tcl = -0.3
)
summary_df_d13C$name <- c("cp1","cp2","cp3")

j = 0

for (name in summary_df_d13C$name) {
  j = j+1
  res <- list_result[[name]]
  
  # Extract columns
  beta  <- res[,1]
  alpha <- res[,2]
  mid_x <- res[,3]
  x1    <- res[,5]
  x2    <- res[,6]
  y1    <- res[,7]
  y2    <- res[,8]
  
  # Median midpoint
  mid_med <- median(mid_x, na.rm = TRUE)
  
  # Define plotting window (±600 years)
  xlim_range <- c(mid_med - 600, mid_med + 600)
  
  # Determine y-range from all posterior endpoints
  ylim_range <- range(c(y1, y2), na.rm = TRUE)
  
  # Create empty plot
  plot(NA,
       xlim = xlim_range,
       ylim = ylim_range,
       xlab = "Age (yr BP)",
       ylab = expression(delta^{13}*C~"\u2030"),
       main = summary_df_d13C[j,"event"])
  
  
  
  # Draw posterior lines
  for (i in 1:nrow(res)) {
    
    lines(c(x1[i], x2[i]),
          c(y1[i], y2[i]),
          col = rgb(0,0,1,0.05))  # transparent blue
  }
  
  # Highlight median line
  med_index <- which.min(abs(mid_x - mid_med))
  
  # Add original data
  lines(GLA.pub$interp_age, GLA.pub$d13C_init,
        col = "grey")
  points(GLA.pub$interp_age, GLA.pub$d13C_init,
         pch = 16, cex = 0.4, col = "black")
  
  df_plot <- summary_df_d13C[j,]
  lines(c(df_plot$start_age_author*1000, df_plot$end_age_author*1000),
        c(df_plot$start_d13C_author,df_plot$end_d13C_author),
        col = "yellow",
        lwd = 2)
  
  
  lines(c(x1[med_index], x2[med_index]),
        c(y1[med_index], y2[med_index]),
        col = "red",
        lwd = 2,lty="dashed")
  
  points(
    df_plot$mid_age_author, df_plot$mid_d13C_author,
    pch = 2, cex = 1, col = "yellow")
  
  arrows(
    x0 = df_plot$mid_age_lwr*1000,
    y0 = df_plot$mid_d13C_author,
    x1 = df_plot$mid_age_upr*1000,
    y1 = df_plot$mid_d13C_author,
    angle = 90,
    code = 3,
    length = 0.05,
    col = "yellow",
    lwd = 1.5
  )
  
  
}

dev.off()



# --------------------
# 3) temporal relationship between indicators
# --------------------

## Duration between Freshening at 17.8k and cooling at 17.0k.

# Author Age Model
author_duration_1 = (summary_df$mid_age_author[2]- summary_df_d13C$mid_age_author[2])*1000 

# Compute CI
mid_duration_1 = mid_d18O_ages_1 - mid_d13C_ages_2
ci_duration_1 = quantile(mid_duration_1, probs = c(0.025, 0.5,0.975), na.rm = TRUE)

# In Main Text
round(author_duration_1) # Duration 
round(ci_duration_1 - author_duration_1) # Duration CI



