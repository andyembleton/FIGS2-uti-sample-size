# Andy Embleton
# Jan-2025
# FIGS2 trial design

library(simstudy)
library(MASS)
library(mvtnorm)
library(ggplot2)

# Setup
set.seed(20250203)
n.sims <- 5000
n <- 50
mu0 <- 2.5
diff <- 1
dropout <- 0.2

correlations <- seq(0.3, 0.7, by = 0.05)
thetas <- seq(3, 5, by = 1) # For reference: variance <- mu0 + ((mu0^2) / theta)

target.diff <- (mu0 - diff) / mu0
mu1 <- mu0 * target.diff

results <- expand.grid(correlation = correlations, theta = thetas, power = NA)

progress = 0
progressmax = n.sims * length(correlations) * length(thetas)
progressbar = txtProgressBar(min = 0, max = progressmax, initial = 1, style = 3) 

for (c in correlations) {
  for (t in thetas) {
    outcomes <- data.frame(pvalues = rep(NA, n.sims))
    for (i in 1:n.sims) {
      progress = progress + 1
      
      # Generate correlated normal data using compound symmetry
      sigma <- matrix(c(1, c, c, 1), nrow = 2)
      normal_data_control <- rmvnorm(n, mean = c(0, 0), sigma = sigma)
      normal_data_active <- rmvnorm(n, mean = c(0, 0), sigma = sigma)
      
      # Transform normal data to uniform
      uniform_data_control <- pnorm(normal_data_control)
      uniform_data_active <- pnorm(normal_data_active)
      
      # Transform uniform data to negative binomial
      baseline_control <- qnbinom(uniform_data_control[, 1], size = t, mu = mu0)
      follow_up_control <- qnbinom(uniform_data_control[, 2], size = t, mu = mu0)
      baseline_active <- qnbinom(uniform_data_active[, 1], size = t, mu = mu0)
      follow_up_active <- qnbinom(uniform_data_active[, 2], size = t, mu = mu1)
      
      df.data.control <- data.frame(id = 1:n, month0 = baseline_control, month6 = follow_up_control, treatment = "control")
      df.data.active <- data.frame(id = (n + 1):(2 * n), month0 = baseline_active, month6 = follow_up_active, treatment = "active")
      
      df.data.both <- rbind(df.data.control, df.data.active)
      
      # Apply random dropout
      dropout_indices <- sample(1:nrow(df.data.both), size = floor(dropout * nrow(df.data.both)))
      df.data.both <- df.data.both[-dropout_indices, ]
      
      fit <- glm(month6 ~ treatment + month0, data = df.data.both, family = 'poisson')
      # negative binomial [glm.nb] does not seem to perform as well

      outcomes$pvalues[i] <- summary(fit)$coefficients[2, 4]
      
      setTxtProgressBar(progressbar, progress)
    }
    power <- length(outcomes$pvalues[outcomes$pvalues < 0.05]) / n.sims
    results[results$correlation == c & results$theta == t, "power"] <- power
  }
}
close(progressbar)

# Plot
ggplot(results, aes(x = correlation, y = power, color = as.factor(theta))) +
  geom_line() +
  scale_x_continuous(breaks = correlations) +
  scale_y_continuous(limits = c(0.75, 0.95), breaks = seq(0.75, 0.95, by = 0.05)) +
  labs(title = bquote("Simulated power ("*mu[0]*"="*.(mu0)*", "*mu[1]*"="*.(mu1)*")") ,
       x = "Correlation",
       y = "Power",
       caption = paste0(n.sims, " simulations, effect size=", target.diff, ", random ", dropout*100, "% dropout, correlation ", min(correlations), " to ", max(correlations), ", theta ", max(thetas), " to ", min(thetas), " (equivalent to variances of ", round(mu0 + ((mu0^2) / max(thetas)),1), " to ", round(mu0 + ((mu0^2) / min(thetas)),1), ")") ,
       color = "Theta") +
  guides(color = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0.85, linetype = "dashed") +
  theme_minimal()