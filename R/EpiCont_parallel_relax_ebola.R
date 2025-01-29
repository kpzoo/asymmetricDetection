library(VGAM)
#library(future.apply)
#library(foreach)
library(parallel)
library(pbapply)
library(zoo) # for the rollsum function
library(EpiControl)

cores=detectCores()-1
cl <- makeCluster(cores)

Epi_pars = data.frame (
  Pathogen = c("COVID-19", "Ebola"),
  R0 = c(2.5, 1.9),
  gen_time = c(6.5, 15.0),
  gen_time_var = c(2.1, 2.1),
  CFR = c(0.0132, 0.5),
  mortality_mean = c(10.0, 10.0),
  mortality_var = c(1.1, 1.1)
)

r_trans_steep <- 1.5  # Growth rate
r_trans_len <- 7  # Number of days for the transition
t0 <- r_trans_len  / 2 # Midpoint of the transition

logistic_function <- function(t, R0, R1, r, t0) {
  K <- R1
  L <- R0
  (K-L) / (1 + exp(-r_trans_steep * (t - t0))) + L
}

I0 <- 10L #initial no. of infections
ndays <- 40L*7L#epidemic length
first_stint <-19L*7L
N <- 1e7 #Total population (if we account for susceptibles)

# for ebola
Noise_pars <- data.frame (
  mtau = 22.34, #Reporting delay mean
  r = 16.0, #Reporting delay variance
  ur_mean = 0.4, #Under reporting, mean/variance
  ur_beta_b = 20.0
)

# for ebola
Action_space <- data.frame (
  NPI = c("No restrictions", "Lockdown"),
  R_coeff = c(1.0, 0.5/1.9), #R0_act = R0 * ctrl_states
  R_beta_a = c(0.0, 5.0), #R0_act uncertainty
  cost_of_NPI = c(0.0, 0.15)
)

Action_space0 <- data.frame (
  NPI = c("No restrictions"),
  R_coeff = c(1.0), #R0_act = R0 * ctrl_states
  R_beta_a = c(0.0), #R0_act uncertainty
  cost_of_NPI = c(0.0)
)

#Sim and control options
cost_sel <- 1L # : bilinear+oveshoot, 2: flat+quadratic for overshoot, 3: ReLu + Overshoot
use_inc <- 1L #1: control for incidence, 0: control for infectiousness
delay_calc_v <- 0L #1: as in ref, 0: from incidence
under_rep_calc <- 3L #1: as in ref, 0: separately, from incidence, 2: same, but using a beta distribution for rho
distr_sel <- 1L #0: Deterministic, 1: Poisson, 2: Binomial
delay_on <- 1L #1: sim with time-delay, 0: no-delay (if 0, set delay_calc_v = 0)
under_rep_on <- 1L #0: no under reporting, 1: calculate with under-reporting

C_target <- 50 #target cases
C_target_pen <- C_target*1.5 #overshoot penalty threshold
R_target <- 1.0
D_target <- 12
D_target_pen <- 50 #max death
#alpha <- 1.3/C_target #~proportional gain (regulates error in cases) covid
alpha = 3.25/C_target #~proportional gain (regulates error in cases) ebola
#beta <- 0.0 #~derivative gain (regulates error in R)
alpha_d <- 0*1.3/D_target
ovp <- 0*5.0 #overshoot penalty
dovp <- 0*10.0 #death overshoot penalty
gamma <- 0.95 #discounting factor

#Simulation parameters
n_ens <- 100L #MC assembly size for 4
sim_ens <- 1000L #assembly size for full simulation

#Frequecy of policy review
rf <- 7L #days 14
R_est_wind <- 5L #rf-2 #window for R estimation
use_S <- 14L

#Prediction window
pred_days <- 14L #12L #14 #21 #12

# Original episim_data
column_names <- c("days", "sim_id", "I", "Lambda", "C", "Lambda_C", "S", "Deaths", "Re", "Rew", "Rest", "R0est", "policy", "R_coeff")

# Create an empty data frame with specified column names
empty_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(empty_df) <- column_names

zero_matrix <- matrix(0, nrow = ndays, ncol = length(column_names))
colnames(zero_matrix) <- column_names

# Combine empty data frame and zero matrix using rbind
episim_data <- rbind(empty_df, zero_matrix)

#initialisation

episim_data['policy'] <- rep(1, ndays)
episim_data['sim_id'] <- rep(1, ndays)
episim_data['days'] <- 1:ndays
episim_data[1,] <- c(1, 1, I0, I0, Noise_pars['ur_mean']*I0, Noise_pars['ur_mean']*I0, N-I0, 0, Epi_pars[2,'R0'], Epi_pars[2,'R0'], 1, 1, 1)

episim_data_ens <- replicate(sim_ens, episim_data, simplify = FALSE)

pb <- txtProgressBar(min = 0, max = sim_ens, initial = 0, style = 3, width = 50, char = "=")
for (ii in 1:sim_ens) {
  episim_data_ens[[ii]]$sim_id <- rep(ii, ndays)
}

# Initialize the progress bar
pb <- txtProgressBar(min = 0, max = sim_ens, initial = 0, style = 3, width = 50, char = "=")

# Ensure each simulation has a unique sim_id
for (ii in 1:sim_ens) {
  episim_data_ens[[ii]]$sim_id <- rep(ii, ndays)
}

episim_data_ens[[1]] <- Epi_MPC_run_wd_K(episim_data_ens[[1]], Epi_pars, Noise_pars, Action_space0, pred_days = pred_days, start_day = 1, n_ens = n_ens, ndays = 146L, R_est_wind = R_est_wind, pathogen = 2, susceptibles = 0, delay = 1, ur = 1, r_dir = 1, N = N)

for (jj in 1:sim_ens) {
  episim_data_ens[[jj]] <- episim_data_ens[[1]]
  episim_data_ens[[jj]]$sim_id <- rep(jj, ndays)
}

#save(episim_data_ens, file = "episim_data_ens_ebola.RData")

load("episim_data_ens_ebola.RData")

#for (jj in 1:100) {
#  episim_data_ens[[jj]]["C"] <- episim_data_ens[[jj]]["I"]
#}

for (jj in 101:sim_ens) {
  episim_data_ens[[jj]] <- episim_data_ens[[1]]
  episim_data_ens[[jj]]$sim_id <- rep(jj, ndays)
}

clusterExport(cl, ls())

clusterEvalQ(cl, episim_data_ens)
clusterEvalQ(cl, {
  library(VGAM)
})

results <- pblapply(1:sim_ens, function(jj) {
  episim_data_ens[[jj]] <- Epi_MPC_run_wd_K(episim_data_ens[[jj]], Epi_pars, Noise_pars, Action_space, pred_days = pred_days, start_day = 145L, n_ens = n_ens, ndays = nrow(episim_data), R_est_wind = R_est_wind, pathogen = 2, susceptibles = 0, delay = 1, ur = 1, r_dir = 1, N = N)
}, cl=cl)

episim_data_ens <- results
stopCluster(cl)

for (jj in 1:sim_ens) {
  episim_data_ens[[jj]] <- head(episim_data_ens[[jj]], -pred_days)
  episim_data_ens[[jj]]["D_roll"] <- rollsum(episim_data_ens[[jj]]["Deaths"], 7, fill = NA)
  episim_data_ens[[jj]]["I_roll"] <- rollsum(episim_data_ens[[jj]]["I"], 7, fill = NA)
  episim_data_ens[[jj]]["D_cum"] <- cumsum(episim_data_ens[[jj]]["Deaths"])
  episim_data_ens[[jj]]["I_cum"] <- cumsum(episim_data_ens[[jj]]["I"])
}

combined_data <- do.call(rbind, episim_data_ens)

library(ggplot2)

# Create a new grouping variable for discontinuities
combined_data <- combined_data %>%
  group_by(sim_id, policy) %>%
  arrange(days) %>%
  mutate(group = cumsum(c(1, diff(days) != 1))) %>%
  ungroup()

combined_data <- combined_data %>%
  arrange(sim_id, days, policy)

policy_labels <- c("1" = "No intervention", "2" = "Lockdown")


# Initialize lists to collect the results
change_1_to_2_list <- list()
change_2_to_1_list <- list()
difference_list <- list()
ratio_list <- list()

# Loop through 'ii' from 1 to 'sim_ens'
for (ii in 1:sim_ens) {
  # Extract the policy vector
  policy_vector <- (episim_data_ens[[ii]]["policy"])$policy

  # Create a Boolean vector where TRUE represents 'policy == 2'
  bool_vector <- policy_vector == 2

  # Count the number of TRUEs (1s)
  true_count <- sum(bool_vector)

  # Count the number of FALSEs (0s)
  false_count <- length(bool_vector) - true_count

  # Calculate the ratio of TRUEs to FALSEs
  ratio <- true_count / (true_count + false_count)

  # Append the ratio to the ratio list
  ratio_list[[ii]] <- ratio

  # Find the index where it first changes from 1 to 2
  change_1_to_2 <- which(policy_vector[-1] == 2 & policy_vector[-length(policy_vector)] == 1)[1] + 1

  # Find the index where it first changes from 2 to 1
  change_2_to_1 <- which(policy_vector[-1] == 1 & policy_vector[-length(policy_vector)] == 2)[1] + 1

  # Calculate the difference
  difference <- change_2_to_1 - change_1_to_2

  # Append the results to the lists
  change_1_to_2_list[[ii]] <- change_1_to_2
  change_2_to_1_list[[ii]] <- change_2_to_1
  difference_list[[ii]] <- difference
}

# Convert lists to numeric vectors
changes <- unlist(change_1_to_2_list)/rf
changes2 <- unlist(change_2_to_1_list)/rf
differences <- unlist(difference_list)/rf

# Load required libraries
library(gridExtra)

# Simulating sample data
# Left panel plot (with continuous lines and custom labels)
p1 <- ggplot(combined_data %>% filter(sim_id == 1)) +
  geom_step(data = subset(combined_data, sim_id != 1), aes(x = days, y = C, color = as.factor(sim_id)), alpha = 0.1) +
  geom_step(aes(x = days, y = C, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0) +
  geom_step(aes(x = days, y = I, color = factor(policy, labels = policy_labels), group = 1), alpha = 1.0, size=0.25) +
  geom_hline(yintercept = C_target, linetype = "dashed", color = "blue", size = 0.25) +
  labs(x = "Days", y = "Reported cases and true infections", color = "Policy") +
  scale_color_manual(values = c("No intervention" = "chartreuse3", "Lockdown" = "red")) +
  guides(color = guide_legend(title = "Policy")) +
  theme(legend.position = c(0.05, 0.95),   # Position the legend in the top left corner
        legend.justification = c(0, 1))

# Right panel histograms
hist1 <- ggplot(data.frame(changes), aes(x = changes)) +
  geom_histogram(aes(y = (..count..) / sum(..count..) * 100),
                 breaks = seq(12 - 0.5, 22 + 0.5, by = 1),
                 fill = "blue", color = "black") +
  ggtitle("LD Start") +
  xlab("LD start in weeks") +
  ylab("Percentage (%)")

hist2 <- ggplot(data.frame(changes2), aes(x = changes2)) +
  geom_histogram(aes(y = (..count..) / sum(..count..) * 100),
                 breaks = seq(14 - 0.5, 32 + 0.5, by = 1),
                 fill = "blue", color = "black") +
  ggtitle("LD End") +
  xlab("LD ending in weeks") +
  ylab("Percentage (%)")

hist3 <- ggplot(data.frame(differences), aes(x = differences)) +
  geom_histogram(aes(y = (..count..) / sum(..count..) * 100),
                 breaks = seq(1 - 0.5, 10 + 0.5, by = 1),
                 fill = "blue", color = "black") +
  ggtitle("First LD Length") +
  xlab("First LD length in weeks") +
  ylab("Percentage (%)")

# Combine histograms into a grid
hist_plots <- grid.arrange(hist1, hist2, hist3, ncol = 1)

# Arrange left and right panels side by side with different widths
grid.arrange(p1, hist_plots, ncol = 2, widths = c(2, 1))  # Left panel is twice as wide as the right one

#save(episim_data_ens, file = "episim_data_ens_relax_no_noise_ebola_1000.RData")

