# Title -------------------------------------------------------------------
# Name: Modelling for Pandemic Preparedness and Responses
# Organization: GWAC & KNUST
# Topic: Mathematical modelling for Mpox disease in West Africa
# Author: Group 1
# Affiliation: German West African Centre for Global Health and Pandemic Prevention & Kwame Nkurmah University of Science and Technology
# Date: September 2024



# Load the necessary library ----------------------------------------------

pacman::p_load(deSolve, 
               rio, 
               janitor, 
               gtsummary, 
               tidyverse)

mpox <- import("mpox_2024.csv")
# Solve the differential equation and plot the graph ----------------------


# 1. Define the SIR Model function
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Differential equations
    # dS/dt = - (β1Is + β2Ir)S/N								    	(1)
    # dIs/dt = β1IsS/N – (λ1 + µ1)Is   								(2)
    # dIr/dt = β2IrS/N – (λ2 + µ2)Ir									(3)
    # dD/dt = µ1Is + µ2Ir                             (4)
    # dR/dt = λ1Is + λ2Ir										          (5)
    
    
    # Rate of change in the susceptible population(1)
    dS <- -(beta_1*I_s + beta_2 * I_r) * S/N
    
    
    # Rate of change from infectious_sexual contact to recovered (2)
    dIs <- (beta_1 * I_s * S)/N - (lambda_1 + mu_1) * I_s
    
    # Rate of change from infectious_direct contact to recovered  (3)
    dIr <- (beta_2 * I_r * S)/N - (lambda_2 + mu_2) * I_r
    
    # Rate of change from infectious to dead                      (4)
    dD <- mu_1 * I_s + mu_2 * I_r
    
    # Rate of recovery  (5)
    dR <- lambda_1 * I_s + lambda_2 * I_r
    return(list(c(dS, dIs, dIr, dD, dR)))
    
  })
}

# Initial state values (S = Susceptible, I = Infectious, R = Recovered)
init_state <- c(S = 99000000,    # Susceptible human population
                I_s = 4467,        # Population who became infectious through sexual contact
                I_r = 1000,        # Population who became infectious through direct contact 
                D = 0,           # Population who die from the infection
                
                R = 0)         # Population  who recovered

                

# Parameters 
parameters <- c(beta_1 = 1.2,   # transmission rate through sexual contact
                beta_2 = 0.7,   # transmission rate through direct contact
                lambda_1 = 0.00247,     # rate of recovery (due to sexual contact)
                lambda_2 = 0.00247,    # rate of recovery (due to direct contact)
                mu_1 = 0.00028,    # death rate from infection due to sexual contact
                mu_2 = 0.00028,     # death rate from infection due to direct contact
                N = 99005467) # Total population

# Time over which to simulate
time <- seq(from = 0,    # from 0 up to 35 at intervals of 1
            to = 252, 
            by = 1)

# Solve the model using ode function from deSolve package
output <- ode(y = init_state,      # inital state variables
              times = time,        # times sequence
              func = sir_model,    # sir_model function as above
              parms = parameters)  # stipulate parameters as above

# Convert output to a dat frame for easier handling
output_df <- as.data.frame(output)

# Export data to Excel
# export(output_df, "mpox_model.xlsx")

# Data visualization -------------------------------------------------------
# Plot the results

# 1. Comparison of model compartments
output_df %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = S, color = "Susceptible"), size = 1.2)+
  geom_line(aes(y = I_s, color = "Infected due to sexual contact"), size = 1.2)+
  geom_line(aes(y = I_r, color = "Infected due to direct contact"), size = 1.2)+
  geom_line(aes(y = D, color = "Dead"), size = 1.2)+
  geom_line(aes(y = R, color = "Recovered"), size = 1.2)+
  labs(x = "Time",
       y = "Population"
       #title = "Comparison of model compartments"
  )+
  scale_color_manual(
    name = "Population", 
    values = c("Susceptible" = "red",
               "Infected due to sexual contact" = "orange",
               "Infected due to direct contact" = "blue", 
               "Dead" = "purple",
               "Recovered" = "green"))+
  
  #theme_classic()+
  theme(legend.position = c(0.8, 0.8))



# Test area ---------------------------------------------------------------

mpox %>% 
  select(Case_status, Location_Admin0) %>% 
  filter(Location_Admin0 == "Democratic Republic of the Congo") %>% 
  tbl_summary(by = Case_status,
              percent = "row") %>% 
  add_overall() %>% 
  bold_labels()

