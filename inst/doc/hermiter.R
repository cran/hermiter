## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(10)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(hermiter)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(magrittr)
library(ggplot2)
library(dplyr)
library(data.table)
library(DT)

## -----------------------------------------------------------------------------
hermite_est <- hermite_estimator(N=10, standardize=TRUE)

## -----------------------------------------------------------------------------
observations <- rlogis(n=1000)
hermite_est <- hermite_estimator(N=10, standardize=TRUE)
hermite_est <- update_batch(hermite_est,observations)

## -----------------------------------------------------------------------------
observations <- rlogis(n=1000)
hermite_est <- hermite_estimator(N=10, standardize=TRUE)
hermite_est <- hermite_est %>% update_batch(observations)

## -----------------------------------------------------------------------------
observations <- rlogis(n=1000)
hermite_est <- hermite_estimator(N=10, standardize=TRUE)
for (idx in c(1:length(observations))) {
  hermite_est <- update_sequential(hermite_est,observations[idx])
}

## -----------------------------------------------------------------------------
observations <- rlogis(n=1000)
hermite_est <- hermite_estimator(N=10, standardize=TRUE)
for (idx in c(1:length(observations))) {
  hermite_est <- hermite_est %>% update_sequential(observations[idx])
}

## -----------------------------------------------------------------------------
x <- seq(-15,15,0.1)
pdf_est <- dens(hermite_est,x)
cdf_est <- cum_prob(hermite_est,x)

p <- seq(0.05,1,0.05)
quantile_est <- quant(hermite_est,p)

## -----------------------------------------------------------------------------
x <- seq(-15,15,0.1)
pdf_est <- hermite_est %>% dens(x)
cdf_est <- hermite_est %>% cum_prob(x)

p <- seq(0.05,0.95,0.05)
quantile_est <- hermite_est %>% quant(p)

## -----------------------------------------------------------------------------
actual_pdf <- dlogis(x)
actual_cdf <- plogis(x)
df_pdf_cdf <- data.frame(x,pdf_est,cdf_est,actual_pdf,actual_cdf)

actual_quantiles <- qlogis(p)
df_quant <- data.frame(p,quantile_est,actual_quantiles)

## -----------------------------------------------------------------------------
ggplot(df_pdf_cdf,aes(x=x)) + geom_line(aes(y=pdf_est, colour="Estimated")) +
  geom_line(aes(y=actual_pdf, colour="Actual")) +
  scale_colour_manual("", 
                      breaks = c("Estimated", "Actual"),
                      values = c("blue", "black")) + ylab("Probability Density")

## -----------------------------------------------------------------------------
ggplot(df_pdf_cdf,aes(x=x)) + geom_line(aes(y=cdf_est, colour="Estimated")) +
  geom_line(aes(y=actual_cdf, colour="Actual")) +
  scale_colour_manual("", 
                      breaks = c("Estimated", "Actual"),
                      values = c("blue", "black")) +
  ylab("Cumulative Probability")

## -----------------------------------------------------------------------------
ggplot(df_quant,aes(x=actual_quantiles)) + geom_point(aes(y=quantile_est),
                                                      color="blue") +
  geom_abline(slope=1,intercept = 0) +xlab("Theoretical Quantiles") +
  ylab("Estimated Quantiles")

## -----------------------------------------------------------------------------
# Prepare Test Data
test_data <- data.frame()
for (i in c(1:5)) {
  exponential_data <- rexp(n=1000)
  logistic_data <- rlogis(n=1000)
  logn_data <- rlnorm(n=1000)
  test_data <- rbind(test_data,data.frame(dist_name=rep("exp",
                length(exponential_data)),idx=i,observations=exponential_data))
  test_data <- rbind(test_data,data.frame(dist_name=rep("logis",
                      length(logistic_data)),idx=i,observations=logistic_data))
  test_data <- rbind(test_data,data.frame(dist_name=rep("lnorm",
                              length(logn_data)),idx=i,observations=logn_data))
}
setDT(test_data)

## -----------------------------------------------------------------------------
# Group observations by distribution and idx and create Hermite estimators
estimates <- test_data[,.(herm_est = list(hermite_estimator(N=10,
             standardize = TRUE) %>% update_batch(observations))),
             by=.(dist_name,idx)]
estimates

## -----------------------------------------------------------------------------
# Group observations by distribution and combine Hermite estimators
combined_estimates <- estimates[,.(herm_comb = list(combine_hermite(herm_est))),
                                by=.(dist_name)]
combined_estimates

## -----------------------------------------------------------------------------
# Estimate probability densities, cumulative probabilities and quantiles
dens_vals <- combined_estimates[,.(dens_est = list(dens(herm_comb[[1]],
                                            c(0.5,1,1.5,2)))),by=.(dist_name)]
cum_prob_vals <- combined_estimates[,.(cum_prob_est = list(cum_prob(herm_comb[[1]],c(0.5,1,1.5,2)))),by=.(dist_name)]
quantile_vals <- combined_estimates[,.(quantile_est = list(quant(herm_comb[[1]],c(0.25,0.5,0.75)))),by=.(dist_name)]

## -----------------------------------------------------------------------------
# Prepare Test Data
test_data <- data.frame()
for (i in c(1:5)) {
  exponential_data <- rexp(n=1000)
  logistic_data <- rlogis(n=1000)
  logn_data <- rlnorm(n=1000)
  test_data <- rbind(test_data,data.frame(dist_name=rep("exp",
                length(exponential_data)),idx=i,observations=exponential_data))
  test_data <- rbind(test_data,data.frame(dist_name=rep("logis",
                      length(logistic_data)),idx=i,observations=logistic_data))
  test_data <- rbind(test_data,data.frame(dist_name=rep("lnorm",
                              length(logn_data)),idx=i,observations=logn_data))
}

## -----------------------------------------------------------------------------
# Group observations by distribution and idx and create Hermite estimators
estimates <- test_data %>% group_by(dist_name,idx) %>% summarise(herm_est = list(hermite_estimator(N=10,standardize = TRUE) %>% update_batch(observations)))
estimates

## -----------------------------------------------------------------------------
# Group observations by distribution and combine Hermite estimators
combined_estimates <- estimates %>% group_by(dist_name) %>% summarise(herm_comb
                                              = list(combine_hermite(herm_est)))
combined_estimates

## -----------------------------------------------------------------------------
# Estimate probability densities, cumulative probabilities and quantiles
dens_vals <- combined_estimates %>%
  rowwise() %>% mutate(dens_est = list(dens(herm_comb,c(0.5,1,1.5,2))))
cum_prob_vals <- combined_estimates %>%
  rowwise() %>% mutate(cum_prob_est = list(cum_prob(herm_comb,c(0.5,1,1.5,2))))
quantile_vals <- combined_estimates %>%
  rowwise() %>% mutate(quantile_est = list(quant(herm_comb,c(0.25,0.5,0.75))))

## -----------------------------------------------------------------------------
# Compute Mean Absolute Error
dens_vals <- dens_vals %>%
  rowwise() %>% mutate(dens_actual = list(do.call(paste0("d",dist_name),
                list(c(0.5,1,1.5,2))))) %>% mutate(mean_abs_error_density =
                                              mean(abs(dens_est-dens_actual)))
cum_prob_vals <- cum_prob_vals %>%
  rowwise() %>% mutate(cum_prob_actual = list(do.call(paste0("p",dist_name),
                list(c(0.5,1,1.5,2)))))%>% mutate(mean_abs_error_cum_prob = mean(abs(cum_prob_est-cum_prob_actual)))
quantile_vals <- quantile_vals %>%
  rowwise() %>% mutate(quantile_actual= list(do.call(paste0("q",dist_name),
            list(c(0.25,0.5,0.75)))))%>% mutate(mean_abs_error_quantiles = mean(abs(quantile_est-quantile_actual)))
mean_abs_error_summary <- data.frame(dist_name=dens_vals$dist_name, mean_abs_error_density=dens_vals$mean_abs_error_density, mean_abs_error_cum_prob=cum_prob_vals$mean_abs_error_cum_prob,
              mean_abs_error_quantiles=quantile_vals$mean_abs_error_quantiles)

## -----------------------------------------------------------------------------
datatable(mean_abs_error_summary) %>% formatRound(columns =c("mean_abs_error_density","mean_abs_error_cum_prob",
                                        "mean_abs_error_quantiles"),digits = 3)

## ----eval=FALSE---------------------------------------------------------------
#  # Not Run. Copy and paste into app.R and run.
#  library(shiny)
#  library(hermiter)
#  library(ggplot2)
#  library(magrittr)
#  
#  ui <- fluidPage(
#      titlePanel("Streaming Statistics Analysis Example: Exponential
#                 i.i.d. stream"),
#      sidebarLayout(
#          sidebarPanel(
#              sliderInput("percentile", "Percentile:",
#                          min = 0.01, max = 0.99,
#                          value = 0.5, step = 0.01)
#          ),
#          mainPanel(
#             plotOutput("plot"),
#             textOutput("quantile_text")
#          )
#      )
#  )
#  
#  server <- function(input, output) {
#      values <- reactiveValues(hermite_est =
#                                   hermite_estimator(N = 10, standardize = TRUE))
#      x <- seq(-15, 15, 0.1)
#      # Note that the stub below could be replaced with code that reads streaming
#      # data from various sources, Kafka etc.
#      read_stream_stub_micro_batch <- reactive({
#          invalidateLater(1000)
#          new_observation <- rexp(10)
#          return(new_observation)
#      })
#      updated_cdf_calc <- reactive({
#          micro_batch <- read_stream_stub_micro_batch()
#          for (idx in seq_along(micro_batch)) {
#              values[["hermite_est"]] <- isolate(values[["hermite_est"]]) %>%
#                  update_sequential(micro_batch[idx])
#          }
#          cdf_est <- isolate(values[["hermite_est"]]) %>%
#              cum_prob(x, clipped = TRUE)
#          df_cdf <- data.frame(x, cdf_est)
#          return(df_cdf)
#      })
#      updated_quantile_calc <- reactive({
#          values[["hermite_est"]]  %>% quant(input$percentile)
#      })
#      output$plot <- renderPlot({
#          ggplot(updated_cdf_calc(), aes(x = x)) + geom_line(aes(y = cdf_est)) +
#              ylab("Cumulative Probability")
#      }
#      )
#      output$quantile_text <- renderText({
#          return(paste(input$percentile * 100, "th Percentile:",
#                       round(updated_quantile_calc(), 2)))
#      })
#  }
#  shinyApp(ui = ui, server = server)

## -----------------------------------------------------------------------------
# Prepare Test Data
num_obs <-2000
test <- rchisq(num_obs,5)
test <- c(test,rlogis(num_obs))
test <- c(test,rnorm(num_obs))

## -----------------------------------------------------------------------------
# Calculate theoretical pdf, cdf and quantile values for comparison
x <- seq(-15,15,by=0.1)
actual_pdf_lognorm <- dchisq(x,5)
actual_pdf_logis <- dlogis(x)
actual_pdf_norm <- dnorm(x)
actual_cdf_lognorm <- pchisq(x,5)
actual_cdf_logis <- plogis(x)
actual_cdf_norm <- pnorm(x)
p <- seq(0.05,0.95,by=0.05)
actual_quantiles_lognorm <- qchisq(p,5)
actual_quantiles_logis <- qlogis(p)
actual_quantiles_norm <- qnorm(p)

## -----------------------------------------------------------------------------
# Construct Hermite Estimator 
h_est <- hermite_estimator(N=20,standardize = T,exp_weight_lambda = 0.005)

## -----------------------------------------------------------------------------
# Loop through test data and update h_est to simulate observations arriving 
# sequentially
count <- 1
res <- data.frame()
res_q <- data.frame()
for (idx in c(1:length(test))) {
  h_est <- h_est %>% update_sequential(test[idx])
  if (idx %% 100 == 0){
    if (floor(idx/num_obs)==0){
      actual_cdf_vals <- actual_cdf_lognorm
      actual_pdf_vals <-actual_pdf_lognorm
      actual_quantile_vals <- actual_quantiles_lognorm
    }
    if (floor(idx/num_obs)==1){
      actual_cdf_vals <- actual_cdf_logis
      actual_pdf_vals <-actual_pdf_logis
      actual_quantile_vals <- actual_quantiles_logis
    }
    if (floor(idx/num_obs)==2){
      actual_cdf_vals <- actual_cdf_norm
      actual_pdf_vals <- actual_pdf_norm
      actual_quantile_vals <- actual_quantiles_norm
    }
    idx_vals <- rep(count,length(x))
    cdf_est_vals <- h_est %>% cum_prob(x, clipped=T)
    pdf_est_vals <- h_est %>% dens(x, clipped=T)
    quantile_est_vals <- h_est %>% quant(p)
    res <- rbind(res,data.frame(idx_vals,x,cdf_est_vals,actual_cdf_vals,
                                pdf_est_vals,actual_pdf_vals))
    res_q <- rbind(res_q,data.frame(idx_vals=rep(count,length(p)),p,
                                    quantile_est_vals,actual_quantile_vals))
    count <- count +1
  }
}
res <- res %>% mutate(idx_vals=idx_vals*100)
res_q <- res_q %>% mutate(idx_vals=idx_vals*100)

## ----eval=FALSE---------------------------------------------------------------
#  # Visualize Results for PDF (Not run, requires gganimate, gifski and transformr
#  # packages)
#  p <- ggplot(res,aes(x=x)) + geom_line(aes(y=pdf_est_vals, colour="Estimated")) + geom_line(aes(y=actual_pdf_vals, colour="Actual")) +
#    scale_colour_manual("",
#                        breaks = c("Estimated", "Actual"),
#                        values = c("blue", "black")) + ylab("Probability Density") +transition_states(idx_vals,transition_length = 2,state_length = 1) +
#    ggtitle('Observation index {closest_state}')
#  anim_save("pdf.gif",p)

## ----eval=FALSE---------------------------------------------------------------
#  # Visualize Results for CDF (Not run, requires gganimate, gifski and transformr
#  # packages)
#  p <- ggplot(res,aes(x=x)) + geom_line(aes(y=cdf_est_vals, colour="Estimated")) + geom_line(aes(y=actual_cdf_vals, colour="Actual")) +
#    scale_colour_manual("",
#                        breaks = c("Estimated", "Actual"),
#                        values = c("blue", "black")) +
#    ylab("Cumulative Probability") +
#    transition_states(idx_vals, transition_length = 2,state_length = 1) +
#    ggtitle('Observation index {closest_state}')
#  anim_save("cdf.gif", p)

## ----eval=FALSE---------------------------------------------------------------
#  # Visualize Results for Quantiles (Not run, requires gganimate, gifski and
#  # transformr packages)
#  p <- ggplot(res_q,aes(x=actual_quantile_vals)) +
#    geom_point(aes(y=quantile_est_vals), color="blue") +
#    geom_abline(slope=1,intercept = 0) +xlab("Theoretical Quantiles") +
#    ylab("Estimated Quantiles") +
#    transition_states(idx_vals,transition_length = 2, state_length = 1) +
#    ggtitle('Observation index {closest_state}')
#  anim_save("quant.gif",p)

## ----eval=FALSE---------------------------------------------------------------
#  citation("hermiter")

