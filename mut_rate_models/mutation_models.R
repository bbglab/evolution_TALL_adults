library(VGAM)
library(scales)



one_step_path = "" # output path for 1 step simulations
linear_path = ""  # output path for linear step simulations
hpsc_file = "" # signature 5 contributions from data of Osorio et al., 2018 Cell PMID: 30485801
divergence_path = "" #output for the divergence results
data_file = "../intermediate_files/signature_counts.tsv"

#####################################
### Plot and simulation functions ###
#####################################


c_data <- read.table(data_file, header = TRUE, sep = "\t")

c_data$DIAGNOSIS_AGE_DAYS <- c_data$DIAGNOSIS_AGE_YEARS * 365.25
 

plot_simulations <- function(c_data, sim_patients, mut_normal = 12 , path, plot_mean = TRUE){
  for (c_patient in names(sim_patients)){
    print(c_patient)
    pdf(file = paste(path , c_patient,"_allsims", ".pdf", sep = ""), width = 5, height = 5)
    par(mfrow = c(5,4), bty = "l", las = 1, mar = c(3,3,3,3), cex = 0.2, lwd = 0.5)
    for(years in names(sim_patients[[c_patient]])){
      plot(NA, ylim = c(0,0.03), main = paste(c_patient, ":  ", years, " years acceleration\nMut_rate: ", mean(sim_patients[[c_patient]][[years]]$mut_rate)),
           xlab = "Mutations", ylab = "Frequency",
           xlim= summary(c(sim_patients[[c_patient]][[years]]$normal_muts + rowSums(sim_patients[[c_patient]][[years]]$Primary),
                           sim_patients[[c_patient]][[years]]$normal_muts + rowSums(sim_patients[[c_patient]][[years]]$Relapse),
                           c_data$TRUNK_SBS5[c_data$PATIENT == c_patient] + c_data$PRIVATE_PRY_SBS5[c_data$PATIENT == c_patient],
                           c_data$TRUNK_SBS5[c_data$PATIENT == c_patient] + c_data$PRIVATE_REL_SBS5[c_data$PATIENT == c_patient]))[c(1,6)])
      polygon(density(sim_patients[[c_patient]][[years]]$normal_muts + rowSums(sim_patients[[c_patient]][[years]]$Primary)), col = "#2c7fb890", border = NA)
      polygon(density(sim_patients[[c_patient]][[years]]$normal_muts + rowSums(sim_patients[[c_patient]][[years]]$Relapse)), col = "#fd8d3c90", border = NA)
      abline(v=c_data$TRUNK_SBS5[c_data$PATIENT == c_patient] + c_data$PRIVATE_PRY_SBS5[c_data$PATIENT == c_patient])
      abline(v=c_data$TRUNK_SBS5[c_data$PATIENT == c_patient] + c_data$PRIVATE_REL_SBS5[c_data$PATIENT == c_patient], lty = 3)
    }
    dev.off()

    years <- sapply(sim_patients[[c_patient]], function(x){
      x_pry <- sum((x$normal_muts + rowSums(x$Primary)) < (c_data$TRUNK_SBS5 + c_data$PRIVATE_PRY_SBS5)[c_data$PATIENT == c_patient])/length(x$normal_muts)
      x_pry <- min(x_pry, 1 - x_pry)
      x_rel <- sum((x$normal_muts + rowSums(x$Relapse)) < (c_data$TRUNK_SBS5 + c_data$PRIVATE_REL_SBS5)[c_data$PATIENT == c_patient])/length(x$normal_muts)
      x_rel <- min(x_rel, 1 - x_rel)
      x_pry * x_rel
    })
    years <- names(years)[grep(max(years, na.rm = TRUE), years)]

    if(length(years) > 1) next
    pdf(file = paste(path , c_patient, ".pdf", sep = ""), width = 5, height = 5)
    par(mfrow = c(1,2), bty = "l", las = 1, mar = c(5,5,5,5), cex = 0.3, lwd = 0.5)
    plot(NA, ylim = c(0,0.03), main = paste(c_patient, ":  ", years, " years acceleration\nMut_rate: ", mean(sim_patients[[c_patient]][[years]]$mut_rate)),
         xlab = "Mutations", ylab = "Frequency",
         xlim= summary(c(sim_patients[[c_patient]][[years]]$normal_muts + rowSums(sim_patients[[c_patient]][[years]]$Primary),
                         sim_patients[[c_patient]][[years]]$normal_muts + rowSums(sim_patients[[c_patient]][[years]]$Relapse),
                         c_data$TRUNK_SBS5[c_data$PATIENT == c_patient] + c_data$PRIVATE_PRY_SBS5[c_data$PATIENT == c_patient],
                         c_data$TRUNK_SBS5[c_data$PATIENT == c_patient] + c_data$PRIVATE_REL_SBS5[c_data$PATIENT == c_patient]))[c(1,6)])
    polygon(density(sim_patients[[c_patient]][[years]]$normal_muts + rowSums(sim_patients[[c_patient]][[years]]$Primary)), col = "#2c7fb890", border = NA)
    polygon(density(sim_patients[[c_patient]][[years]]$normal_muts + rowSums(sim_patients[[c_patient]][[years]]$Relapse)), col = "#fd8d3c90", border = NA)
    abline(v=c_data$TRUNK_SBS5[c_data$PATIENT == c_patient] + c_data$PRIVATE_PRY_SBS5[c_data$PATIENT == c_patient])
    abline(v=c_data$TRUNK_SBS5[c_data$PATIENT == c_patient] + c_data$PRIVATE_REL_SBS5[c_data$PATIENT == c_patient], lty = 3)

    c_pat <- c_data[c_data$PATIENT == c_patient,]
    plot(NA, ma = c_patient, xlab = "Years", ylab = "Mutations",
         xlim = c(0,c_pat$DIAGNOSIS_AGE_DAYS + c_pat$PRIMARY_TO_RELAPSE_AGE_DAYS)/365.25,
         ylim = c(0,c_pat$TRUNK_SBS5 + max(c_pat$PRIVATE_PRY_SBS5, c_pat$PRIVATE_REL_SBS5)))
    abline(1,mut_normal, lty = 2)
    abline(h=c_pat$TRUNK_SBS5, lty = 3)
    for(i in 1:nrow(sim_patients[[c_patient]][[years]]$Primary)){
      y_offset <- c_pat$DIAGNOSIS_AGE_DAYS/365.25 - as.numeric(years)
      lines(x= c(0,y_offset),
            y = c(0,sim_patients[[c_patient]][[years]]$normal_muts[i]), col = "#00683720")
      lines(x = y_offset + (0:ncol(sim_patients[[c_patient]][[years]]$Primary))/365.25,
            y = c(0, sapply(1:ncol(sim_patients[[c_patient]][[years]]$Primary), function(x) sum(sim_patients[[c_patient]][[years]]$Primary[i, 1:as.numeric(x)]))) +
                    sim_patients[[c_patient]][[years]]$normal_muts[i], col = "#2c7fb820")
      lines(x = y_offset + (0:ncol(sim_patients[[c_patient]][[years]]$Relapse))/365.25,
            y = c(0, sapply(1:ncol(sim_patients[[c_patient]][[years]]$Relapse), function(x) sum(sim_patients[[c_patient]][[years]]$Relapse[i, 1:as.numeric(x)]))) +
                    sim_patients[[c_patient]][[years]]$normal_muts[i], col = "#fd8d3c20")
    }
    if(plot_mean == TRUE){
      y_offset <-  c_pat$DIAGNOSIS_AGE_DAYS/365.25 - as.numeric(years)
      lines(x = y_offset + (0:ncol(sim_patients[[c_patient]][[years]]$Primary))/365.25,
            y = c(0, sapply(1:ncol(sim_patients[[c_patient]][[years]]$Primary), function(x) sum(colMeans(sim_patients[[c_patient]][[years]]$Primary)[1:as.numeric(x)]))) +
                  mean(sim_patients[[c_patient]][[years]]$normal_muts), col = "#2c7fb8", lwd = 2)
      lines(x = y_offset + (0:ncol(sim_patients[[c_patient]][[years]]$Relapse))/365.25,
            y = c(0, sapply(1:ncol(sim_patients[[c_patient]][[years]]$Relapse), function(x) sum(colMeans(sim_patients[[c_patient]][[years]]$Relapse)[1:as.numeric(x)]))) +
                  mean(sim_patients[[c_patient]][[years]]$normal_muts), col = "#fd8d3c", lwd = 2)
    }
    points(x= c(c_pat$DIAGNOSIS_AGE_DAYS/365.25, (c_pat$DIAGNOSIS_AGE_DAYS + c_pat$PRIMARY_TO_RELAPSE_AGE_DAYS)/365.25),
           y = c(c_pat$PRIVATE_PRY_SBS5, c_pat$PRIVATE_REL_SBS5) + c_pat$TRUNK_SBS5, col = c("#2c7fb8", "#fd8d3c"), pch = 16)

    n_mut <- c_pat$TRUNK_SBS5 - mean(sim_patients[[c_patient]][[years]]$normal_muts)
    i <- 0
    while(n_mut > 0){
      i <- i + 1
      n_mut <- n_mut - mean(c(sim_patients[[c_patient]][[years]]$Primary[,i], sim_patients[[c_patient]][[years]]$Relapse[,i]))
    }
    c_intersect <- ncol(sim_patients[[c_patient]][[years]]$Primary) - i

    text(x = (c_pat$DIAGNOSIS_AGE_DAYS - c_intersect)/365.25,
         y = c_pat$TRUNK_SBS5, labels = paste(round(c_intersect, digits = 2), "days")
    )
    dev.off()
  }
}

one_step <- function(c_data, mut_normal = 12, years_range = 10, iter = 100, path = ".",
                     rho = 0.0002, draw_plot = TRUE ){
  sim_patients <- apply(c_data[!is.na(c_data$PRIMARY_TO_RELAPSE_AGE_DAY),], 1, function(x){
    print(x["PATIENT"])
    out <- list()
    for(years in seq(0.1,years_range,0.1)){
      out[[as.character(years)]][["normal_muts"]] <- rbetabinom(iter, round(as.numeric(x["DIAGNOSIS_AGE_DAYS"]) - 365.25 * years), mut_normal/365.25, rho)
      remaining_shared <- as.numeric(x["TRUNK_SBS5"]) - as.numeric(out[[as.character(years)]][["normal_muts"]])
      remaining_shared <- pmax(0, remaining_shared)
      mut_rate <- ((as.numeric(x["PRIVATE_PRY_SBS5"]) + remaining_shared) / (years * 365.25) +
                  (as.numeric(x["PRIVATE_REL_SBS5"]) + remaining_shared) /(years * 365.25 + as.numeric(x["PRIMARY_TO_RELAPSE_AGE_DAYS"]))) / 2
      mut_rate[mut_rate < mut_normal/365.25 ] <- mut_normal/365.25
      out[[as.character(years)]][["mut_rate"]] <- mut_rate * 365.25
      c_factor <- ceiling(max(mut_rate))
      out[[as.character(years)]][["Primary"]] <- vapply(1:round(365.25 * years), function(y){
        rbetabinom(iter,c_factor,mut_rate / c_factor, rho)
      }, numeric(iter))

      out[[as.character(years)]][["Relapse"]] <- vapply(1:round(365.25 * years + as.numeric(x["PRIMARY_TO_RELAPSE_AGE_DAYS"])), function(y){
        rbetabinom(iter, c_factor, mut_rate / c_factor, rho)
      }, numeric(iter))
    }
    out
  })

  names(sim_patients) <- c_data$PATIENT[!is.na(c_data$PRIMARY_TO_RELAPSE_AGE_DAY)]
  if(draw_plot == TRUE) {
    plot_simulations(c_data, sim_patients, mut_normal, path)
  }
  return(sim_patients)
}

linear_growth <- function(c_data, mut_normal = 12, years_range = 10, iter = 100, path = ".",
                          rho = 0.0002, draw_plot = TRUE){
  sim_patients <- apply(c_data[!is.na(c_data$PRIMARY_TO_RELAPSE_AGE_DAY),], 1, function(x){
    print(x["PATIENT"])
    out <- list()
    for(years in seq(0.1,years_range,0.1)){
      out[[as.character(years)]][["normal_muts"]] <- rbetabinom(iter, round(as.numeric(x["DIAGNOSIS_AGE_DAYS"]) - 365.25 * years), mut_normal/365.25, rho)
      remaining_shared <- as.numeric(x["TRUNK_SBS5"]) - as.numeric(out[[as.character(years)]][["normal_muts"]])
      remaining_shared <- pmax(0, remaining_shared)

      mut_rate <- (2 * (as.numeric(x["PRIVATE_PRY_SBS5"]) + remaining_shared)/ (mut_normal * years *(years * 365.25 + 1)) +
                     2 * (as.numeric(x["PRIVATE_REL_SBS5"]) + remaining_shared) / ((mut_normal/365.25) * (years * 365.25 + as.numeric(x["PRIMARY_TO_RELAPSE_AGE_DAYS"])) * (years * 365.25 + as.numeric(x["PRIMARY_TO_RELAPSE_AGE_DAYS"]) + 1)))/2

      out[[as.character(years)]][["mut_rate"]] <- mut_rate ## increment rate per day -> mut_normal * (1+rt)
      c_factor <- ceiling(max((mut_normal/365.25)*(1 + mut_rate * as.numeric((365.25 * years + as.numeric(x["PRIMARY_TO_RELAPSE_AGE_DAYS"]))))))
      out[[as.character(years)]][["Primary"]] <- vapply(1:floor(365.25 * years), function(y){
        rbetabinom(iter,c_factor,((mut_normal/365.25)*(1 + mut_rate * as.numeric(y)))/c_factor, rho)
      }, numeric(iter))
      out[[as.character(years)]][["Relapse"]] <- vapply(1:floor(365.25 * years + as.numeric(x["PRIMARY_TO_RELAPSE_AGE_DAYS"])), function(y){
        rbetabinom(iter,c_factor,((mut_normal/365.25)*(1 + mut_rate * as.numeric(y)))/c_factor, rho)
      }, numeric(iter))
    }
    out
  })
  names(sim_patients) <- c_data$PATIENT[!is.na(c_data$PRIMARY_TO_RELAPSE_AGE_DAY)]
  if(draw_plot == TRUE){
  plot_simulations(c_data, sim_patients, mut_normal, path)
  }
  return(sim_patients)
}


##########################
### Running the models ###
##########################

hpsc <- read.table(hpsc_file, header = TRUE, stringsAsFactors = FALSE)

hpsc_corr <- function(rho, hpsc){
  c_lm <- summary(lm(SBS5_total~Age, data = hpsc))
  c_sbs5 <- vapply(hpsc$Age, function(x){
    rbetabinom(100, round(as.numeric(x)*365.25), c_lm$coefficients["Age","Estimate"]/365.25, rho) + c_lm$coefficients["(Intercept)","Estimate"]
  }, numeric(100))
  abs(c_lm$adj.r.squared - mean(apply(c_sbs5, 1, function(x) summary(lm(as.numeric(x)~hpsc$Age))$adj.r.squared)))
}

c_opt <- lapply(1:500, function(x) optimize(hpsc_corr, hpsc = hpsc,  lower = 0.00001, upper = 0.1, tol = 0.00001))

c_rho <- data.frame(minimum = sapply(c_opt, function(x) x$minimum),
                    objective = sapply(c_opt, function(x) x$objective))
c_rho <- c_rho[order(c_rho$objective),]

c_linear <- linear_growth(c_data, draw_plot = TRUE, rho = c_rho$minimum[1], path = linear_path)
c_1step <- one_step(c_data, draw_plot = TRUE, rho = c_rho$minimum[1], path = one_step_path)

save(c_linear, c_1step, file = paste(divergence_path, "simulations.RData", sep="/"))

##################################
### Estimating divergence time ###
##################################


divergence_time <- list()

for (c_model in c("linear", "1step")){
  if(c_model == "linear"){
    sim_patients <- c_linear
  }else{
    sim_patients <- c_1step
  }
  for(i in as.character(c_data[!is.na(c_data$PRIMARY_TO_RELAPSE_AGE_DAY),"PATIENT"])){
    print(i)
    c_patient <- c_data[c_data$PATIENT == i,]
    divergence_time[[i]][[c_model]] <- data.frame(
                              Patient = i,
                              Model = c_model,
                              Years_acc = names(sim_patients[[i]]),
                              Likelihood = sapply(sim_patients[[i]], function(x) {
                                p_density <- density(x$normal_muts + rowSums(x$Primary))
                                if(c_patient$TRUNK_SBS5 + c_patient$PRIVATE_PRY_SBS5 > max(p_density$x) | c_patient$TRUNK_SBS5 + c_patient$PRIVATE_PRY_SBS5 < min(p_density$x)){
                                  return(0)
                                }else{
                                  p_density <- p_density$y[which.min(abs(p_density$x - c_patient$TRUNK_SBS5 - c_patient$PRIVATE_PRY_SBS5))]
                                }
                                r_density <- density(x$normal_muts + rowSums(x$Relapse))
                                if(c_patient$TRUNK_SBS5 + c_patient$PRIVATE_REL_SBS5 > max(r_density$x) | c_patient$TRUNK_SBS5 + c_patient$PRIVATE_REL_SBS5 < min(r_density$x)){
                                  return(0)
                                }else{
                                  r_density <- r_density$y[which.min(abs(r_density$x - c_patient$TRUNK_SBS5 - c_patient$PRIVATE_REL_SBS5))]
                                }
                                return(p_density * r_density)
                              }),
                              Divergence = sapply(sim_patients[[i]], function(x) {
                                                    if(c_patient$TRUNK_SBS5 < mean(x$normal_muts)){ #during nomal mut_rate phase
                                                      return(c_patient$DIAGNOSIS_AGE_DAYS - c_patient$TRUNK_SBS5 * 365.25 / 12)
                                                    }
                                                    n_mut <- c_patient$TRUNK_SBS5 - mean(x$normal_muts)
                                                    j <- 0
                                                    while(n_mut > 0 && ncol(x$Primary) > j + 1){
                                                      j <- j + 1
                                                      n_mut <- n_mut - mean(c(x$Primary[,j], x$Relapse[,j]))
                                                    }
                                                    return((ncol(x$Primary) - j))
                                                  }),
                              stringsAsFactors = FALSE
                          )
  }
}

divergence_time <- lapply(divergence_time, function(x) do.call(rbind, x))

for (i in seq_along(divergence_time)){
  divergence_time[[i]]$Likelihood[divergence_time[[i]]$Divergence == 1] <- 0
  divergence_time[[i]]$Likelihood[divergence_time[[i]]$Likelihood < sort(divergence_time[[i]]$Likelihood,decreasing = TRUE)[20] ] <- 0
  divergence_time[[i]]$Likelihood <- divergence_time[[i]]$Likelihood / sum(divergence_time[[i]]$Likelihood)
  divergence_time[[i]]$Likelihood[is.nan(divergence_time[[i]]$Likelihood)] <- 0
  divergence_time[[i]]$cex <- divergence_time[[i]]$Likelihood / max(divergence_time[[i]]$Likelihood)
}

divergence_time[["PAT5"]] <- NULL

c_means <- sapply(divergence_time, function(x) weighted.mean(x$Divergence, x$Likelihood))
c_variance <- sapply(names(divergence_time), function(x) sum(divergence_time[[x]]$Likelihood * ((divergence_time[[x]]$Divergence - c_means[x])^2)))
c_sd <- 1.96 * sqrt(c_variance / sapply(names(divergence_time), function(x) (sum(divergence_time[[x]]$Likelihood > 0) - 1)/sum(divergence_time[[x]]$Likelihood > 0)))
c_sd[is.nan(c_sd)] <- 0

c_models <- data.frame(Model = divergence_time[[1]]$Model, Years_acc = divergence_time[[1]]$Years_acc, Freq = 0, Patients = 0)

for (i in names(divergence_time)){
   c_models$Freq <- c_models$Freq + divergence_time[[i]]$Likelihood
   c_models$Patients <- c_models$Patients + (divergence_time[[i]]$Likelihood > 0)
}

for (i in names(divergence_time)){
  divergence_time[[i]]$Likelihood <- divergence_time[[i]]$Likelihood * c_models$Freq
  divergence_time[[i]]$Likelihood[divergence_time[[i]]$Likelihood < max(divergence_time[[i]]$Likelihood)/5 ] <- 0
  divergence_time[[i]]$Likelihood <- divergence_time[[i]]$Likelihood / sum(divergence_time[[i]]$Likelihood)
  divergence_time[[i]]$cex <- divergence_time[[i]]$Likelihood / max(divergence_time[[i]]$Likelihood)
}

svg(filename =  paste(divergence_path,"divergence_plot.svg", sep = "/"), width = 8, height = 5)
par(mfrow = c(1,1), bty = "n", las = 2, mar = c(4,4,4,4), cex.axis = 0.8, lwd = 1)
plot(NA, xlim = c(0,length(divergence_time)), ylim = c(900,0), xaxt = "n", ylab = "Days pre-diagnosis", xlab = "")
axis(3, at= 1:length(divergence_time), labels = names(divergence_time))
points(x = unlist(lapply(1:length(divergence_time), function(x) as.numeric(x) + runif(nrow(divergence_time[[as.numeric(x)]]), -0.3, 0.3))),
       y = unlist(lapply(divergence_time, function(x) x$Divergence)),
       pch = 16,
       cex = unlist(lapply(divergence_time, function(x) x$cex))*1.6,
       col = alpha(c(linear = "red", '1step' = "#006837")[unlist(lapply(divergence_time, function(x) x$Model))], 0.6)
      )
c_means <- sapply(divergence_time, function(x) weighted.mean(x$Divergence, x$Likelihood))
c_variance <- sapply(names(divergence_time), function(x) sum(divergence_time[[x]]$Likelihood * ((divergence_time[[x]]$Divergence - c_means[x])^2)))
c_sd <- 1.96 * sqrt(c_variance / sapply(names(divergence_time), function(x) (sum(divergence_time[[x]]$Likelihood > 0) - 1)/sum(divergence_time[[x]]$Likelihood > 0)))
c_sd[is.nan(c_sd)] <- 0
segments(x0 = 1:length(divergence_time) - 0.3, x1 = 1:length(divergence_time) + 0.3, y0 = c_means, y1 = c_means, lty = 12, lwd = 2, cex = 0.5)

dev.off()

svg(filename = paste(divergence_path,"divergence_plot_rev.svg", sep = "/"), width = 8, height = 5)
plot(NA, ylim = c(0,length(divergence_time)), xlim = c(0, 900), yaxt = "n", xlab = "Days pre-diagnosis", ylab = "")
axis(2, at= 1:length(divergence_time), labels = names(divergence_time), las = 2)
points(y = unlist(lapply(1:length(divergence_time), function(x) as.numeric(x) + runif(nrow(divergence_time[[as.numeric(x)]]), -0.3, 0.3))),
       x = unlist(lapply(divergence_time, function(x) x$Divergence)),
       pch = 16,
       cex = unlist(lapply(divergence_time, function(x) x$cex))*1.6,
       col = alpha(c(linear = "red", '1step' = "#006837")[unlist(lapply(divergence_time, function(x) x$Model))], 0.6)
)
segments(y0 = 1:length(divergence_time) - 0.3, y1 = 1:length(divergence_time) + 0.3, x0 = c_means, x1 = c_means, lty = 12, lwd = 2, cex = 0.5)

dev.off()

save(c_means, c_sd, divergence_time, file = paste(divergence_path, "results_divergence_time.RData", sep = "/"))

