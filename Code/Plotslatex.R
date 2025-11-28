
ad1 <- evaluate_samplers(posterior_name = "earnings-logearn_interaction", chains = 2, main_iter = 20000, adapter = "adapter1", warmup_iter = 2000)
ad2 <- evaluate_samplers(posterior_name = "earnings-logearn_interaction", chains = 2, main_iter = 20000, adapter = "adapter2", warmup_iter = 2000)



plot_scenario_6panels <- function(posterior_name, res_ad1, res_ad2,
                                  file = paste0(posterior_name, "_sixpanel.pdf")) {
  
  # color map by sampler name (adjust if needed)
  all_names <- unique(c(colnames(res_ad1$esjd_curves),
                        colnames(res_ad1$dt_curves),
                        colnames(res_ad1$mse_curves),
                        colnames(res_ad2$esjd_curves),
                        colnames(res_ad2$dt_curves),
                        colnames(res_ad2$mse_curves)))
  base_map <- c(bimodal = "#9467bd", barker = "#1f77b4",
                amala   = "#2ca02c",  amh    = "#d62728",
                hmc     = "#7f7f7f")
  cols <- base_map[match(all_names, names(base_map))]
  if (anyNA(cols)) cols[is.na(cols)] <- seq_along(which(is.na(cols)))
  
  plot_panel <- function(M, title, ylab) {
    if (is.null(M) || ncol(M) == 0) { plot.new(); title(main = paste(title, "(n/a)")); return() }
    keep <- match(colnames(M), all_names)
    clr  <- cols[keep]
    matplot(M, type = "l", lty = 1, col = clr,
            xlab = "Markov chain iteration", ylab = ylab, main = title)
    if (ncol(M) <= 6) legend("topright", legend = colnames(M), lty = 1, col = clr, cex = 0.8, bty = "n")
  }
  
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  pdf(file, width = 9, height = 6)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1) + 0.2, oma = c(0, 0, 2, 0))
  
  # Row 1: Adapter 1
  plot_panel(res_ad1$esjd_curves, "Adapter 1 – ESJD", "ESJD (running avg)")
  plot_panel(res_ad1$dt_curves,   "Adapter 1 – DT",   expression(d[t]~"(log-variance distance)"))
  plot_panel(res_ad1$mse_curves,  "Adapter 1 – MSE",  "MSE 1st moment")
  
  # Row 2: Adapter 2
  plot_panel(res_ad2$esjd_curves, "Adapter 2 – ESJD", "ESJD (running avg)")
  plot_panel(res_ad2$dt_curves,   "Adapter 2 – DT",   expression(d[t]~"(log-variance distance)"))
  plot_panel(res_ad2$mse_curves,  "Adapter 2 – MSE",  "MSE 1st moment")
  
  mtext(posterior_name, outer = TRUE, line = 0.1, cex = 1.2)
  dev.off()
}

eli1 <- evaluate_samplers(posterior_name = "earnings-logearn_interaction", chains = 2, main_iter = 20000, adapter = "adapter1", warmup_iter = 2000)
eli2 <- evaluate_samplers(posterior_name = "earnings-logearn_interaction", chains = 2, main_iter = 20000, adapter = "adapter2", warmup_iter = 2000)

eli1 <- evaluate_samplers(posterior_name = "eight_schools-eight_schools_noncentered", chains = 2, main_iter = 20000, adapter = "adapter1", warmup_iter = 2000)
eli2 <- evaluate_samplers(posterior_name = "eight_schools-eight_schools_noncentered", chains = 2, main_iter = 20000, adapter = "adapter2", warmup_iter = 2000)

gg1 <- evaluate_samplers(posterior_name = "mesquite-logmesquite_logvolume", chains = 2, main_iter = 20000, adapter = "adapter1", warmup_iter = 2000)
gg2 <- evaluate_samplers(posterior_name = "mesquite-logmesquite_logvolume", chains = 2, main_iter = 20000, adapter = "adapter2", warmup_iter = 2000)


plot_scenario_6panels("eight_schools-eight_schools_centered", ad1, ad2)
plot_scenario_6panels("earnings-logearn_interaction", eli1, eli2)
plot_scenario_6panels("mesquite-logmesquite_logvolume", gg1, gg2)
