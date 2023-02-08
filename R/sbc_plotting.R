source("R/5sm_data_generating_functions.R")

sbc_results_const_advi_m <- readRDS("data/sbc_results_constant_m.RDS")
sbc_results_const_advi_m %>% plot_sbc_bars(rates = "const", remove_pars = "S|total|raw") -> sbc_bars_const_m
ggsave(filename = "final_plots/sbc_advi_density_const_m.png", sbc_bars_const_m[[3]], device = "png", dpi = 600)
ggsave(filename = "final_plots/sbc_const_m.png", sbc_bars_const_m[[2]], device = "png", dpi = 600)


sbc_results_const_advi_w <- readRDS("data/sbc_results_constant_w.RDS")
sbc_results_const_advi_w %>% plot_sbc_bars(rates = "const", remove_pars = "S|total|raw") -> sbc_bars_const_w
ggsave(filename = "final_plots/sbc_advi_density_const_w.png", sbc_bars_const_w[[3]], device = "png", dpi = 600)
ggsave(filename = "final_plots/sbc_const_w.png", sbc_bars_const_w[[2]], device = "png", dpi = 600)

sbc_results_const_hmc_m <- readRDS("data/sbc_results_constant_m_hmc.RDS")
sbc_results_const_hmc_m %>% plot_sbc_bars(rates = "const", remove_pars = "S|total|raw") -> sbc_bars_const_m_hmc
ggsave(filename = "final_plots/sbc_hmc_density_const_m.png", sbc_bars_const_m_hmc[[3]], device = "png", dpi = 600)
ggsave(filename = "final_plots/sbc_hmc_const_m.png", sbc_bars_const_m_hmc[[2]], device = "png", dpi = 600)

sbc_results_const_hmc_w <- readRDS("data/sbc_results_constant_w_hmc.RDS")
sbc_results_const_hmc_w %>% plot_sbc_bars(rates = "const", remove_pars = "S|total|raw") -> sbc_bars_const_w_hmc
ggsave(filename = "final_plots/sbc_hmc_density_const_w.png", sbc_bars_const_w_hmc[[3]], device = "png", dpi = 600)
ggsave(filename = "final_plots/sbc_hmc_const_w.png", sbc_bars_const_w_hmc[[2]], device = "png", dpi = 600)


sbc_results_const_advi_m_meanfield <- readRDS("data/sbc_results_constant_m_meanfield.RDS")
sbc_results_const_advi_m_meanfield %>% 
  plot_sbc_bars(rates = "const", remove_pars = "S|total|raw") -> 
  sbc_results_const_advi_m_meanfield_bars

ggsave(filename = "final_plots/sbc_const_advi_m_meanfield_density.png", sbc_results_const_advi_m_meanfield[[3]], device = "png", dpi = 600)
ggsave(filename = "final_plots/sbc_const_advi_m_meanfield_sbc.png", sbc_results_const_advi_m_meanfield_bars[[2]], device = "png", dpi = 600)


sbc_results_const_advi_w_meanfield <- readRDS("data/sbc_results_constant_m_meanfield.RDS")
sbc_results_const_advi_w_meanfield %>% 
  plot_sbc_bars(rates = "const", remove_pars = "S|total|raw") -> 
  sbc_results_const_advi_w_meanfield_bars

ggsave(filename = "final_plots/sbc_const_advi_w_meanfield_density.png", sbc_results_const_advi_w_meanfield_bars[[3]], device = "png", dpi = 600)
ggsave(filename = "final_plots/sbc_const_advi_w_meanfield_sbc.png", sbc_results_const_advi_w_meanfield_bars[[2]], device = "png", dpi = 600)




ggpubr::ggarrange(sbc_results_const_advi_m_meanfield_bars[[2]], sbc_bars_const_m[[2]], sbc_bars_const_m_hmc[[2]], nrow = 3) ->
  all_sbc_men
ggsave(filename = "final_plots/all_sbc_men.png", all_sbc_men, device = "png", dpi = 600)


ggpubr::ggarrange(sbc_results_const_advi_w_meanfield_bars[[2]], sbc_bars_const_w[[2]], sbc_bars_const_w_hmc[[2]], nrow = 3) ->
  all_sbc_women
ggsave(filename = "final_plots/all_sbc_men.png", all_sbc_women, device = "png", dpi = 600)


