k_fit_brms_GBER

indirect <- c(summary(posterior_samples_GBs$b_logER_scalelog_streamflow * posterior_samples_GBs$b_scaleparint3m_scalelog_streamflow * posterior_samples_GBs$b_scalelaketempC_scalelog_streamflow))


posterior_samples_GBs <- posterior_samples(k_fit_brms_GBER)

# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GBs))

# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GBs %>%
  select(b_logER_logERlag, b_logER_scalepar_int_3m, b_logER_scalelake_tempC, b_scalelaketempC_scalelog_streamflow,
         b_scaleparint3m_scalelog_streamflow,b_logER_scalelog_streamflow) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')


indirect <- c(summary(posterior_samples_BWs$b_logER_scalelog_streamflow * posterior_samples_BWs$b_scaleparint3m_scalelog_streamflow * posterior_samples_BWs$b_scalelaketempC_scalelog_streamflow))


posterior_samples_BWs <- posterior_samples(k_fit_brms_BWER)
# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BWs))

# Prepare the posterior samples data frame with shore levels
posterior_samples_BWq <- posterior_samples_BWs %>%
  select(b_logER_logERlag, b_logER_scalepar_int_3m, b_logER_scalelake_tempC, b_scalelaketempC_scalelog_streamflow,
         b_scaleparint3m_scalelog_streamflow,b_logER_scalelog_streamflow) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW <- melt(posterior_samples_BWq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')


SEM_df <- rbind(posterior_samples_melted_GB, posterior_samples_melted_BW)


# Plot histograms with credible intervals and mean estimates
site_semplot<-ggplot(SEM_df, aes(x = Estimate, fill = shore, color = shore)) +
  facet_wrap(~ Coefficient, scales = 'free') +
  geom_density(alpha = 0.75, linewidth = 0.25) +
  scale_fill_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  labs(title = 'Posterior Distributions of Coefficients by Site',
       x = 'Estimate',
       y = 'Density') +
  geom_vline(xintercept=0, linetype = "solid", size = 0.25, color = "black") 

# ggsave(plot = site_semplot, filename = paste("./figures/NS24_metab_streamER_semplot1.png",sep=""),width=6.5,height=3.75,dpi=300)



######## GPPP
indirect <- c(summary(posteriors_gb$b_logGPP_scalelog_streamflow * posteriors_gb$b_scalelaketempC_scalelog_streamflow * posteriors_gb$b_scaleparint3m_scalelog_streamflow))




posterior_samples_GBs_gpp <- posterior_samples(k_fit_brms_gb)

# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GBs_gpp))

# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq_gpp <- posterior_samples_GBs_gpp %>%
  select(b_logGPP_logGPPlag, b_logGPP_scalepar_int_3m, b_logGPP_scalelake_tempC, b_scalelaketempC_scalelog_streamflow,
         b_scaleparint3m_scalelog_streamflow, b_logGPP_scalelog_streamflow) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB_gpp <- melt(posterior_samples_GBq_gpp, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')


posterior_samples_BWs <- posterior_samples(k_fit_brms)

indirect <- c(summary(posterior_samples_BWs$b_logGPP_scalelog_streamflow * posterior_samples_BWs$b_scalelaketempC_scalelog_streamflow * posterior_samples_BWs$b_scaleparint3m_scalelog_streamflow))


# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BWs))

# Prepare the posterior samples data frame with shore levels
posterior_samples_BWq_gpp <- posterior_samples_BWs %>%
  select(b_logGPP_logGPPlag, b_logGPP_scalepar_int_3m, b_logGPP_scalelake_tempC, b_scalelaketempC_scalelog_streamflow,
         b_scaleparint3m_scalelog_streamflow, b_logGPP_scalelog_streamflow) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW_gpp <- melt(posterior_samples_BWq_gpp, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')


SEM_df_gpp <- rbind(posterior_samples_melted_GB_gpp, posterior_samples_melted_BW_gpp)


# Plot histograms with credible intervals and mean estimates
site_semplot_gpp<-ggplot(SEM_df_gpp, aes(x = Estimate, fill = shore, color = shore)) +
  facet_wrap(~ Coefficient, scales = 'free') +
  geom_density(alpha = 0.75, linewidth = 0.25) +
  scale_fill_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  labs(title = 'Posterior Distributions of Coefficients by Site',
       x = 'Estimate',
       y = 'Density') +
  geom_vline(xintercept=0, linetype = "solid", size = 0.25, color = "black") 

# ggsave(plot = site_semplot_gpp, filename = paste("./figures/NS24_metab_streamGPP_semplot1.png",sep=""),width=6.5,height=3.75,dpi=300)
