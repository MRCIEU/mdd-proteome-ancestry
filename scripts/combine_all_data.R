#load the results for all ancestries into one data frame
#european ancenstry(eur)
load(file.path("results", "eur.rdata"))
eur <- res
eur$pop <- "EUR"

#african ancestry(afr)
load(file.path("results", "afr.rdata"))
afr<- res
afr$pop <- "AFR"

#central south asian ancestry(csa)
load(file.path("results", "csa.rdata"))
csa<- res
csa$pop <- "CSA"

#east asian ancestry(eas)
load(file.path("results", "eas.rdata"))
eas<- res
eas$pop <- "EAS"

#combine all data frames 

res_comb <- bind_rows(eur,eas,csa,afr)

head(res_comb)
View(res_comb)

#plot results of combined data

ggplot(res_comb, aes(y=pop, x=b, color=pop)) +
  geom_point() +
  geom_errorbarh(aes(xmin=b-se*1.96, xmax=b+se*1.96)) +
  geom_vline(xintercept=0, linetype="dashed") +
  labs(title = "Combined forest plot of potential MDD drug targets across ancestry",
       x = "beta_IV (logodds ratio per SD)",
       y = "Proteins") +
  facet_grid(exposure ~ .)+
  theme_minimal() +  # Apply a minimal theme
  theme(
    plot.title = element_text(hjust = 0.5), # Center the plot title
    axis.text.y=element_blank(), # remove y axis text
    strip.text.y=element_text(angle=0)
  )


# Display and save results

ggsave(file=file.path("results", "combined.png"))

