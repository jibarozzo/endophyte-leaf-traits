####################### RELATIVE ABUNDANCE  in E+ and E- ########
lifespan_names <- c(
  `CORDAL` = "CORDAL - 303 leaf days",
  `HYBAPR` = "HYBAPR - 393 leaf days",
  `VIROSU` = "VIROSU - 645 leaf days",
  `CALOLO` = "CALOLO - 860 leaf days", 
  `GUSTSU` = "GUSTSU -  1105 leaf days", 
  `HEISCO` = "HEISCO -  1750 leaf days"
)

########### gloming taxa and merging data to make plot ##

View(sample_data(ps2s))
#make a new phyloseq object that groups all taxa by class
ps3 <- tax_glom(ps2s, "class") 
#make it relative
ps0 <- transform_sample_counts(ps3, function(x) x / sum(x))
#merge samples so that they are all grouped/binned by species
ps1 <- merge_samples(ps0, "Species")
# repair factors
sample_data(ps1)$Species <- levels(sample_data(ps0)$Species)[get_variable(ps1, "Species")]
sample_data(ps1)$TimePoint <- levels(sample_data(ps0)$TimePoint)[get_variable(ps1, "TimePoint")]
#make relative again.
ps4 <- transform_sample_counts(ps1, function(x) x / sum(x))

rel_abundplot <- plot_bar(ps4, "Species", "Abundance") + geom_bar(aes(fill=class), stat="identity", position="stack")+
  labs(x = "Species", y = "Relative abundance (%)", fill="Class", labeller=as_labeller(lifespan_names)) + scale_x_discrete(limits=c("CORDAL", "HYBAPR", "VIROSU", "CALOLO", "GUSTSU", "HEISCO"))
rel_abundplot
ggsave("Aim2_RelativeAbundace_4April2022.jpeg", plot=rel_abundplot,dpi=600, units=c("mm"), width=190, height=180)
