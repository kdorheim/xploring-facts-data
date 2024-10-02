# Running Hector with parameter samples from Matilda
library(hector)
library(ggplot2)
library(ncdf4)
library(tidyverse)
library(grid)
library(gridExtra)

# Load in parameter csv
params_all <- read.csv(file.path("data","hector_params.csv"))
colnames(params_all) <- list("beta","q10_rh","npp_flux0","alpha","diff","S")

# Get 200 random rows from params_all
rand_rows <- sample.int(nrow(params_all),200) #(replace=FALSE)
params_200 <- params_all[rand_rows,]
## Load in 200 params that work for certain
#params_200 <- read.csv(file.path("analysis","sample_params.csv"))
#rownames(params_200) <- params_200$X
#params_200$X <- NULL

# Prepare hector core
ini_file <- system.file("input/hector_ssp585.ini",package="hector")
core <- newcore(ini_file)

# Loop through Hector runs with all 200 samples
# PROBLEM: some samples will error out (Assertion failed: Flux and pool values may not be negative in ?)
# Temp solution: just get a new set of random samples, since almost all run just fine. Possible to-do item: get which exact samples don't work?
start.time <- Sys.time()
# Create empty dataframe in which we'll store outputs
h_results <- data.frame(scenario = character(),
                        year = double(),
                        variable = character(),
                        value = double(),
                        units = character(),
                        sample = integer())
# TODO (in general throughout code): switch for loops for lapply/similar functions?
for (sample in 1:200) {
    for (parameter in colnames(params_200)) {
        setvar(core, NA, parameter, params_200[[parameter]][sample], getunits(parameter))
    }
    reset(core)
    run(core)
    # Get results
    result <- fetchvars(core, 1750:2300, list(GMST(),HEAT_FLUX()))
    # Convert heat flux to ocean heat content
    result[result$variable==HEAT_FLUX(),]$value <- cumsum(result[result$variable==HEAT_FLUX(),]$value)*5.10065e14*0.71*31556930
    result[result$variable==HEAT_FLUX(),]$units <- "J"
    result[result$variable==HEAT_FLUX(),]$variable <- "ohc"
    result[["sample"]] <- sample
    # Append to results df
    h_results <- bind_rows(h_results,result)
}
h_results$model <- "hector"
end.time <- Sys.time()
total.time <- end.time-start.time

# Rename variables for plotting
h_plot_results <- h_results
h_plot_results[h_plot_results$variable=="gmst",]$variable <- "Global Mean Surface Temperature (degC)"
h_plot_results[h_plot_results$variable=="ohc",]$variable <- "Ocean Heat Content (J)"

h_climate <- ggplot(h_plot_results) +
    aes(x = year, y = value, group = sample) +
    geom_line(color="black",alpha=0.1) +
    facet_wrap(~variable, scales = "free_y") +
    ylab(NULL) +
    theme_bw()

# Histogram
h_results_2020 <- h_results[h_results$year==2020,]
h_results_2050 <- h_results[h_results$year==2050,]
h_results_2100 <- h_results[h_results$year==2100,]

h_gmst_2020 <- h_results_2020[h_results_2020$variable=="gmst",]
h_gmst_2050 <- h_results_2050[h_results_2050$variable=="gmst",]
h_gmst_2100 <- h_results_2100[h_results_2100$variable=="gmst",]

h_ohc_2020 <- h_results_2020[h_results_2020$variable=="ohc",]
h_ohc_2050 <- h_results_2050[h_results_2050$variable=="ohc",]
h_ohc_2100 <- h_results_2100[h_results_2100$variable=="ohc",]

# GMST histograms
h_gmst_2020_hist <- ggplot(h_gmst_2020,aes(x=value)) +
    geom_histogram(color="black",fill="gray") +
    labs(title="Hector GMST Values: 2020",x="degC", y = "Count") +
    theme_bw()

h_gmst_2050_hist <- ggplot(h_gmst_2050,aes(x=value)) +
    geom_histogram(color="black",fill="gray") +
    labs(title="Hector GMST Values: 2050",x="degC", y = "Count") +
    theme_bw()

h_gmst_2100_hist <- ggplot(h_gmst_2100,aes(x=value)) +
    geom_histogram(color="black",fill="gray") +
    labs(title="Hector GMST Values: 2100",x="degC", y = "Count") +
    theme_bw()

# Ocean heat content histograms
h_ohc_2020_hist <- ggplot(h_ohc_2020,aes(x=value)) +
    geom_histogram(color="black",fill="gray") +
    labs(title="Hector OHC Values: 2020",x="J", y = "Count") +
    theme_bw()

h_ohc_2050_hist <- ggplot(h_ohc_2050,aes(x=value)) +
    geom_histogram(color="black",fill="gray") +
    labs(title="Hector OHC Values: 2050",x="J", y = "Count") +
    theme_bw()

h_ohc_2100_hist <- ggplot(h_ohc_2100,aes(x=value)) +
    geom_histogram(color="black",fill="gray") +
    labs(title="Hector OHC Values: 2100",x="J", y = "Count") +
    theme_bw()

# Display all histograms
h_hist_grid <- grid.arrange(h_gmst_2020_hist,h_gmst_2050_hist,h_gmst_2100_hist,
                            h_ohc_2020_hist,h_ohc_2050_hist,h_ohc_2100_hist,
                            nrow = 2)
# Print histogram grid as pdf
pdf(file.path("analysis","plots","hector_climate_histograms.pdf"), height = 11, width = 8.5, paper = "letter")
grid.draw(h_hist_grid)
dev.off()

# Get FaIR climate to compare
nc_gmst <- nc_open(file.path("data","OUTPUTS","fair_results","tlm.offline","offline_gsat.nc"))
f_gmst <- ncvar_get(nc_gmst,"surface_temperature")
samples <- ncvar_get(nc_gmst,"samples")
years <- ncvar_get(nc_gmst,"years")
nc_close(nc_gmst)

nc_ohc <- nc_open(file.path("data","OUTPUTS","fair_results","tlm.offline","offline_ohc.nc"))
f_ohc <- ncvar_get(nc_ohc,"ocean_heat_content")
nc_close(nc_ohc)

# Get values in dataframe format to work with ggplot
colnames(f_gmst) <- samples
colnames(f_ohc) <- samples
df_gmst <- as.data.frame(f_gmst)
df_ohc <- as.data.frame(f_ohc)
df_gmst$year <- years
df_ohc$year <- years

f_gmst_long <- df_gmst %>%
    pivot_longer(
        cols = -year,
        names_to = "sample",
        values_to = "value"
    )
f_gmst_long$variable <- "gmst"
f_gmst_long$model <- "fair"

f_ohc_long <- df_ohc %>%
    pivot_longer(
        cols = -year,
        names_to = "sample",
        values_to = "value"
    )
f_ohc_long$variable <- "ohc"
f_ohc_long$model <- "fair"

# Put FaIR climate results into one dataframe and plot
f_results <- bind_rows(f_gmst_long,f_ohc_long)
f_results$sample <- as.integer(f_results$sample)
all_results <- bind_rows(h_results,f_results)

# Rename variables for plotting
f_plot_results <- f_results
f_plot_results[f_plot_results$variable=="gmst",]$variable <- "Global Mean Surface Temperature (degC)"
f_plot_results[f_plot_results$variable=="ohc",]$variable <- "Ocean Heat Content (J)"

f_climate <- ggplot(f_plot_results) +
    aes(x = year, y = value, group = sample) +
    geom_line(color="black",alpha=0.1) +
    facet_wrap(~variable, scales = "free_y") +
    theme_bw()

# Gather F and H together
all_plot_results <- bind_rows(h_plot_results,f_plot_results)
all_plot_results <- filter(all_plot_results,between(year,1750,2300)) # filter to Hector's shorter range

all_climate <- ggplot(all_plot_results,aes(x=year,y=value,group=interaction(sample,model),color=model)) +
    geom_line(alpha=0.3) +
    facet_wrap(~variable,scales="free_y") +
    theme_bw()
ggsave(file.path("analysis","plots","hector_v_fair.png"))

# Histogram comparing F and H for certain years
all_gmst_2020 <- filter(all_results,variable=="gmst") %>% filter(year==2020)
all_gmst_2050 <- filter(all_results,variable=="gmst") %>% filter(year==2050)
all_gmst_2100 <- filter(all_results,variable=="gmst") %>% filter(year==2100)

all_ohc_2020 <- filter(all_results,variable=="ohc") %>% filter(year==2020)
all_ohc_2050 <- filter(all_results,variable=="ohc") %>% filter(year==2050)
all_ohc_2100 <- filter(all_results,variable=="ohc") %>% filter(year==2100)

# GMST histograms
all_gmst_2020_hist <- ggplot(all_gmst_2020, aes(x=value, color=model, fill=model)) +
    geom_histogram(alpha=0.2, position="identity") +
    labs(title = "Global Mean Surface Temperature Comparison: 2020",x="degC",y="Count") +
    theme_bw()
all_gmst_2050_hist <- ggplot(all_gmst_2050, aes(x=value, color=model, fill=model)) +
    geom_histogram(alpha=0.2, position="identity") +
    labs(title = "Global Mean Surface Temperature Comparison: 2050",x="degC",y="Count") +
    theme_bw()
all_gmst_2100_hist <- ggplot(all_gmst_2100, aes(x=value, color=model, fill=model)) +
    geom_histogram(alpha=0.2, position="identity") +
    labs(title = "Global Mean Surface Temperature Comparison: 2100",x="degC",y="Count") +
    theme_bw()

# OHC histograms
all_ohc_2020_hist <- ggplot(all_ohc_2020, aes(x=value, color=model, fill=model)) +
    geom_histogram(alpha=0.2, position="identity") +
    labs(title = "Ocean Heat Content Comparison: 2020",x="degC",y="Count") +
    theme_bw()
all_ohc_2050_hist <- ggplot(all_ohc_2050, aes(x=value, color=model, fill=model)) +
    geom_histogram(alpha=0.2, position="identity") +
    labs(title = "Ocean Heat Content Comparison: 2050",x="degC",y="Count") +
    theme_bw()
all_ohc_2100_hist <- ggplot(all_ohc_2100, aes(x=value, color=model, fill=model)) +
    geom_histogram(alpha=0.2, position="identity") +
    labs(title = "Ocean Heat Content Comparison: 2100",x="degC",y="Count") +
    theme_bw()

comparison_hist_grid <- grid.arrange(all_gmst_2020_hist,all_gmst_2050_hist,all_gmst_2100_hist,
                                     all_ohc_2020_hist,all_ohc_2050_hist,all_ohc_2100_hist,
                                     nrow = 2)
# Print histogram grid as pdf
pdf(file.path("analysis","plots","comparison_climate_histograms.pdf"), height = 11, width = 8.5, paper = "letter")
grid.draw(comparison_hist_grid)
dev.off()


# Load in FaIR parameters for comparison with Hector's
nc_fairparam <- nc_open(file.path("analysis","fair_ar6_climate_params_v4.0.nc"))
fair_param_names <- names(nc_fairparam$var)
f_ecs <- ncvar_get(nc_fairparam, "F2x")
nc_close(nc_fairparam)

h_ecs <- params_all$S

# Put the two sets in one dataframe and plot histogram
h_ecs_df <- data.frame(ecs=h_ecs,model="hector")
f_ecs_df <- data.frame(ecs=f_ecs,model="fair")
ecs_all <- bind_rows(h_ecs_df,f_ecs_df)


ecs_hist <- ggplot(ecs_all, aes(x=ecs, color=model, fill=model)) +
    geom_histogram(alpha=0.2, position="identity") +
    ggtitle("Equilibrium Climate Sensitivity Parameter: Hector vs FaIR") +
    theme_bw()
ggsave(file.path("analysis","plots","ecs_histogram.png"))
