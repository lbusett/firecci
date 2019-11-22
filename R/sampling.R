# Script used for sampling the tessels to be used for validation.
#
# Data used are the results of the analysis done using the "plot_longseries.R"
# script, for a given length.
#
# The script first computes the number of sample to be taken for each stratum
# (biom x fire intensity) following the approach in Padilla et al., (2017)
# Then, it identifies a minimum value for the time series length that allows
# to get a smapling ratio above the specified one, for each stratum.
# Finally, it filters the tessels based on the determined ts_length threshsold,
# and then performs a random sampling in each stratum, with N equal to nh.


library(dplyr)
library(ggplot2)
library(data.table)

in_datafile <- "data/results_sampling/data_16.RData"
in_data     <- get(load(in_datafile))
out_sampfile <- "data/results_sampling/samp_16.RData"
min_lgt     <- 0
sampling_ratio <- 5
frac_limit <- 10
# fix the seed to allow reproducibility - change this for a new
# sampling
set.seed(42)

# Remove tessels almost completely in water
land_vector <- rnaturalearth::ne_download(scale    = 50,
                                          type     = 'land',
                                          category = "physical") %>%
    sf::st_as_sf()

intersected <- in_data %>%
    dplyr::mutate(totarea = as.numeric(sf::st_area(.))) %>%
    sf::st_intersection(land_vector)

fracland <- intersected %>%
    dplyr::mutate(int_area = as.numeric(sf::st_area(.))) %>%
    dplyr::group_by(pathrow) %>%
    dplyr::mutate(fracland = 100 * sum(int_area) / totarea[[1]]) %>%
    dplyr::select(pathrow, fracland) %>%
    sf::st_drop_geometry()

in_data <- in_data %>%
    dplyr::left_join(fracland) %>%
    dplyr::filter(fracland >= frac_limit) %>%
    dplyr::filter(!is.na(fracland))

# Identify number of cases to be sampled, by BIOME and FIRE INTENSITY -----
# (Following Padella et al., 2017)
n_samples   <- 100
nsamp_tbl <- in_data %>%
    dplyr::select(pathrow, biome_code, tot_ba, tot_area, ts_length, FI) %>%
    dplyr::filter(ts_length >= min_lgt) %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(biome_code, FI) %>%
    dplyr::summarize(Nh = dplyr::n(), BAh = mean(tot_ba, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(NhBah = Nh * (BAh ^ 0.5)) %>%
    dplyr::mutate(proportion = NhBah / sum(NhBah)) %>%
    dplyr::mutate(nh = round(n_samples * proportion))

#iteration #1
n_samples <- 96
nsamp_tbl <- in_data %>%
    dplyr::select(pathrow, biome_code, tot_ba, tot_area, ts_length, FI) %>%
    dplyr::filter(ts_length >= min_lgt) %>%
    dplyr::mutate(tot_ba = ifelse(biome_code == "Medit. Forests", 0 , tot_ba)) %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(biome_code, FI) %>%
    dplyr::summarize(Nh = dplyr::n(), BAh = mean(tot_ba, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(NhBah = Nh * (BAh ^ 0.5)) %>%
    dplyr::mutate(proportion = NhBah / sum(NhBah)) %>%
    dplyr::mutate(nh = n_samples * proportion) %>%
    dplyr::mutate(nh = round(ifelse(nh == 0, 2, nh)))

# after identification of number of samples, do the actual sampling -----

# First of all, let's see how long the TS can be in order to have at least a sampling_ratioX
# selection (i.e. if I have to sample 2, and a sampling_ratio of 5, I need at least a
# population of 10)


find_length <- function(in_data, nsamp_tbl, lgt, ratio){

    # in_data2 <- dplyr::left_join(in_data, samp_data)
    subdata <- in_data %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(ts_length >= lgt) %>%
        dplyr::group_by(biome_code, FI) %>%
        dplyr::tally() %>%
        dplyr::left_join(nsamp_tbl) %>%
        dplyr::mutate(lgt = lgt, ratio = n/nh)
    subdata
}


out_lgt <- lapply(as.list(seq(1,200,1)),
                  FUN=function(x) find_length(in_data, nsamp_tbl, x, sampling_ratio)) %>%
    data.table::rbindlist()

plot_ratios <- ggplot2::ggplot(out_lgt) +
    theme_light() +
    geom_line(aes(x = lgt, y = ratio)) +
    facet_grid(FI~biome_code, scales = "free_y") +
    geom_hline(aes(yintercept = sampling_ratio), col = "red", linetype = "dashed")+
    geom_text(aes(x =150, y = 30, label = glue::glue("{nh} out of: {Nh}")), col = "blue")

# Compute the minimum ts_length that guarantees a sampling ratio above
# sampling_ratio
my_lgt <- out_lgt %>%
    dplyr::filter(ratio < sampling_ratio) %>%
    dplyr::summarize(minlgt = min(lgt))


# Now subset based on sampling_ratio, and extract the sample:

# subset based on computed min_length
sub_data <- in_data %>%
    dplyr::filter(ts_length >= my_lgt$minlgt)

# do the sampling for each group and regroup the rezsults at the end
ind <- 1
sampled_tessels <- list()
for (biome in unique(sub_data$biome_code)) {
    for (flev in (unique(sub_data$FI))) {
        n_samp <- nsamp_tbl %>%
            dplyr::filter(biome_code == biome & FI == flev) %>%
            .[["nh"]]
        data_to_samp <- sub_data %>%
            dplyr::filter(biome_code == biome & FI == flev)
        sample <- dplyr::sample_n(data_to_samp, n_samp)
        sampled_tessels[[ind]] <- sample %>%
            sf::st_cast("MULTIPOLYGON")
        ind <- ind + 1
    }
}
sampled_tessels <- data.table::rbindlist(sampled_tessels) %>%
    sf::st_as_sf()

tmap::tmap_mode("view")
mapplot <- tmap::tm_shape(sampled_tessels) +
    tmap::tm_polygons(col = "biome_code") +
    tmap::tm_shape(land_vector) +
    tmap::tm_borders(col = "red")

mapplot_facet <- tmap::tm_shape(sampled_tessels) +
    tmap::tm_polygons(col = "biome_code") +
    tmap::tm_facets("FI") +
    tmap::tm_borders()

save(nsamp_tbl, sampled_tessels, out_lgt,
     plot_ratios,
     mapplot, mapplot_facet,
     file = out_sampfile)

# Sampling qithout square root correction ----
#
# n_samples   <- 100
# samp_data_nosqrt <- in_data %>%
#     dplyr::select(pathrow, biome_code, tot_ba, tot_area, ts_length, FI) %>%
#     dplyr::filter(ts_length >= min_lgt) %>%
#     sf::st_drop_geometry() %>%
#     dplyr::group_by(biome_code, FI) %>%
#     dplyr::summarize(Nh = dplyr::n(), BAh = mean(tot_ba, na.rm = TRUE)) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(NhBah = Nh * (BAh)) %>%
#     dplyr::mutate(proportion = NhBah / sum(NhBah)) %>%
#     dplyr::mutate(nh = round(100 * proportion))
#
# n_samples = 100-12
# samp_data_nosqrt <- in_data %>%
#     dplyr::select(pathrow, biome_code, tot_ba, tot_area, ts_length, FI) %>%
#     dplyr::filter(ts_length >= min_lgt) %>%
#     dplyr::mutate(tot_ba = ifelse(biome_code %in% c("Medit. Forests", "Boreal Forest", "Temp. Forest"), 0 , tot_ba)) %>%
#     sf::st_drop_geometry() %>%
#     dplyr::group_by(biome_code, FI) %>%
#     dplyr::summarize(Nh = dplyr::n(), BAh = mean(tot_ba, na.rm = TRUE)) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(NhBah = Nh * (BAh)) %>%
#     dplyr::mutate(proportion = NhBah / sum(NhBah)) %>%
#     dplyr::mutate(nh = round(n_samples * proportion)) %>%
#     dplyr::mutate(nh = ifelse(nh == 0, 2, nh))
#
# save(samp_data, samp_data_nosqrt, file = out_sampfile)

# samp_data <- in_data %>%
#     dplyr::group_by(biome_code, FI) %>%
#     sf::st_drop_geometry() %>%
#     dplyr::summarize(Nh = n(),
#                      cum_ba = sum(tot_ba)) %>%
#     dplyr::mutate(prop_fact = Nh * sqrt(cum_ba)) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(wgt = prop_fact / sum(prop_fact)) %>%
#     dplyr::select(1,2,3,5,6)
#
# aaa <- in_data %>%
#     dplyr::left_join(samp_data) %>%
#     sf::st_drop_geometry()
# #
# # st <- aaa %>%
# dplyr::group_by(biome_code, FI) %>%
#     # mutate(num_rows=n()) %>%
#     dplyr::sample_n(100 * wgt, weight=wgt)
