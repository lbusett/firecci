# Script to create plots of variability of long series as a function
# of biome and "allowed hole length"

library(magrittr)
lgt      <- 16
year     <- 2018
plot_list <- list()
out_plotfile  <- "data/results_sampling/plots_16.RData"
out_datafile  <- "data/results_sampling/data_16.RData"

out_plotfile  <- "data/results_sampling/plots_16_tttt.RData"
out_datafile  <- "data/results_sampling/data_16_tttt.RData"

in_ba_file <- file.path(here::here(), ("data/burnt_area_L8_tessels.gpkg"))

if (lgt == 16) {
    in_lengths_file <- file.path(here::here(), ("data-raw/World2018_Landsat/WorldLandsat2018_16Days_2.shp"))
} else {
    if (lgt == 32) {
        in_lengths_file <- file.path(here::here(), ("data-raw/World2018_Landsat/WorldLandsat2018_32Days_2.shp"))
    } else {
        in_lengths_file <- file.path(here::here(), ("data-raw/World2018_Landsat/WorldLandsat2018_48Days.shp"))
    }
}

ba_data  <- sf::st_read(in_ba_file, stringsAsFactors = FALSE)
plot_list$tot_ba <- tmap::tm_shape(ba_data) +
    tmap::tm_polygons(col = "tot_ba", style = "log10_pretty")

# seq_data <- sf::st_read(in_lengths_file, stringsAsFactors = FALSE) %>%
#     dplyr::select(pathrow, TimeSerie, N.Images, CloudCover, Start_Date, End_Date) %>%
#     dplyr::mutate(pathrow = stringr::str_pad(pathrow, 6, "left", 0))
seq_data <- sf::st_read(in_lengths_file, stringsAsFactors = FALSE) %>%
    dplyr::select(pathrow, TimeSerie, N.Images, CloudCover, Start_Date, End_Date) %>%
    dplyr::mutate(pathrow = stringr::str_pad(pathrow, 6, "left", 0))
names(seq_data) <- c("pathrow", "ts_length", "n_images", "cloudcov", "startdate", "enddate", "geometry")

# join the two datasets and compute cumulate burnt areas columns

in_data <- ba_data %>%
    # Join with sequences and add infpo on Biome ----
    dplyr::left_join(sf::st_drop_geometry(seq_data)) %>%
    dplyr::mutate(biome_code = as.character(dplyr::recode_factor(biome_code ,"0" = "Other",
                                                                 "1" = "Trop. Forest",
                                                                 "2" = "Temp. Forest",
                                                                 "3" = "Boreal Forest",
                                                                 "4" = "Savannas",
                                                                 "5" = "Temp. Grasslands",
                                                                 "6" = "Medit. Forests"))) %>%
    dplyr::mutate(biome_code = forcats::fct_relevel(biome_code, c("Boreal Forest",
                                                                  "Temp. Forest",
                                                                  "Temp. Grasslands",
                                                                  "Medit. Forests",
                                                                  "Savannas",
                                                                  "Trop. Forest",
                                                                  "Other"))) %>%
    # sort for readability
    dplyr::mutate(lat = sf::st_coordinates(sf::st_geometry(sf::st_centroid(.)))[,2]) %>%
    dplyr::arrange(biome_code, ts_length) %>%
    dplyr::mutate(ts_length = ifelse(ts_length > 200, 200, ts_length)) %>%
    # Compute cumulated burnt area by biome (on sequence) and normalize by ----
    # total cumulated
    dplyr::group_by(biome_code) %>%
    dplyr::arrange(biome_code, tot_ba) %>%
    dplyr::mutate(cum_ba = cumsum(tot_ba)) %>%
    dplyr::mutate(norm_cum_ba = cum_ba / max(cum_ba)) %>%
    # add the number of thiessen in biome and "posityion" along the sequence of cumulates ----
    dplyr::add_tally(name = "N") %>%
    dplyr::mutate("N_seq" = seq_along(biome_code)) %>%
    dplyr::ungroup() %>%
    sf::st_as_sf() %>%
    dplyr::select( dplyr::everything(), geom)


# identify thresholds for distinction between high and low BA ("position" and threshold ----
# of area ----
thresh_data <- in_data %>%
    dplyr::group_by(biome_code) %>%
    dplyr::summarize(thresh_pos = dplyr::first(which(norm_cum_ba > 0.2)),
                     thresh_ba  = tot_ba[dplyr::first(which(norm_cum_ba > 0.2))]) %>%
    sf::st_drop_geometry()

# Join the dateaset with the thresholds, and assign FI on basis of tot_ba ----
# above or below thresh_ba
in_data <- in_data %>%
    dplyr::left_join(thresh_data) %>%
    dplyr::mutate(FI = ifelse(tot_ba > thresh_ba, "High", "Low")) %>%
    dplyr::mutate(ts_length = ifelse(is.na(ts_length), 0 , ts_length))

plot_list$plot_thresh <- ggplot2::ggplot(in_data) +
    ggplot2::geom_line(ggplot2::aes(y = norm_cum_ba, x = N_seq)) +
    ggplot2::facet_wrap(~biome_code, scales = "free_x") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0.2),
                        col = "red", linetype = "dashed", alpha = 0.3) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = thresh_pos),
                        col = "red", linetype = "dashed", alpha = 0.3) +
    ggplot2::theme_light() +
    ggplot2::xlab("N") + ggplot2::ylab("Normalized Cumulated BA") +
    ggplot2::geom_text(data = thresh_data, ggplot2::aes(x = thresh_pos, y = 0.2,
                                                        label = format(thresh_ba, digits = 3)),
                       nudge_y = +0.1, color = "blue", size = 2.5) +
    ggplot2::ggtitle("Normalized Cumulated Burnt Area vs Biome",
                     subtitle = "numbers correspond to BA @ the 20th percentile of the distribution")

# plot maps of high/low FI as a function of biome
plot_list$plot_FI <- tmap::tm_shape(in_data) +
    tmap::tm_polygons(col = "FI", palette = "-Dark2")
# +
#     tmap::tm_facets("biome_code")

#  create plot of histogram of length of series, by biome

hist_data <- in_data %>%
    dplyr::group_by(biome_code, ts_length) %>%
    dplyr::tally()


plot_list$plot_hist <- ggplot2::ggplot(in_data) +
    ggplot2::geom_histogram(ggplot2::aes(x = ts_length),
                            binwidth = 16, color = "black", position = "dodge") +
    ggplot2::stat_bin(binwidth=16, geom="text", colour="blue", size=2.5,
                      ggplot2::aes(x = ts_length, label=..count.., y=1.1*(..count..))) +
    # ggplot2::geom_col(ggplot2::aes(x = as.factor(n_images), y = n)) +
    ggplot2::facet_wrap(~biome_code) +
    ggplot2::theme_light() +
    ggplot2::ggtitle("Number of tessels as a function of continuous series length") +
    ggplot2::xlab("Length of continuous series") +
    ggplot2::ylab(glue::glue("Number of tessels (year = {year} - maxdist = {lgt})")) +
    ggplot2::scale_x_continuous(breaks = seq(0,200,16)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))

plot_list$plot_hist_facet <- ggplot2::ggplot(in_data) +
    ggplot2::geom_histogram(ggplot2::aes(x = ts_length, fill = FI),
                            binwidth = 16, color = "black", position = "dodge") +
    ggplot2::stat_bin(binwidth=16, geom="text", colour="black", size=2.5,
                      ggplot2::aes(x = ts_length, label=..count.., group=FI, y=1.2*(..count..))) +
    # ggplot2::geom_col(ggplot2::aes(x = as.factor(n_images), y = n)) +
    ggplot2::facet_grid(biome_code~FI, scale = "free_y") +
    # ggplot2::geom_text(ggplot2::aes(x = biome_code)) +
    ggplot2::theme_light() +
    ggplot2::ggtitle("Number of tessels as a function of continuous series length") +
    ggplot2::xlab("Length of continuous series") +
    ggplot2::ylab(glue::glue("Number of tessels (year = {year} - maxdist = {lgt})")) +
    ggplot2::scale_x_continuous(breaks = seq(0,200,16)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))+
    ggplot2::theme(legend.position = "None")


# All FI togheter ----
#
#  create plot of number of tessels above a given length
tallyfunc <- function(x, thresh) {
    x_sub <- dplyr::filter(x, ts_length > thresh)
    # %>%
    #     sf::st_drop_geometry()
    out <- x_sub %>%
        dplyr::group_by(biome_code) %>%
        dplyr::tally() %>%
        dplyr::mutate(thresh = thresh) %>%
        sf::st_cast("MULTIPOLYGON")
    out
}
x = as.list(seq(0, 200, 16))
tallydata <- lapply(x, FUN = function(x) tallyfunc(in_data, x)) %>%
    data.table::rbindlist() %>%
    dplyr::filter(thresh > 0) %>%
    sf::st_as_sf()

# tallydata$thresh_fct <- glue::glue("lgt = {tallydata$thresh}")

plot_list$mapFI <- tmap::tm_shape(tallydata) +
    tmap::tm_polygons(col = "biome_code") +
    tmap::tm_facets("thresh")

tallydata <- sf::st_drop_geometry(tallydata)
plot_list$plot_number <- ggplot2::ggplot(tallydata) +
    ggplot2::geom_line(ggplot2::aes(x = thresh, y = n)) +
    ggplot2::facet_wrap(~biome_code) +
    ggplot2::theme_light() +
    ggplot2::ggtitle("Number of tessels with continuous length above a given threshold",
                     subtitle = "Each panel indicates how many tessels could be \"available\" for sampling
                     when using the corresponding threshold on the minimum length of the serie") +
    ggplot2::xlab("Length of continuous series") +
    ggplot2::ylab(glue::glue("Number of tessels with length above threshold (year = {year} - maxdist = {lgt})")) +
    ggplot2::scale_x_continuous(breaks = seq(0,200,16)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))

plot_list$plot_number_bar <- ggplot2::ggplot(tallydata) +
    ggplot2::geom_col(ggplot2::aes(x = biome_code, y = n)) +
    ggplot2::geom_text(ggplot2::aes(x = biome_code, y = 1.2*n, label = n),
                       color = "blue", size = 2.5) +
    ggplot2::facet_wrap(~thresh, ncol = 3, scales = "free_y") +
    ggplot2::theme_light() +
    ggplot2::ggtitle("Number of tessels with continuous length above a given threshold - All Data",
                     subtitle = "Each panel indicates how many tessels could be \"available\" for sampling
                     when using the corresponding threshold on the minimum length of the serie") +
    ggplot2::ylab(glue::glue("Number of tessels with length above threshold (year = {year} - maxdist = {lgt})")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,  hjust = 1)) +
    ggplot2::xlab("") +
    ggplot2::theme(legend.position = "None")

plot_list$thmap  <- tmap::tm_shape(tallydata) +
    tmap::tm_polygons(col = "biome_code") +
    tmap::tm_facets("thresh")

# Analysis by FI ----

tallyfunc <- function(x, thresh) {
    x_sub <- dplyr::filter(x, ts_length > thresh) %>%
        sf::st_drop_geometry()
    out <- x_sub %>%
        dplyr::group_by(biome_code, FI) %>%
        dplyr::tally() %>%
        dplyr::mutate(thresh = thresh)
    out
}
x = as.list(seq(0, 200, 16))
tallydata <- lapply(x, FUN = function(x) tallyfunc(in_data, x)) %>%
    data.table::rbindlist() %>%
    dplyr::filter(thresh > 0)

plot_list$plot_number_FI <- ggplot2::ggplot(tallydata) +
    ggplot2::geom_line(ggplot2::aes(x = thresh, y = 1.2*n, color = FI)) +
    ggplot2::facet_grid(FI~biome_code, scales = "free_y") +
    ggplot2::theme_light() +
    ggplot2::ggtitle("Number of tessels with continuous length above a given threshold",
                     subtitle = glue::glue("Indicates how many tessels could be \"available\" for sampling when using different",
                                           " thresholds on the minimum length of the serie")) +
    ggplot2::xlab("Length of continuous series") +
    ggplot2::ylab(glue::glue("Number of tessels with length above threshold (year = {year} - maxdist = {lgt})")) +
    ggplot2::scale_x_continuous(breaks = seq(0,200,16)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))

plot_list$plot_number_bar_FI_H <- ggplot2::ggplot(dplyr::filter(tallydata, FI == "High")) +
    ggplot2::geom_col(ggplot2::aes(x = biome_code, y = n)) +
    ggplot2::geom_text(ggplot2::aes(x = biome_code, y = 1.2 * n, label = n),
                       color = "blue", size = 2.5) +
    ggplot2::facet_wrap(~thresh, ncol = 3, scales = "free_y") +
    ggplot2::theme_light() +
    ggplot2::ggtitle("Number of tessels with continuous length above a given threshold - High FI",
                     subtitle = glue::glue("Each panel indicates how many tessels could be \"available\" for sampling when using the corresponding",
                                           " threshold on the minimum length of the serie")) +
    ggplot2::ylab(glue::glue("Number of tessels with length above threshold (year = {year} - maxdist = {lgt})")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,  hjust = 1)) +
    ggplot2::xlab("") +
    ggplot2::theme(legend.position = "None")

plot_list$plot_number_bar_FI_L <- ggplot2::ggplot(dplyr::filter(tallydata, FI == "Low")) +
    ggplot2::geom_col(ggplot2::aes(x = biome_code, y = n)) +
    ggplot2::geom_text(ggplot2::aes(x = biome_code, y = 1.2 * n, label = n),
                       color = "blue", size = 2.5) +
    ggplot2::facet_wrap(~thresh, ncol = 3, scales = "free_y") +
    ggplot2::theme_light() +
    ggplot2::ggtitle("Number of tessels with continuous length above a given threshold - Low FI",
                     subtitle = glue::glue("Each panel indicates how many tessels could be \"available\" for sampling when using the corresponding",
                                           " threshold on the minimum length of the serie")) +
    ggplot2::ylab(glue::glue("Number of tessels with length above threshold (year = {year} - maxdist = {lgt})")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,  hjust = 1)) +
    ggplot2::xlab("") +
    ggplot2::theme(legend.position = "None")

save(plot_list, file = out_plotfile)
save(in_data,   file = out_datafile)

# plot_number_2_FI_L + p200lot_number_2_FI_H
#
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

