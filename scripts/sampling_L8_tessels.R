# script for sampling firecci tessels in order to decide which are to
# be analyzed

# Step 1: Extract the fraction of burnt area for each tessel ----
#

plan_sampling <- drake_plan(
    # library(firecci),
    in_ba_file = file_in(file.path(here::here(), "data-raw/2018-Annual-ESACCI-L4_FIRE-BA-MODIS-fv5.1_km2.tif")),
    in_tessels = file_in(file.path(here::here(), ("data-raw/L8_tessels.gpkg"))),
    # out_file   = file.path(here::here(), ("data/L8_tessels_with_ba.gpkg")),
    intensity  = firecci::extract_intensity(in_ba_file,
                                            in_tessels,
                                            file.path(here::here(), ("data/L8_tessels_with_ba.gpkg"))),
    outfile = {
        out_file = file.path(here::here(), ("data/L8_tessels_with_ba.gpkg"))
        sf::st_write(intensity, file_out(out_file))
        },
    report = rmarkdown::render(
        knitr_in("rmds/report_sampling.Rmd"),
        output_file = file_out(file.path(here::here(), "reports/report.html"))
        # output_file = file_out(file.path(system.file("reports/", package = "firecci"),
        #                                  "report_sampling.html")),
        #                        quiet = TRUE
        # )

    )
)

ba_data  <- sf::st_read(out_file)
seq_data <- sf::st_read("data-raw/World2018_Landsat/WorldLandsat2018_16Days_2.shp") %>%
    dplyr::select(pathrow, TimeSerie, N.Images, CloudCover, Start_Date, End_Date) %>%
    dplyr::mutate(pathrow = stringr::str_pad(pathrow, 6, "left", 0))
names(seq_data) <- c("pathrow", "ts_length", "n_images", "cloudcov", "startdate", "enddate", "geometry")


in_data <- ba_data %>%
    dplyr::left_join(st_drop_geometry(seq_data)) %>%
    dplyr::mutate(biome_code = as.character(dplyr::recode_factor(biome_code ,"0" = "Desert",
                                                    "1" = "Tropical Forests",
                                                    "2" = "Temperate Broadleaved",
                                                    "3" = "Coniferous Forest",
                                                    "4" = "Savannas",
                                                    "5" = "Temp. Grasslands",
                                                    "6" = "Mediterranean"))) %>%
    dplyr::mutate(lat = st_coordinates(st_geometry(sf::st_centroid(.)))[,2]) %>%
    # dplyr::filter(!is.na(n_images)) %>%
    dplyr::arrange(biome_code, ts_length)


tochange <- which(in_data$biome_code == "Desert" & in_data$lat > 40)
in_data$biome_code[tochange] <- "Tundra"

    # dplyr::mutate(biome_code = ifelse(biome_code == "Tundra and Desert" & lat > 40,
    #                                   "Tundra", biome_code))

out_data <- in_data %>%
    dplyr::group_by(biome_code, ts_length) %>%
    dplyr::tally()

ggplot2::ggplot(out_data) +
    ggplot2::geom_line(aes(x = ts_length, y = n)) +
    ggplot2::facet_wrap(~biome_code) +
    ggplot2::ggtitle("Number of L8 tessels as a function of Length of continuous serie and Olsen Biome - 16days separation") +
    ggplot2::theme_light() +
    ggplot2::xlab("Length of continuos serie [Days]") +
    ggplot2::ylab("Number of L8 Tessels") +
    ggplot2::geom_vline(aes(xintercept = 120), col = "red", linetype = "dashed")

thresh <- 120
myfunc <- function(x, thresh) {
    x_sub <- dplyr::filter(x, ts_length > thresh)
    # %>%
    #     st_drop_geometry()
    out <- x_sub %>%
        dplyr::group_by(biome_code) %>%
        dplyr::tally() %>%
        dplyr::mutate(thresh = thresh) %>%
        sf::st_cast("MULTIPOLYGON")
    out
}
x = as.list(seq(0, 320, 16))
out_temp <- lapply(x, FUN = function(x) myfunc(in_data, x))

    out <- data.table::rbindlist(out_temp)%>%
        sf::st_as_sf() %>%
        sf::st_cast("POLYGON")

tm_shape(out) + tm_polygons(col = "biome_code") + tm_facets("thresh")

ggplot2::ggplot(out_temp) +
    ggplot2::geom_line(aes(x = thresh, y = n)) +
    ggplot2::facet_wrap(~biome_code) +
    ggplot2::ggtitle("Number of L8 tessels as a function of Length of continuous serie and Olsen Biome - 16days separation") +
    ggplot2::theme_light() +
    ggplot2::xlab("Length of continuos serie [Days]") +
    ggplot2::ylab("Number of L8 Tessels") +
    ggplot2::geom_vline(aes(xintercept = 120), col = "red", linetype = "dashed") +
    ylim(0,300)

config <- drake_config(plan_sampling)
vis_drake_graph(config)
make(plan_sampling)

intensity <- readd("intensity")

aa <- in_data %>%
    dplyr::filter(ts_length > 80)
