s2_ba_func <- function(){

    in_ba_file <-file.path(file.path(here::here(),
                                     "data-raw/2018-Annual-ESACCI-L4_FIRE-BA-MODIS-fv5.1_km2.tif"))
    in_tessels <- file.path(here::here(), ("data/S2tessels/AFRICA_S2_tessels.gpkg"))

intensity  = firecci::extract_intensity(
    in_ba_file,
    in_tessels,
    file.path(here::here(), ("data/S2tessels/AFRICA_S2_tessels_with_ba.gpkg")))
}
