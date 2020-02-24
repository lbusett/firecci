sample_s2_tiles <- function() {
    require(dplyr)
    in_tilesfile <- "data/S2tiles/AFRICA_S2_tiles_small_with_ba.gpkg"
    tilesdata <- sf::st_read(in_tilesfile, stringsAsFactors = FALSE)

    # tm_shape(tilesdata) + tm_polygons(alpha = 0)
    # tt = "34MEU"

    result = list()
    for (tt in seq_along(unique(tilesdata$tile))) {
        # for (tt in seq_along(1:10)) {
        print(tt)
        curtile <- unique(tilesdata$tile)[tt]
        tiledata <- tilesdata %>%
            dplyr::filter(tile_id == curtile)

        aoi =  sf::st_geometry(sf::st_sf(sf::st_buffer(sf::st_centroid(tiledata[1,]),
                                                       0.0001)),
                               nQuadSegs = 40) %>%
            sf::st_buffer(0) %>%
            lwgeom::st_make_valid()

        s2_meta <- getSpatialData::getSentinel_query(
            c("2019-01-01", "2020-01-01"),
            aoi = aoi,
            platform = "Sentinel-2", username = "user", password = "user")

        s2_meta <- s2_meta %>%
            dplyr::filter(processinglevel == "Level-1C",
                          cloudcoverpercentage<=30,
                          tileid == curtile) %>%
            dplyr::select(tileid, datatakesensingstart, cloudcoverpercentage) %>%
            dplyr::mutate(date = as.Date(datatakesensingstart)) %>%
            dplyr::rename(clouds = cloudcoverpercentage) %>%
            dplyr::select(tileid, date, clouds) %>%
            dplyr::arrange(date)

        diffdata <- as.numeric(c(-difftime(s2_meta$date[1:(length(s2_meta$date)-1)] ,
                                           s2_meta$date[2:length(s2_meta$date)]), as.numeric(NA)))
        s2_meta$belowthresh <- diffdata <= diff_threshold
        s2_meta$diffdata    <- diffdata
        rle <- rle(s2_meta$belowthresh)

        lgts   <- list()
        starts <- list()
        ends   <- list()
        dates   <- list()
        for (seq_in in seq_along(rle$lengths)) {
            if (is.na(rle$values[seq_in]) || rle$values[seq_in] == FALSE) {
                lgts[[seq_in]]   <- -999
                starts[[seq_in]] <- -999
                ends[[seq_in]]   <- -999
            } else {
                # browser()
                if (seq_in == 1) {
                    start_ind  <- 1
                } else {
                    start_ind  <- cumsum(rle$lengths)[seq_in - 1] + 1
                }
                end_ind    <- cumsum(rle$lengths)[seq_in] + 1
                starts[[seq_in]] <- s2_meta$date[start_ind]
                ends[[seq_in]]   <- s2_meta$date[end_ind]
                lgts[[seq_in]]   <- as.numeric(ends[[seq_in]] -  starts[[seq_in]])
                dates[[seq_in]]  <- s2_meta$date[start_ind:end_ind]
            }
        }

        max_lgt       <- max(unlist(lgts))
        if (max_lgt != -999) {
            which_max_lgt <- which.max(unlist(lgts))
            dates_max_lgt <- as.Date(dates[which_max_lgt][[1]])

            out_data <- tiledata %>%
                tibble::as_tibble() %>%
                dplyr::mutate(start_date    = min(dates_max_lgt) ,
                              end_date      = max(dates_max_lgt),
                              lgt           = max_lgt,
                              n_dates_max_lgt = length(dates_max_lgt),
                              dates_max_lgt = list(dates_max_lgt)) %>%
                dplyr::select(everything(), attr(tilesdata, "sf_column"))

        } else {
            out_data <- tiledata %>%
                tibble::as_tibble() %>%
                dplyr::mutate(start_date    = NA ,
                              end_date      = NA,
                              lgt           = 0,
                              n_dates_max_lgt = 0,
                              dates_max_lgt = list(NA)) %>%
                dplyr::select(everything(), attr(tilesdata, "sf_column"))
        }
        result[[tt]] <- out_data
    }
    result_out <- data.table::rbindlist(result) %>%
        tibble::as.tibble()
    sf::st_write(result_out,
                 "data/S2tiles/AFRICA_S2_single_tiles_10k_choosen_length.gpkg",
                 delete_layer = TRUE)
    save(result_out, file = "data/S2tiles/AFRICA_S2_single_tiles_10k_choosen_length.RData")
}
