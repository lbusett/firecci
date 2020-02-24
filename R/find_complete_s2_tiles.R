find_complete_s2_tiles <- function() {

    in_tiles <- readRDS("data-raw/s2_tiles.rds")
    in_footprints <- sf::st_read("data/S2tessels/AFRICA_S2_orbits.gpkg")

    int_tiles <- sf::st_intersects(in_tiles, in_footprints)
    lgt_int <- lengths(int_tiles)
    afr_tiles <- in_tiles[which(lgt_int > 0), ]

    centroids <- sf::st_coordinates(sf::st_centroid(afr_tiles))
    which_cent_ok <- which(centroids[,2] < 25)
    afr_tiles_centers <- afr_tiles[which_cent_ok,]  %>%
        sf::st_centroid() %>%
        dplyr::mutate(X = sf::st_coordinates(sf::st_centroid(.))[,1],
                      Y = sf::st_coordinates(sf::st_centroid(.))[,2])  %>%
        dplyr::mutate(zone = substring(tile_id, 1,2)) %>%
        dplyr::mutate(N_S = dplyr::if_else(Y > 0, 6, 7)) %>%
        dplyr::mutate(EPSG = as.numeric(paste0("32", N_S, zone)))

    square_func <- function(x) {
        out_square <- x %>%
            sf::st_transform(x$EPSG) %>%
            sf::st_buffer(50000.00, nQuadSegs = 100) %>%
            sf::st_bbox() %>%
            sf::st_as_sfc() %>%
            sf::st_cast("LINESTRING") %>%
            sf::st_segmentize(100) %>%
            sf::st_transform(4326) %>%
            sf::st_cast("POLYGON")
        out_square <- sf::st_sf(st_drop_geometry(x), geometry= st_sfc(out_square)) %>%
            dplyr::select(-X,-Y,-zone,-N_S)
        out_square
    }

    out_squares <- list()
    for (square in 1:dim(afr_tiles2)[1]) {
        out_squares[[square]] <- square_func(afr_tiles_centers[square,])
    }

    out_squares <- data.table::rbindlist(out_squares) %>%
        sf::st_as_sf()

    sing_10x10 = st_filter(out_squares, in_footprints, .predicate = st_within)
    sf::st_write(sing_10x10, "data/S2tiles/AFRICA_S2_single_tiles_10k.gpkg",
                 delete_dsn = TRUE)

    # find self intersections ----
    out_ints <- list()
    for (square in 1:dim(sing_10x10)) {
        print(square)
        ints <- sf::st_intersection(sing_10x10[square,], sing_10x10) %>%
            dplyr::mutate(area = sf::st_area(.)) %>%
            dplyr::arrange(area) %>%
            dplyr::mutate(int_area = as.numeric(area/1e+10)) %>%
            dplyr::filter(int_area > 0.1)
        out_ints[[square]] <- ints
    }
    out_ints_test <- data.table::rbindlist(out_ints) %>%
        sf::st_as_sf()

    out_ints_test <- out_ints_test %>%
        sf::st_collection_extract("POLYGON")

    sf::st_write(out_ints_test, "data/S2tiles/AFRICA_S2_single_tiles_10k_ints.gpkg",
                 delete_dsn = TRUE)

    out_ints_test <- sf::st_read("data/S2tiles/AFRICA_S2_single_tiles_10k_ints.gpkg")
    hist(out_ints_test$area)

    singles_complete <- dplyr::filter(out_ints_test, area > 9e09)
    singles_ints     <- dplyr::filter(out_ints_test, area < 9e09) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(tiletot = paste(sort(c(tile_id, tile_id.1)), collapse = "")) %>%
        dplyr::group_by(tiletot) %>%
        dplyr::slice(sample(1:2,1))

    singles_complete_final <- dplyr::filter(singles_complete,
                                            !(tile_id %in% unique(singles_ints$tile_id)))
    sf::st_write(singles_complete_final, "data/S2tiles/AFRICA_S2_single_tiles_10k_choosen.gpkg",
                 delete_dsn = TRUE)

    # Compute Burnt Area and Burnt Area fraction over tessels ----
    in_ba_file <- file.path(file.path(here::here(),
                                      "data-raw/2018-Annual-ESACCI-L4_FIRE-BA-MODIS-fv5.1_km2.tif"))
    in_tiles_file     <- file.path(here::here(),
                                   "data/S2tiles/AFRICA_S2_single_tiles_10k_choosen.gpkg")
    out_tiles_ba_file <- file.path(here::here(),
                                   "data/S2tiles/AFRICA_S2_tiles_small_with_ba.gpkg")

    intensity  = firecci::extract_intensity(
        in_ba_file,
        in_tiles_file,
        out_tiles_ba_file)

    # Now get info about biome and burnt area ----
    biome_file <- file.path(file.path(here::here(),
                                      "data/olson_biomes_reshuffled.gpkg"))
    biomes <- sf::st_read(biome_file) %>%
        dplyr::filter(!is.na(biome_name))

    in_tiles <- sf::st_read(out_tiles_ba_file)

    # intersect tiles with biome map ----
    int_biomes <- in_tiles %>%
        sf::st_intersection(biomes) %>%
        dplyr::mutate(int_area = as.numeric(sf::st_area(.)))

    # for each tessel, get the biome with greater area ----
    int_biomes2 <- int_biomes %>%
        dplyr::group_by(tile_id, biome_name) %>%
        dplyr::mutate(lc_area = sum(as.numeric(int_area))) %>%
        dplyr::group_by(tile_id) %>%
        dplyr::arrange(tile_id) %>%
        dplyr::top_n(1, lc_area) %>%
        dplyr::ungroup() %>%
        dplyr::select(tile_id,biome_name) %>%
        sf::st_drop_geometry()


    tiles_with_biome <- in_tiles %>%
        dplyr::left_join(int_biomes2)


    # get intersecting tessels
    africa_file <- file.path(here::here(), "data-raw/Africa.shp")
    africa_sub <- sf::st_read(africa_file, quiet = TRUE) %>%
        sf::st_set_crs(. ,4326) %>%
        sf::st_buffer(0) %>%
        dplyr::group_by() %>%
        dplyr::filter(ID != 689 & ID != 690 & !(ID %in% 27:41)) %>%
        dplyr::summarize() %>%
        sf::st_zm()
    which_intersect_land <- sf::st_intersects(tiles_with_biome, africa_sub)
    tiles_with_biome <- tiles_with_biome[which(lengths(which_intersect_land) != 0), ]
    # pl_voronoi <- tm_shape(voronoi) + tm_borders(col = "red")

    # remove tessels with a very small amount of land ----
    tiles_w_fc <- tiles_with_biome %>%
        sf::st_intersection(africa_sub) %>%
        dplyr::group_by(tile_id) %>%
        dplyr::summarise(tot_area = 10000*first(tot_area)) %>% #this "dissolves" the intersection by ID
        dplyr::mutate(area_int = sf::st_area(.)) %>%
        dplyr::mutate(fc = as.numeric(area_int / tot_area))
    which_above_10 <- which(tiles_w_fc$fc > 10)
    tiles_final  <- tiles_with_biome[which_above_10, ]
    tiles_final$fc <- tiles_w_fc$fc[which_above_10]
    tiles_final <- tiles_final %>%
        dplyr::select(1,2,3,4,5,7)

    if(dim(tiles_with_biome)[1] != dim(in_tiles)[1]) {
        stop("Something strange happened!")
    } else {
        out_tiles_file <- file.path(here::here(),
                                    "data/S2tiles/AFRICA_S2_tiles_ba_biome.gpkg")
        sf::st_write(tiles_with_biome, out_tiles_file, delete_dsn = TRUE)
    }

    # tile = 100
    # ints <- afr_tiles[tile, ] %>%
    #     dplyr::mutate(totarea = sf::st_area(.)) %>%
    #     sf::st_within(in_footprints) %>%
    #     dplyr::mutate(intarea = sf::st_area(.))
    #
    # int_tiles <- afr_tiles %>%
    #     dplyr::mutate(totarea = sf::st_area(.)) %>%
    #     sf::st_intersection(in_footprints) %>%
    #     dplyr::mutate(intarea = sf::st_area(.)) %>%
    #     dplyr::group_by(tile_id) %>%
    #     dplyr::add_count()
    #
    # sing_tiles <- int_tiles %>%
    #     dplyr::filter(n == 1) %>%
    #     sf::st_as_sf()


}
