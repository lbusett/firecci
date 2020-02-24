#' @title create_s2_tessel
#' @description Create datasets to be used for tessellation of S2 acquisitions
#'  over Africa ----
#' @return The function creates the datasets needed to perform tessellation of
#'  S2 acquisitions, with a latitudinal spacing equal to that used on L8 tessellation
#' @rdname create_s2_tessel
#' @export
#' @author Lorenzo Busetto, phD (2019)
#' @importFrom sf st_read st_zm st_set_crs st_buffer st_intersects st_intersection
#'  st_as_sfc st_bbox st_write st_cast st_geometry st_union st_graticule
#' @importFrom dplyr group_by filter summarize select mutate summarise_all
#' first rename summarise arrange row_number
#' @importFrom data.table rbindlist
#'

create_s2_tessel_small <- function(outdir) {
    require(tmap)
    if (missing("outdir")) {
        outdir <- file.path(here::here(), "data/S2tessels_small/")
        dir.create(outdir)
    } else {
        if(dir.exists(dirname(outdir))) {
            dir.create(outdir)
        } else {
            stop("`outdir` must specify the path to an existing folder",
                 "or to a non-existing subfolder of an existing folder!")
        }
    }

    # get a polygon vector of Africa ----
    africa_file <- file.path(here::here(), "data-raw/Africa.shp")
    africa <- sf::st_read(africa_file, quiet = TRUE) %>%
        sf::st_set_crs(. ,4326) %>%
        sf::st_buffer(0) %>%
        dplyr::group_by() %>%
        dplyr::filter(ID != 689 & ID != 690) %>%
        dplyr::summarize() %>%
        sf::st_zm()

    # Identify the S2 orbit footprints above africa
    fprints_file <- file.path(here::here(),
                              "data-raw/Sentinel-2A_MP_ACQ_KML_20190912T120000_20190930T150000.kml")
    fprints <- sf::st_read(fprints_file, "NOMINAL") %>%
        sf::st_zm()
    sub_fprints <- dplyr::filter(fprints, lengths(sf::st_intersects(fprints, africa)) > 0)
    sub_fprints <- sf::st_intersection(sub_fprints,
                                       sf::st_as_sfc(sf::st_bbox(africa) + c(-2,-2,+2,+2))) %>%
        dplyr::select(1:3,19,20,21,22) %>%
        dplyr::mutate(orbit = as.numeric(as.character(OrbitRelative))) %>%
        dplyr::group_by(orbit) %>%
        dplyr::summarise_all(dplyr::first)

    pl_foot <- tm_shape(sub_fprints) + tm_borders(col = "darkgreen")

    sf::st_write(sub_fprints,
                 file.path(here::here(), "data/S2tessels_small/AFRICA_S2_orbits.gpkg"),
                 delete_layer = TRUE)


    # Identify the S2 tiles footprints above africa ----
    tiles <- readRDS("data-raw/s2_tiles.rds")
    int_tiles <- sf::st_filter(tiles, sub_fprints)
    afr_tiles <- sf::st_filter(int_tiles, africa)
    sf::st_write(afr_tiles,
                 file.path(here::here(), "data/S2tessels_small/AFRICA_S2_tiles.gpkg"),
                 delete_layer = TRUE)

    # Import S2 orbit tracks and subset over Africa ----

    tracks <- sf::st_read(file.path(here::here(),
                                    "data-raw/S2A_rel_orbit_groundtrack_10Sec.kml")) %>%
        sf::st_zm() %>%
        dplyr::filter(Relative_Orbit %in% unique(sub_fprints$orbit)) %>%
        dplyr::rename(orbit = Relative_Orbit) %>%
        dplyr::mutate(orbit = as.numeric(as.character(orbit))) %>%
        sf::st_intersection(sf::st_as_sfc(sf::st_bbox(africa) + c(-2,-2,+2,+2))) %>%
        sf::st_cast("LINESTRING")

    # identify and remove "spurious" horizontal lines ----
    diffs <- lapply(sf::st_geometry(tracks),
                    FUN = function(x) data.frame(diff = sf::st_bbox(x)$ymin - sf::st_bbox(x)$ymax)) %>%
        data.table::rbindlist()
    which_not_horiz <- which(diffs < -0.15)

    # union the resulting segments
    tracks <- tracks[which_not_horiz,] %>%
        dplyr::group_by(orbit) %>%
        dplyr::summarise(geometry = sf::st_union(geometry))

    pl_track <- tm_shape(tracks) + tm_lines(col = "blue", lty = "dashed")

    # save ----
    sf::st_write(tracks,
                 file.path(here::here(),
                           "data/S2tessels_small/AFRICA_S2_tracks.gpkg"),
                 delete_layer = TRUE)

    # sample points along tracks ----
    #
    # Build lines along which to sample based on L8 path/rows, based on centroids  ----
    # of L8 tessels: this will allow to get the SAME sampling strategy on both
    # L8 and S2
    L8_lnestring <- sf::st_read(file.path(here::here(),
                                          "data-raw/L8_tessels.gpkg")) %>%
        dplyr::select(-area) %>%
        dplyr::mutate(l8row = substring(pathrow,4,6)) %>%
        sf::st_centroid() %>%
        dplyr::group_by(l8row) %>%
        dplyr::summarise() %>%
        sf::st_cast("LINESTRING") %>%
        dplyr::ungroup()

    # sample S2 orbits on L8 "rows" to get centroids ----
    points <- sf::st_intersection(tracks, L8_lnestring) %>%
        dplyr::select(orbit, l8row) %>%
        dplyr::group_by(orbit) %>%
        dplyr::arrange(orbit) %>%
        dplyr::mutate(nn = dplyr::row_number()) %>%
        dplyr::mutate(id = paste(orbit, nn, sep = "_")) %>%
        dplyr::select(-nn) %>%
        sf::st_union()

    pl_points <- tm_shape(points) + tm_dots(size = 0.01)

    outfile <- file.path(here::here(),
                         "data/S2tessels_small/AFRICA_S2_points.gpkg")
    sf::st_write(points, outfile, delete_layer = TRUE)

    # Create tessels ----

    envelope <- sf::st_convex_hull(points)
    voronoi  <- sf::st_voronoi(points, envelope) %>%
        sf::st_cast() %>%
        sf::st_intersection(envelope)
    voronoi <- voronoi %>%
        sf::st_sf() %>%
        dplyr::mutate(ID = 1:length(voronoi))

    # now get a "meaningful ID". Use the TM row and the S2 orbit ----
    # We derive them by bringing back the intersection points between
    # L8 rows and S2 orbits, and joining the attributes
    points_to_join <- sf::st_intersection(tracks, L8_lnestring) %>%
        dplyr::select(orbit, l8row)
    joined_pts_tbl <- voronoi %>%
        sf::st_intersection(points_to_join) %>%
        dplyr::arrange(ID) %>%
        sf::st_drop_geometry()

    voronoi <- voronoi %>%
        dplyr::left_join(joined_pts_tbl, by = "ID") %>%
        dplyr::mutate(ID_PR = paste(orbit, l8row, sep = "_")) %>%
        dplyr::select(everything(), geometry)

    # split tessels in two halves, to get "similar" area with L8 tessels ----
    small_voronoi <- voronoi %>%
        lwgeom::st_split(tracks) %>%
        sf::st_collection_extract("POLYGON") %>%
        dplyr::group_by(ID_PR) %>%
        dplyr::mutate(n = letters[1:dplyr::n()]) %>%
        dplyr::arrange(ID_PR) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(ID_PR_small = paste(ID_PR, n, sep = "_")) %>%
        dplyr::mutate(ID = 1:dplyr::n()) %>%
        dplyr::mutate(area = sf::st_area(.)) %>%
        dplyr::select(ID, area, orbit, l8row, ID_PR, ID_PR_small)

    voronoi <- small_voronoi
    # remove voronois not intersecting land ----

    # create a new "africa" vector, after removing cabo verde
    africa_sub <- sf::st_read(africa_file, quiet = TRUE) %>%
        sf::st_set_crs(. ,4326) %>%
        sf::st_buffer(0) %>%
        dplyr::group_by() %>%
        dplyr::filter(ID != 689 & ID != 690 & !(ID %in% 27:41)) %>%
        dplyr::summarize() %>%
        sf::st_zm()

    # get intersecting tessels
    which_intersect_land <- sf::st_intersects(voronoi, africa_sub)
    voronoi <- voronoi[which(lengths(which_intersect_land) != 0), ]
    pl_voronoi <- tm_shape(voronoi) + tm_borders(col = "red")
    outfile <- file.path(here::here(),
                         "data/S2tessels_small/AFRICA_S2_tessels_small_full.gpkg")
    sf::st_write(voronoi, outfile, delete_layer = TRUE)

    # compute centroids of voronois, to be able to subset below 25°N
    # Then, assign a unique identifier to be used later and compute
    # total area
    centroids <- sf::st_coordinates(sf::st_centroid(voronoi))
    which_cent_ok <- which(centroids[,2] < 25)
    voronoi <- voronoi[which_cent_ok,]

    # remove tessels with a very small amount of land ----
    vor_w_fc <- voronoi %>%
        sf::st_intersection(africa_sub) %>%
        dplyr::group_by(ID) %>%
        dplyr::summarise(area = first(area)) %>% #this "dissolves" the intersection by ID
        dplyr::mutate(area_int = sf::st_area(.)) %>%
        dplyr::mutate(fc = as.numeric(area_int / area))
    pl_vor_w_fc <- tm_shape(vor_w_fc) + tm_polygons(col = "fc", style = "fixed",
                                                    breaks = seq(0,1,0.1),
                                                    palette = "-Reds")

    which_above_10 <- which(vor_w_fc$fc > 0.10)
    voronoi_final  <- voronoi[which_above_10, ]
    voronoi_final$fc <- vor_w_fc$fc[which_above_10]

    pl_voronoi_sub <- tm_shape(voronoi_final) + tm_borders(col = "red")

    # save the final tessels ----
    outfile <- file.path(here::here(),
                         "data/S2tessels_small/AFRICA_S2_tessels_small.gpkg")
    sf::st_write(voronoi_final, outfile, delete_layer = TRUE)

    #save to R object
    outfile <- file.path(here::here(),
                         paste0("data/S2tessels_small/AFRICA_S2_tessels_small", Sys.Date(), ".RData"))
    save(pl_foot, pl_track, pl_points, pl_vor_w_fc, pl_voronoi, pl_voronoi_sub,
         voronoi_final, file = outfile)

    # Compute Burnt Area and Burnt Area fraction over tessels ----
    in_ba_file <- file.path(file.path(here::here(),
                                      "data-raw/2018-Annual-ESACCI-L4_FIRE-BA-MODIS-fv5.1_km2.tif"))
    in_tessels_file     <- file.path(here::here(),
                                     "data/S2tessels_small/AFRICA_S2_tessels_small.gpkg")
    out_tessels_ba_file <- file.path(here::here(),
                                     "data/S2tessels_small/AFRICA_S2_tessels_small_with_ba.gpkg")

    intensity  = firecci::extract_intensity(
        in_ba_file,
        in_tessels_file,
        out_tessels_ba_file)

    # Associate Biome to tessels ----

    # to eventually rebuild the "reshuffled" biome gpkg:
    biome_file <- file.path(file.path(here::here(),
                                      "data-raw/olson_biomes/official/wwf_terr_ecos.shp"))
    biomes <- sf::st_read(biome_file)
    biomes <- biomes %>%
        dplyr::mutate(biome_name = dplyr::case_when(
            BIOME %in% c(1, 2, 3)    ~ "Trop_Forest",
            BIOME %in% c(4, 5)       ~ "Temp_Forest",
            BIOME %in% c(6)          ~ "Boreal_Forest",
            BIOME %in% c(7)          ~ "Trop_Sav",
            BIOME %in% c(8, 9,10)    ~ "Temp_Grass",
            BIOME %in% c(11, 13, 14) ~ "Other",
            BIOME %in% c(12)         ~ "Medit_Forest",
            TRUE ~ as.character(NA)
        )) %>%
        sf::st_intersection(sf::st_as_sfc(sf::st_bbox(africa))) %>%
        dplyr::group_by(biome_name) %>%
        dplyr::summarise() %>%
        dplyr::ungroup()

    out_biome_file <- file.path(file.path(here::here(),
                                          "data/olson_biomes_reshuffled.gpkg"))
    sf::st_write(biomes, out_biome_file, delete_layer = TRUE)
    biome_file <- file.path(file.path(here::here(),
                                      "data/olson_biomes_reshuffled.gpkg"))

    biomes <- sf::st_read(biome_file) %>%
        dplyr::filter(!is.na(biome_name))

    in_tessels <- sf::st_read(out_tessels_ba_file)

    # intersect tessels with biome map ----
    int_biomes <- in_tessels %>%
        sf::st_intersection(biomes) %>%
        dplyr::mutate(int_area = as.numeric(sf::st_area(.)))

    # for each tessel, get the biome with greater area ----
    int_biomes2 <- int_biomes %>%
        dplyr::group_by(ID_PR_small, biome_name) %>%
        dplyr::mutate(lc_area = sum(as.numeric(int_area))) %>%
        dplyr::group_by(ID_PR_small) %>%
        dplyr::arrange(ID_PR_small) %>%
        dplyr::top_n(1, lc_area) %>%
        dplyr::ungroup() %>%
        dplyr::select(ID, ID_PR_small, biome_name) %>%
        sf::st_drop_geometry()


    tessels_with_biome <- in_tessels %>%
        dplyr::left_join(int_biomes2)  %>%
        dplyr::select(1,3,4,5,6,8,fc)

    if(dim(tessels_with_biome)[1] != dim(in_tessels)[1]) {
        stop("Something strange happened!")
    } else {
        out_tessels_file <- file.path(here::here(),
                                      "data/S2tessels_small/AFRICA_S2_tessels_small_ba_biome.gpkg")
        sf::st_write(tessels_with_biome, out_tessels_file, delete_dsn = TRUE)
    }

    # finally, create an associated file with 10*10 Km polygon over each tessels centroid ----

    tessels_centroids <- tessels_with_biome %>%
        sf::st_centroid() %>%
        dplyr::mutate(X = sf::st_coordinates(.)[,1],
                      Y = sf::st_coordinates(.)[,2])  %>%
        dplyr::mutate(zone = (floor((X + 180)/6) %% 60) + 1) %>%
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
    for (square in 1:dim(tessels_centroids)[1]) {
        out_squares[[square]] <- square_func(tessels_centroids[square,])
    }

    out_squares <- data.table::rbindlist(out_squares) %>%
        sf::st_as_sf() %>%
        dplyr::select(-fc)

    out_squares_file <- file.path(here::here(), ("data/S2tessels_small/AFRICA_S2_squares_smalltessels.gpkg"))
    sf::st_write(out_squares, out_squares_file, delete_dsn = TRUE)

    # old - based on fixed dist

    # points <- sf::st_intersection(tracks, sf::st_graticule(lat = seq(-90, 90, lat_dist),
    #                                                        lon = c(-180,180))) %>%
    #     dplyr::select(orbit, degree, type) %>%
    #     dplyr::group_by(orbit) %>%
    #     dplyr::arrange(orbit) %>%
    #     dplyr::mutate(nn = dplyr::row_number()) %>%
    #     dplyr::mutate(id = paste(orbit, nn, sep = "_")) %>%
    #     dplyr::select(-nn)

    # outfile <- file.path(outdir, paste0("AFRICA_S2_points_", lat_dist, "deg.gpkg"))
    # sf::st_write(points, outfile, delete_dsn = TRUE)

}
