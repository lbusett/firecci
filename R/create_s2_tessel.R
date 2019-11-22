#' @title create_s2_tessel
#' @description Create datasets to be used for tessellation of S2 acquisitions
#'  over Africa ----
#' @param lat_dist `numeric [°]` Desired istance between tessel centroids (in degrees)
#' @return The function creates the datasets needed to perform tessellation of
#'  S2 acquisitions, with a given latitudinal spacing. Results are saved in
#'  data/S2tessels.
#' @details DETAILS
#' @examples
#' \dontrun{
#' create_s2_tessel(3)
#' }
#' @rdname create_s2_tessel
#' @export
#' @author Lorenzo Busetto, phD (2019)
#' @importFrom sf st_read st_zm st_set_crs st_buffer st_intersects st_intersection
#'  st_as_sfc st_bbox st_write st_cast st_geometry st_union st_graticule
#' @importFrom dplyr group_by filter summarize select mutate summarise_all
#' first rename summarise arrange row_number
#' @importFrom data.table rbindlist
#'

create_s2_tessel <- function(lat_dist, outdir = NULL) {

    if (is.null(outdir)) {
        outdir <- file.path(system.file(package = "firecci"), "data/S2tessels")
        dir.create(outdir)
    } else {
        if(dir.exists(dirname(outdir))) {
            dir.create(outdir)
        } else {
            stop("`outdir` must specify the path to an existing folder",
                 "or to a non-existing subfolder of an existing folder!")
        }
    }

    fprints_file <- system.file("/data-raw/Sentinel-2A_MP_ACQ_KML_20190912T120000_20190930T150000.kml",
                                package = "firecci")
    fprints <- sf::st_read(fprints_file, "NOMINAL") %>%
        sf::st_zm()

    africa_file <- system.file("/data-raw/Africa.shp",
                               package = "firecci")
    africa <- sf::st_read(africa_file, quiet = TRUE) %>%
        sf::st_set_crs(. ,4326) %>%
        sf::st_buffer(0) %>%
        dplyr::group_by() %>%
        dplyr::filter(ID != 689 & ID != 690) %>%
        dplyr::summarize() %>%
        sf::st_zm()

    sub_fprints <- dplyr::filter(fprints, lengths(sf::st_intersects(fprints, africa)) > 0) # will call filter.sf

    sub_fprints <- sf::st_intersection(sub_fprints, sf::st_as_sfc(sf::st_bbox(africa) + c(-2,-2,+2,+2))) %>%
        # dplyr::filter(Name != "22182-3") %>%
        dplyr::select(1:3,19,20,21,22) %>%
        dplyr::mutate(orbit = as.numeric(as.character(OrbitRelative))) %>%
        dplyr::group_by(orbit) %>%
        dplyr::summarise_all(dplyr::first)
    # dplyr::filter(!orbit %in% c(109, 77, 66,23))

    sf::st_write(sub_fprints,
                 file.path(outdir, "AFRICA_S2_orbits.gpkg"), delete_dsn = TRUE)

    # Import S2 orbit tracks and subset over Africa ----
    #
    tracks <- sf::st_read("data-raw/S2A_rel_orbit_groundtrack_10Sec.kml") %>%
        sf::st_zm() %>%
        dplyr::filter(Relative_Orbit %in% unique(sub_fprints$orbit)) %>%
        dplyr::rename(orbit = Relative_Orbit) %>%
        dplyr::mutate(orbit = as.numeric(as.character(orbit))) %>%
        sf::st_intersection(sf::st_as_sfc(sf::st_bbox(africa) + c(-2,-2,+2,+2))) %>%
        sf::st_cast("LINESTRING")

    # identify horizontal lines ----
    diffs <- lapply(sf::st_geometry(tracks),
                    FUN = function(x) data.frame(diff = sf::st_bbox(x)$ymin - sf::st_bbox(x)$ymax)) %>%
        data.table::rbindlist()
    which_horiz <- which(diffs < -0.1)

    # save ----
    tracks <- tracks[which_horiz,] %>%
        dplyr::group_by(orbit) %>%
        dplyr::summarise(geometry = sf::st_union(geometry))

    sf::st_write(tracks,
                 file.path(outdir, "AFRICA_S2_tracks.gpkg"),
                 delete_dsn = TRUE)

    # sample points along tracks ----

    points <- sf::st_intersection(tracks, sf::st_graticule(lat = seq(-90, 90, lat_dist),
                                                           lon = c(-180,180))) %>%
        dplyr::select(orbit, degree, type) %>%
        dplyr::group_by(orbit) %>%
        dplyr::arrange(orbit) %>%
        dplyr::mutate(nn = dplyr::row_number()) %>%
        dplyr::mutate(id = paste(orbit, nn, sep = "_")) %>%
        dplyr::select(-nn)

    outfile <- file.path(outdir, paste0("AFRICA_S2_points_", lat_dist, "deg.gpkg"))
    sf::st_write(points, outfile, delete_dsn = TRUE)

    #Tiles ----

    tiles <- readRDS("data-raw/s2_tiles.rds")

    tiles_africa <- dplyr::filter(tiles, lengths(
        sf::st_intersects(tiles, sf::st_as_sfc(sf::st_bbox(africa) + c(-2,-2,+2,+2)))) > 0)
    tiles_africa <- tiles_africa %>%
        dplyr::mutate(utm_zone = substring(tile_id, 1,2))

    sf::st_write(tiles_africa, file.path(outdir,
                                         "AFRICA_S2_tiles.gpkg"),
                 delete_dsn = TRUE)

}
