#' @title extract_intensity
#' @description Compute total burnt area for each polygon of a tesselation, from
#'  a yearly fireCCI aggregated product ----
#' @param in_ba_file `character` Full path to a FireCCI product
#' @param in_tessels `character` Full path to a polygon tesselation file
#' @param out_file `character` Full path to file where results will be saved (if NULL,
#' results are only saved in the output object)
#' @return The function returns a `sf` object identical to `in_tessels`, with an
#'  additional `tot_ba` column, computed as the sum of Burnt Area of the cells
#'  of the raster intersecting each tessel, multiplied by the percentage of area included
#'  in the cell.
#'  If `out_file` is not NULL, results are saved also in the selected filename.
#' @details DETAILS
#' @examples
#' \dontrun{
#' extract_intensity(in_ba_file = "data-raw/2018-Annual-ESACCI-L4_FIRE-BA-MODIS-fv5.1_km2.tif",
#'                   In_tessels = "data-raw/L8_tessels.gpkg",
#'                   out_file = NULL
#' }
#' @rdname extract_intensity
#' @export
#' @author Lorenzo Busetto, phD (2019)
#' @importFrom sf st_read st_zm st_set_crs st_buffer st_intersects st_intersection
#'  st_as_sfc st_bbox st_write st_cast st_geometry st_union st_graticule
#' @importFrom dplyr group_by filter summarize select mutate summarise_all
#' first rename summarise arrange row_number
#' @importFrom data.table rbindlist
#'
extract_intensity <- function(in_ba_file,
                              in_tessels,
                              out_file = NULL) {

    in_ba   <- raster::raster(in_ba_file, quiet = TRUE)
    names(in_ba) <- "BA"
    in_tess <- sf::st_read(in_tessels, quiet = TRUE)
    in_tess <- dplyr::mutate(in_tess,
                             tot_ba = as.numeric(NA))
    for (hex in seq_len(dim(in_tess)[1])) {
        # for (hex in 1:100) {
        print(hex)
        test_ext <- raster::extract(in_ba,
                                    in_tess[hex,],
                                    weights = TRUE,
                                    normalizeWeights = FALSE,
                                    df = TRUE)
        out_val <- sum(test_ext[["BA"]] * test_ext[["weight"]], na.rm = TRUE)
        in_tess[hex,"tot_ba"] <- out_val
    }

    in_tess$tot_area <- as.numeric(sf::st_area(in_tess)) / 1000000
    # in_tess$frac_ba <- as.numeric(in_tess$tot_ba / area)
    in_tess <- in_tess %>%
        dplyr::mutate(frac_ba = tot_ba / tot_area) %>%
        dplyr::select(pathrow, biome_code, tot_ba, frac_ba, tot_area)

    if (!is.null(out_file)) {
        sf::st_write(in_tess, out_file, delete_layer = TRUE)
    }
    in_tess
}
