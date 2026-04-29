library(Require)
#this figure was added last and separately
Require::Require("reproducible")
Require(sf)
Require(ggplot2)
Require(ggspatial)
Require("tmap")

##version 1 ####
data("World", package = "tmap")


sa <- reproducible::prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip", 
                               destinationPath = "inputs", 
                               fun = "terra::vect"
)
sa$REGION_NOM <- NULL #ditch french -it causes the invalid multibyte

targetCRS <- terra::crs("EPSG:3348")# (NAD83(CSRS) / Statistics Canada Lambert) are commonly used for large areas of Canada. 
# EPSG:3348 (NAD83(CSRS) / Statistics Canada Lambert) are commonly used for large areas of Canada. 
sa <- sa[sa$ECOREGION %in% c(137,138),] |> 
  terra::project(targetCRS)  |>
  sf::st_as_sf()

bbox <- st_bbox(sa)
Canada  <- prepInputs(url = 'https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000b21a_e.zip', 
                      destinationPath = "inputs") |>
  st_transform(targetCRS)
can_bbox <- st_bbox(Canada)
# 
# na <- World[World$name %in% c("United States of America", "Canada"),] |>
#   st_transform(targetCRS)
# na <- st_crop(na, can_bbox)

na <- prepInputs(url = paste0("https://www.cec.org/files/atlas_layers/0_reference/",
                              "0_01_political_boundaries/politicalboundaries_shapefile.zip"),
                 destinationPath = "inputs", 
                 projectTo = targetCRS)
na <- st_crop(na, can_bbox)
# Build graticules in lon/lat, then project
graticules <- st_graticule(
  lon = seq(-160, -40, by = 10),
  lat = seq(40,  80,  by = 10),
  crs = st_crs(4326)
) |>
  st_transform(targetCRS) |>
  st_crop(can_bbox)



# 3. Build the figure
MapFigure <- ggplot() +
  scale_fill_brewer(
    palette = "Set2",
    name = "Ecoregion"
  ) +
  geom_sf(
    data = na,
    fill = "grey93",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_sf(
    data = graticules,
    color = "grey85",
    linewidth = 0.3,
    linetype = "dotted"
  ) +
  # Study areas
  geom_sf(
    aes(fill = REGION_NAM),
    data = sa,
    color = "black",
    linewidth = 0.4,
    alpha = 0.85
  ) +
  annotation_scale(location = "bl", width_hint = 0.25, 
                   bar_cols = c("darkgrey", "grey95")) +
  annotation_north_arrow(
    location = "tl",
    height = unit(1.2, "cm"),
    width  = unit(0.7, "cm"),
    style = north_arrow_orienteering(fill = c("white", "darkgrey"))
  ) +
  coord_sf(crs = targetCRS, expand = FALSE) +
  theme_classic() +
  theme(
    # axis.text   = element_blank(),
    # axis.ticks  = element_blank(),
    # axis.title  = element_blank(),
    legend.position = c(0.82, 0.78),
    legend.background = element_rect(
      fill = alpha("white", 0.8),
      color = NA
    )
)

ggsave(
  plot = MapFigure,
  filename = "manuscript_figures/studyAreaMap.png",
  width  = 170,
  height = 130,
  units  = "mm"
)

