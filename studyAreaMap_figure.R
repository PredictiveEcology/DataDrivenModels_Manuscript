library(Require)
#this figure was added last and separately
Require::Require("reproducible")
Require("sf")
Require("terra")
Require("ggplot2")
Require("ggspatial")
Require("tmap")
Require("tidyterra")
Require("scales")
Require("PredictiveEcology/climateData")
Require("shadowtext")
Require("grid")
##version 1 ####
data("World", package = "tmap")

options("reproducible.cachePath" = "cache")

#get data
# 
# sa <- reproducible::prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip", 
#                                destinationPath = "inputs", 
#                                fun = "terra::vect"
# )
# 
# sa$REGION_NOM <- NULL #ditch french -it causes the invalid multibyte
# 

# ---- Projection ----
targetCRS <- st_crs("EPSG:3348")  # Canada Albers Equal Area

# ---- Study area ----
sa <- reproducible::prepInputs(
  url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip",
  destinationPath = "inputs",
  fun = "terra::vect"
)
sa <- sa[sa$ECOREGION %in% c(137, 138),] 

# ---- Political boundaries ----
na <- prepInputs(
  url = paste0(
    "https://www.cec.org/files/atlas_layers/0_reference/",
    "0_01_political_boundaries/politicalboundaries_shapefile.zip"
  ),
  destinationPath = "inputs"
)

ca <- na[na$COUNTRY == "CAN",]
# ---- Define a Canada-centered bbox
can_bbox <- st_bbox(ca)

na <- na[na$COUNTRY %in% c("USA", "CAN"),]
na <- na[!na$NAME_En %in% c("Hawaii", "Puerto Rico"),]

na <- st_crop(na, can_bbox)

sa <- st_as_sf(sa) |>
  st_transform(st_crs(na))
sa_bbox <- st_bbox(sa)
sa_inset <- postProcess(na, cropTo = sa) |>
  sf::st_union()

# ---- Simple graticules (cheap & correct) ----
graticules <- st_graticule(
  lat = seq(40, 80, by = 10),
  lon = seq(-160, -40, by = 20),
  crs = 4326
) |>
  st_transform(crs(na)) |>
  st_crop(can_bbox)


# ---- Plot ----
insetMap <- ggplot() +
  geom_sf(data = na,
          fill = "grey93", 
          color = "grey70", 
          linewidth = 0.2) +
  geom_sf(data = graticules, color = "grey85", 
          linewidth = 0.3, linetype = "dotted") +
  geom_sf(data = sa_inset, color = "black", alpha = 0.8, linewidth = 0.4) +
  coord_sf(xlim = c(can_bbox["xmin"], can_bbox["xmax"]),
           ylim = c(can_bbox["ymin"], can_bbox["ymax"]),
           expand = FALSE) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    plot.background = element_rect(fill = "white", colour = NA)
)

insetMap
# ggsave(
#   plot = insetMap,
#   filename = "manuscript_figures/studyAreaMap.png",
#   width  = 170,
#   height = 130,
#   units  = "mm"
# )

#inset map could use some text saying "Canada"

# 1. Ensure everything is in Web Mercator -------------------------------

# sa_3857       <- st_transform(sa, 3857)
# sal_bbox_3857 <- st_bbox(st_transform(sal, 3857))

# 2. Build map ------------------------------------------------------------
targetCRS_WebMercator <- terra::crs("EPSG:3857") 
sa <- st_transform(sa, targetCRS_WebMercator)

sal <- st_buffer(sa, 200000)
sal_bbox <- st_bbox(sal)

mainMap <- ggplot() +
  
  coord_sf(
    crs = 3857,
    xlim = sal_bbox[c("xmin", "xmax")],
    ylim = sal_bbox[c("ymin", "ymax")],
    expand = FALSE
  ) +
  
  # Basemap tiles
  basemaps::basemap_gglayer(
    sal_bbox,
    map_service = "carto",
    map_type = "voyager_no_labels"
    ) +
  scale_fill_identity(guide = "none") + 
  # Study area outline
  geom_sf(
    data = sa,
    fill = NA,
    alpha = 0.6,
    color = "grey48",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_shadowtext(
    data = sa,
    inherit.aes = FALSE,
    aes(label = REGION_NAM, 
        geometry = geometry),
    stat = "sf_coordinates",
    colour = "darkgreen",
    bg.colour = "white",
    bg.r = 0.15,
    size = 4,

  ) +
  # Scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.25
  ) +
  theme_classic() +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )


#make the inset

bb  <- sal_bbox   # any sf object already used in mainMap
xspan <- bb["xmax"] - bb["xmin"]
yspan <- bb["ymax"] - bb["ymin"]

inset_xmin <- bb["xmin"] + 0.01 * xspan
inset_xmax <- bb["xmin"] + 0.35 * xspan

inset_ymin <- bb["ymax"] - 0.35 * yspan
inset_ymax <- bb["ymax"] - 0.01 * yspan



finalMap <- mainMap +
  annotation_custom(
    grob = ggplotGrob(insetMap),
    xmin = inset_xmin,
    xmax = inset_xmax,
    ymin = inset_ymin,
    ymax = inset_ymax
  )

finalMap

ggsave(
  plot = mainMap,
  filename = "manuscript_figures/studyAreaMap.png",
  width  = 170,
  height = 130,
  units  = "mm"
)
