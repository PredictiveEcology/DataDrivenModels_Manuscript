repos <- c("predictiveecology.r-universe.dev", getOption("repos"))
# Need the latest version
if (tryCatch(packageVersion("SpaDES.project") < "0.1.1", error = function(x) TRUE)) {
  install.packages(c("SpaDES.project", "Require"), repos = repos)
}
#need cmake
projPath <- "~/git/DataDrivenModels_Manuscript"
customOpts <- list(gargle_oauth_email = "ianmseddy@gmail.com",
             gargle_oauth_cache = "~/google_drive_cache",
             gargle_oauth_client_type = "web",
             reproducible.inputPaths = "~/data"
             )
#usethis failed to install first time around...not sure why - 

studyAreaEcozone <- "Boreal Shield"
#can also be Montane Cordillera or Boreal PLains - note the typo


inSim <- SpaDES.project::setupProject(
  packages = c("usethis", "googledrive", "httr2", "RCurl", "XML", "bcdata"),
  useGit=  TRUE,
  paths = list(projectPath = projPath,
               modulePath = file.path("modules"),
               cachePath = file.path("cache"),
               scratchPath = tempdir(),
               inputPath = file.path("inputs"),
               outputPath = file.path("outputs")
  ),
  modules = c("PredictiveEcology/Biomass_speciesFactorial@development"
              # , "PredictiveEcology/Biomass_speciesParameters"
              # , "PredictiveEcology/Biomass_core"
              ),
  options = customOpts, 
  times = list(start = 2011, end = 2021),
  studyArea = {
    sa <- reproducible::prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip", 
                                   destinationPath = paths$inputPath, 
                                   fun = "terra::vect"
                                   )
    targetCRS <- terra::crs("EPSG:3348")
    # EPSG:3348 (NAD83(CSRS) / Statistics Canada Lambert) are commonly used for large areas of Canada. 
    # c("Clear Hills Upland") Boreal Plains - Ecoregion 137
    # c("Fraser Plateau") Montane Cordillera Ecoregion 202
    # c("Lake of the Woods") Boreal Shield Ecoregion 91
    ecodistrictName <- switch(studyAreaEcozone, 
                              "Boreal Shield" = {"Lake of the Woods"},
                              "Montane Cordillera" = {"Fraser Plateau"},
                              "Boreal PLains" = {"Clear Hills Upland"})
    sa <- sa[sa$REGION_NAM == ecodistrictName,] |>
      terra::project(targetCRS) 
    #consider buffering? not sure if we have disturbances yet
    #sa <- terra::buffer(sa, 5000)
    return(sa)
  },  
  rasterToMatch = {
    rtm <- terra::rast(sa, res = c(250, 250), vals = 1) |>
    reproducible::postProcess(maskTo = sa)
  },
  studyAreaANPP = {
    ecozones <- reproducible::prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip", 
                            destinationPath = paths$inputPath)
    ecozones <- ecozones[ecozones$ZONE_NAME == studyAreaEcozone]
  },
  sppEquiv = {
    spp <- LandR::speciesInStudyArea(studyArea = studyArea,
                                     dPath = "inputs",)
    sppEquiv <- LandR::sppEquivalencies_CA[KNN %in% spp$speciesList,]
    sppEquiv <- sppEquiv[LANDIS_traits != "",]
    sppEquiv <- sppEquiv[grep(pattern = "Spp", x = sppEquiv$KNN, invert = TRUE),]
  }
)
# 8: FADM - Timber Supply Area (TSA) (multiple, wms, kml, oracle_sde)
# ID: 8daa29da-d7f4-401c-83ae-d962e3a28980
