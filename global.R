repos <- c("predictiveecology.r-universe.dev", getOption("repos"))
# Need the latest version
if (tryCatch(packageVersion("SpaDES.project") < "0.1.1", error = function(x) TRUE)) {
  install.packages(c("SpaDES.project", "Require"), repos = repos)
}
#need cmake

#### Settings users will want to set themselves ####
projPath <- "~/git/DataDrivenModels_Manuscript"
customOpts <- list(gargle_oauth_email = "ianmseddy@gmail.com",
             gargle_oauth_cache = "~/google_drive_cache",
             gargle_oauth_client_type = "web"
             )
 
##### Set up ####
studyAreaEcozone <- "Montane Cordillera"
studyAreaName <- "Fraser_Basin_district980"


simFocal <- SpaDES.project::setupProject(
  packages = c("usethis", "googledrive", "httr2", "RCurl", "XML", "bcdata"),
  useGit=  TRUE,
  require = c("PredictiveEcology/reproducible@AI (>= 2.1.2.9050)", 
              "PredictiveEcology/SpaDES.core@box (>= 2.1.5.9022)", 
              "PredictiveEcology/SpaDES.experiment (>= 0.0.2.9005)"),
  paths = list(projectPath = projPath,
               modulePath = file.path("modules"),
               cachePath = file.path("cache"),
               scratchPath = tempdir(),
               inputPath = file.path("inputs"),
               outputPath = file.path("outputs")),
  modules = c("PredictiveEcology/Biomass_speciesFactorial@development"
              , "PredictiveEcology/Biomass_borealDataPrep@development"
              , "PredictiveEcology/Biomass_speciesParameters@development"
              , "PredictiveEcology/Biomass_core@development"
  ),
  options = c(customOpts, 
              'spades.saveFileExtensions' = data.frame(fun = c("saveRDS", "writeRaster"), 
                                                       package = c("base", "terra"), 
                                                       exts = c("rds", "tif"))),
  times = list(start = 2011, end = 2031),
  studyArea = {
    sa <- reproducible::prepInputs(url = "hhttps://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip", 
                                   destinationPath = paths$inputPath, 
                                   fun = "terra::vect"
    )
    targetCRS <- terra::crs("EPSG:3348")
    # EPSG:3348 (NAD83(CSRS) / Statistics Canada Lambert) are commonly used for large areas of Canada. 
    
    sa <- sa[sa$ECODISTRIC == 980,] |>
      terra::project(targetCRS) 
    #980 is in the Fraser Basin ecoregion - it conveniently has no Douglas-fir
    return(sa)
  },  
  rasterToMatch = {
    rtm <- terra::rast(sa, res = c(250, 250), vals = 1) |>
      reproducible::postProcess(maskTo = sa)
  },
  studyAreaANPP = {
    ecozones <- reproducible::prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip", 
                                         destinationPath = paths$inputPath)
    ecozones <- ecozones[ecozones$ZONE_NAME == studyAreaEcozone,]
  },
  sppEquiv = {
    spp <- LandR::speciesInStudyArea(studyArea = studyArea,
                                     dPath = "inputs",)
    sppEquiv <- LandR::sppEquivalencies_CA[KNN %in% spp$speciesList,]
    sppEquiv <- sppEquiv[LANDIS_traits != "",]
    sppEquiv <- sppEquiv[grep(pattern = "Spp", x = sppEquiv$KNN, invert = TRUE),]
  },
  argsForFactorial = {
    a <- list(cohortsPerPixel = 1:2,
              growthcurve = seq(0.65, 0.85, 0.02),
              mortalityshape = seq(19, 25, 2), 
              longevity = seq(150, 600, 50), 
              mANPPproportion = seq(3.3, 6.3, 0.3)) 
  },
  params = list(
    .globals = list(
      dataYear = 2011,
      .plots = "png",
      .studyAreaName = studyAreaName,
      minCohortBiomass = 1 #for plotting purposes - so many growth curves lose their tail if we remove B < 10
    ), 
    Biomass_speciesParameters = list(
      standAgesForFitting = c(21, 121), 
      quantileAgeSubset = 99
      
    )
  ), 
  outputs = rbind(
      data.frame(objectName = "pixelGroupMap", fun = "terra::writeRaster", 
                 saveTime = c(seq(times$start, times$end, 5)), 
                 arguments = list(overwrite = TRUE), ext = "tif"), 
      data.frame(objectName = "cohortData", fun = "saveRDS", 
                 saveTime = c(seq(times$start, times$end, 5)), 
                 arguments = list(overwrite = TRUE), ext = "tif"),
      data.frame(objectName = c("species"), 
                 fun = "saveRDS", saveTime = times$end, 
                 arguments = list(overwrite = TRUE), ext = "rds"), 
      data.frame(objectName = c("speciesEcoregion"), 
                 fun = "saveRDS", saveTime = times$end, 
                 arguments = list(overwrite = TRUE), ext = "rds")
    )
)

#experiment args - run name, outputs, 


simFocalInit <- do.call(simInit, simFocal)



