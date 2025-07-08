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


simProject <- SpaDES.project::setupProject(
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
  times = list(start = 2011, end = 2091),
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
    #buffer it for dispersal
    sa <- terra::buffer(sa, 5000)
    return(sa)
  },  
  rasterToMatch = {
    rtm <- terra::rast(sa, res = c(250, 250), vals = 1) |>
      reproducible::postProcess(maskTo = sa)
  },
  # studyAreaANPP = {
  #   ecozones <- reproducible::prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip", 
  #                                        destinationPath = paths$inputPath)
  #   ecozones <- ecozones[ecozones$ZONE_NAME == studyAreaEcozone,]
  # },
  ecoregionLayer = {
    subzone <- bcdata::bcdc_get_data("f358a53b-ffde-4830-a325-a5a03ff672c3") |>
      reproducible::postProcess(to = studyArea)
    return(subzone)
  },
  sppEquiv = {
    spp <- LandR::speciesInStudyArea(studyArea = studyArea,
                                     dPath = "inputs",)
    sppEquiv <- LandR::sppEquivalencies_CA[KNN %in% spp$speciesList,]
    sppEquiv <- sppEquiv[LANDIS_traits != "",]
    sppEquiv <- sppEquiv[grep(pattern = "Spp", x = sppEquiv$KNN, invert = TRUE),]
    sppEquiv[LandR == "Popu_bal", LandR := "Popu_tre"]
    #there is insignificant balsam poplar in the montane cordillera, so combine it with aspen
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
      .useCache = ".inputObjects",
      dataYear = 2011,
      .plots = "png",
      .studyAreaName = studyAreaName,
      minCohortBiomass = 1 #for plotting purposes - so many growth curves lose their tail if we remove B < 10
    ), 
    Biomass_speciesParameters = list(
      standAgesForFitting = c(21, 121), 
      quantileAgeSubset = 99
    ), 
    Biomass_core = list(
      useCache = ".inputObjects"
    ), 
    Biomass_borealDataPrep = list(
      ecoregionLayerField = "SUBZONE"
    )
  ), 
  outputs =  {
    outputs <- rbind(
      data.table(objectName = "pixelGroupMap", 
                 saveTime = c(seq(times$start, times$end, 5)), 
                 exts = ".tif", fun = "writeRaster", package = "terra"), 
      data.table(objectName = "cohortData", 
                 saveTime = c(seq(times$start, times$end, 5))), 
      data.table(objectName = "speciesEcoregion", 
                 saveTime = times$end), 
      data.table(objectName = "species", 
                 saveTime = times$end),
      fill = TRUE
    )
    outputs[is.na(fun), c("exts", "fun", "package") := .("rds", "saveRDS", "base")]
    outputs
  }
)

#####experiment args #### 
inSim <- do.call(what = SpaDES.core::simInit, simProject)

SpaDES.experiment::experiment(inSim, replicates = 3, dirPrefix = "focalFitting_all")

inSim@params$Biomass_speciesParameters$speciesFittingApproach <- "single"

SpaDES.experiment::experiment(inSim, replicates = 3, dirPrefix = "singleFitting_all")


