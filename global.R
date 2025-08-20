repos <- c("predictiveecology.r-universe.dev", getOption("repos"))
# Need the latest version
if (tryCatch(packageVersion("SpaDES.project") < "0.1.1", error = function(x) TRUE)) {
  install.packages(c("SpaDES.project", "Require"), repos = repos)
}
#needed cmake

#### Settings users will want to set themselves ####

projPath <- "~/git/DataDrivenModels_Manuscript"

options(gargle_oauth_email = "ianmseddy@gmail.com")

if (grep("Windows", osVersion)) {
  options(gargle_oauth_client_type = "installed", 
          gargle_oauth_cache = "../../google_drive_cache")
} else {
  options(gargle_oauth_cache = "~/google_drive_cache",
          gargle_oauth_client_type = "web")
}

##### Set up ####
studyAreaEcozone <- "Montane Cordillera"
studyAreaName <- "Fraser_Basin_district980"

projPath <- getwd()
simProject <- SpaDES.project::setupProject(
  packages = c("usethis", "googledrive", "httr2", "RCurl", "XML", "bcdata"),
  useGit = TRUE,
  require = c("PredictiveEcology/reproducible@AI (>= 2.1.2.9050)",
              "PredictiveEcology/SpaDES.core@box (>= 2.1.5.9022)",
              "PredictiveEcology/SpaDES.experiment (>= 0.0.2.9005)"),
  paths = list(projectPath = projPath,
               modulePath = file.path("modules"),
               cachePath = file.path("cache"),
               scratchPath = tempdir(),
               inputPath = file.path("inputs"),
               outputPath = file.path("outputs")),
  modules = c(
    "PredictiveEcology/Biomass_speciesFactorial@development"
    , "PredictiveEcology/Biomass_borealDataPrep@development"
    , "PredictiveEcology/Biomass_speciesParameters@fixingSingle"
    , "PredictiveEcology/Biomass_core@development"
  ),
  times = list(start = 2011, end = 2091),
  studyArea = {
    sa <- reproducible::prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip", 
                                   destinationPath = paths$inputPath, 
                                   fun = "terra::vect"
    )
    targetCRS <- terra::crs("EPSG:3348")
    # EPSG:3348 (NAD83(CSRS) / Statistics Canada Lambert) are commonly used for large areas of Canada. 
    
    sa <- sa[sa$ECODISTRIC == 986,] |> #980 is like 90% pine
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
  studyAreaANPP = {
    ecozones <- reproducible::prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
                                         destinationPath = paths$inputPath)
    ecozones <- ecozones[ecozones$ZONE_NAME == studyAreaEcozone,]
  },
  ecoregionLayer = {
    subzone <- Cache(bcdata::bcdc_get_data, 
                     record = "f358a53b-ffde-4830-a325-a5a03ff672c3", 
                     userTags = c("BECs")) |>
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
              growthcurve = seq(0.1, 0.9, 0.05),
              mortalityshape = seq(21, 25, 2), 
              longevity = seq(150,550, 50), 
              mANPPproportion = seq(3.3, 6.6, 0.3)) 
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
      quantileAgeSubset = 99, 
      .useCache = ".inputObjects"
      
    ), 
    Biomass_core = list(
      .useCache = ".inputObjects"
    ),
    Biomass_borealDataPrep = list(
      ecoregionLayerField = "SUBZONE"
    )
  ), 
  outputs =  {
    outputs <- rbind(
      data.table(objectName = "pixelGroupMap", 
                 saveTime = c(seq(times$start, times$end, 5)), 
                 exts = ".tif", fun = "writeRaster", package = "raster"), 
      data.table(objectName = "cohortData", 
                 saveTime = c(seq(times$start, times$end, 5))), 
      data.table(objectName = "speciesEcoregion", 
                 saveTime = times$end), 
      data.table(objectName = "species", 
                 saveTime = times$end),
      fill = TRUE
    )
    outputs[is.na(fun), c("exts", "fun", "package") := .("rds", "saveRDS", "base")]
    outputs <- as.data.frame(outputs)
    # outputs$arguments <- list("writeRaster" = list("overwrite" = TRUE))
    outputs
  }
)

#####experiment args #### 
inSim <- do.call(simInit, simProject)

SpaDES.experiment::experiment(inSim, replicates = 3, dirPrefix = "focalFitting_MC")


inSim$params$Biomass_speciesParameters$speciesFittingApproach <- "single"

SpaDES.experiment::experiment(inSim, replicates = 3, dirPrefix = "singleFitting_MC")


#
simProject$params$Biomass_speciesParameters$speciesFittingApproach <- "focal"


# cdRep1 <- readRDS("outputs/focalFitting_all/rep1/cohortData_year2051.rds")
# cdRep2 <- readRDS("outputs/focalFitting_all/rep2/cohortData_year2051.rds") 
# pgRep1 <- rast("outputs/focalFitting_all/rep1/pixelGroupMap_year2051.tif")
# pgRep2 <- rast("outputs/focalFitting_all/rep2/pixelGroupMap_year2051.tif")
# cdRep1[pixelGroup %in% pgRep1[30000],]
# cdRep2[pixelGroup %in% pgRep2[30000],]
