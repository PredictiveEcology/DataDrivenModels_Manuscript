repos <- c("predictiveecology.r-universe.dev", getOption("repos"))
# Need the latest version
if (tryCatch(packageVersion("SpaDES.project") < "0.1.1", error = function(x) TRUE)) {
  install.packages(c("SpaDES.project", "Require"), repos = repos)
  install.packages(c("reproducible", "SpaDES.core"), repos = repos)
}
#needed cmake

#### Settings users will want to set themselves ####

#added googledrive json to ~/.secrets for non-interactive use



##### Set up ####
studyAreaEcozone <- "Montane Cordillera"
studyAreaName <- "Fraser_Basin_district980"
#studyAreaName Should be BowronValley_986 but I forgot to update this
#

projPath <- getwd()

simProject <- SpaDES.project::setupProject(
  packages = c("usethis", "googledrive", "httr2", "RCurl", "XML", "bcdata"),
  useGit = FALSE, #tried this to keep it from updating 
  require = c("reproducible",
              "PredictiveEcology/SpaDES.experiment@development (>= 0.0.2.9005)",
              "PredictiveEcology/LandR@development (>= 1.1.5.9091)"),
  options = list(reproducible.inputPaths = "~/data", 
                 gargle_oauth_email = "ianmseddy@gmail.com", #TODO: use a file (this is here from debugging tmux)
                 gargle_oauth_cache = ".secrets",
                 reproducible.useMemoise = FALSE, ##limit RAM use due to immense size of factorial object
                 terra.memfrac = 0, #limit RAM use due to immense size of factorial object
                 gargle_oauth_client_type = "web"),
  paths = list(projectPath = projPath,
               modulePath = file.path("modules"),
               cachePath = file.path("cache"),
               scratchPath = tempdir(),
               inputPath = file.path("inputs"),
               outputPath = file.path("outputs")),
  modules = c(
    "PredictiveEcology/Biomass_speciesFactorial@development"
    , "PredictiveEcology/Biomass_borealDataPrep@development"
    , "PredictiveEcology/Biomass_speciesParameters@development"
    , "PredictiveEcology/Biomass_core@development"
  ),
  times = list(start = 2020, end = 2120),
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
    subzone$fullZoneSbZone <- paste(subzone$ZONE, subzone$SUBZONE)
    subzone$subzoneNum <- as.numeric(as.factor(subzone$fullZoneSbZone))
    return(subzone)
  },
  sppEquiv = {
    species <- LandR::speciesInStudyArea(studyArea = studyArea, dPath = paths$inputPath, sppEquivCol = "LandR")
    sppEquiv <- LandR::sppEquivalencies_CA[LandR %in% species$speciesList,]
    sppEquiv <- sppEquiv[LANDIS_traits != "",]
    sppEquiv[LandR == "Popu_bal", LandR := "Popu_tre"]
  },
  sppEquivLong = {
    #this object dictates the biomass equation to determine AGB from DBH + Height via the PSP column
    #TODO: this will eventually use a dedicated biomass equation column 
    sppEquivLong <- LandR::sppEquivalencies_CA
    sppEquivLong[LandR == "Pice_eng_gla", PSP := "engelmann spruce"]
  },
  argsForFactorial = {
    a <- list(cohortsPerPixel = 1:2,
              growthcurve = seq(0, 1, 0.05),
              mortalityshape = seq(21, 25, 2),
              longevity = seq(150, 550, 50),
              mANPPproportion = seq(2.8, 6.4, 0.4))
  },
  params = list(
    .globals = list(
      dataYear = 2020,
      .plots = "png",
      .studyAreaName = studyAreaName,
      initialB = 20, #this is necessary if including growthcurves > 0.9
      minCohortBiomass = 9 
    ), 
    Biomass_speciesParameters = list(
      standAgesForFitting = c(20, 125), 
      PSPdataTypes = c("BC", "NFI", "AB"),
      minDBH = 5,
      quantileAgeSubset = 99
    ), 
    Biomass_borealDataPrep = list(
      ecoregionLayerField = "subzoneNum"
    )
  ), 
  outputs =  {
    outputs <- rbind(
      data.table(objectName = "pixelGroupMap", saveTime = c(seq(times$start, times$end, 5)), 
                 exts = ".tif", fun = "writeRaster", package = "terra"), 
      data.table(objectName = "cohortData", saveTime = c(seq(times$start, times$end, 5))), 
      data.table(objectName = "speciesEcoregion", saveTime = times$end), 
      data.table(objectName = "ecoregion", saveTime = times$end), 
      data.table(objectName = "species", saveTime = times$end),
      data.table(objectName = "ecoregionMap", saveTime = times$end, exts = ".tif", 
                 fun = "writeRaster", package = "terra"),
      fill = TRUE
    )
    outputs[is.na(fun), c("exts", "fun", "package") := .("rds", "saveRDS", "base")]
    outputs <- as.data.frame(outputs)
    # outputs$arguments <- list("writeRaster" = list("overwrite" = TRUE))
    outputs$arguments <- list(overwrite = TRUE)
    return(outputs)
  }
)

#####experiment args #### 
inSim <- do.call(simInit, simProject)
#save the static objects
terra::writeRaster(inSim$rstLCC, "outputs/rstLCC.tif", overwrite = TRUE)
#originally disturbed forest landcover (LCC 240) are reclassified in ecoregionRst  
ecoregionKey <- as.data.table(inSim$ecoregionLayer)
ecoregionKey <- unique(ecoregionKey[,.(ZONE, SUBZONE, subzoneNum)])
saveRDS(ecoregionKey, "outputs/ecoregionKey.rds") 
#join this with cats(sim$ecoregionMap) by subzoneNum = ecoregionName to get the zone and subzone of all ecoregions

#TODO: we do not need 3 replicates from BSP
focalSims <- SpaDES.experiment::experiment(inSim, replicates = 3, dirPrefix = "focalFitting_MC")
#run single
inSim@params$Biomass_speciesParameters$speciesFittingApproach <- "single"
singleSims <- SpaDES.experiment::experiment(inSim, replicates = 3, dirPrefix = "singleFitting_MC")

