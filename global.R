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
studyAreaEcozone <- c("Montane Cordillera", "Boreal PLain") #note this capital L in PLains is necessary!
#studyAreaName Should be BowronValley_986 

studyAreaName = "ClearHillsUpland137_PeaceLowland138"
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
               outputPath = file.path("outputs/focal")),
  modules = c(
    "PredictiveEcology/Biomass_speciesFactorial@development"
    , "PredictiveEcology/Biomass_borealDataPrep@development"
    , "PredictiveEcology/Biomass_speciesParameters@development"
    , "PredictiveEcology/Biomass_core@development"
    , "PredictiveEcology/Biomass_yieldTables@main"
  ),
  times = list(start = 2020, end = 2120),
  studyArea = {
    sa <- reproducible::prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip", 
                                   destinationPath = paths$inputPath, 
                                   fun = "terra::vect"
    )
    targetCRS <- terra::crs("EPSG:3348")
    # EPSG:3348 (NAD83(CSRS) / Statistics Canada Lambert) are commonly used for large areas of Canada. 
    sa <- sa[sa$ECOREGION %in% c(137,138),] |> 
      terra::project(targetCRS) 
    #980 is in the Fraser Basin ecoregion - it conveniently has no Douglas-fir
    #buffer it for dispersal
    sa <- terra::aggregate(sa)
    return(sa)
  },  
  rasterToMatch = {
    rtm <- terra::rast(sa, res = c(250, 250), vals = 1) |>
      reproducible::postProcess(maskTo = sa)
  },
  studyAreaANPP = {
    ecozones <- reproducible::prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
                                         destinationPath = paths$inputPath)
    ecozones <- ecozones[ecozones$ZONE_NAME %in% studyAreaEcozone,]
  },
  sppEquiv = {
    species <- LandR::speciesInStudyArea(studyArea = studyArea, dPath = paths$inputPath, sppEquivCol = "LandR")
    sppEquiv <- LandR::sppEquivalencies_CA[LandR %in% species$speciesList,]
    sppEquiv <- sppEquiv[LANDIS_traits != ""]
    #this combines engelmann with hybrid spruce for simplicity
    sppEquiv[LandR == "Pice_eng_gla", LandR := "Pice_eng"]
    sppEquiv <- sppEquiv[!LandR %in% c("Pinu_ban", "Betu_pap", "Lari_lar", "Pinu_ban")]
    #the initial landscape is ~ 42% aspen, 17% black spruce, 18% white spruce, 
    # 12% pine, 4% engelmann, 2.5% balsam poplar, and 2.4% larch, 1.4% birch, 
    #0.5% jack pine and 0.01% subalpine fir. So remove those under 2.5%
    sppEquiv
  },
  sppEquivLong = {
    #this object dictates the biomass equation to determine AGB from DBH + Height via the PSP column
    sppEquivLong <- LandR::sppEquivalencies_CA
    #this gives hybrid spruce the equation of engelmann 
    sppEquivLong[Latin_full == "Picea engelmannii x glauca", c("PSP") := .("engelmann spruce")]
    sppEquivLong
  },
  argsForFactorial = {
    a <- list(cohortsPerPixel = 1:2,
              growthcurve = seq(0, 0.95, 0.05), #19
              mortalityshape = seq(21, 25, 2),#3
              longevity = seq(150, 500, 50), #7
              mANPPproportion = seq(2.8, 6.4, 0.4)) #9
  },
  params = list(
    .globals = list(
      dataYear = 2020,
      .plots = "png",
      .studyAreaName = studyAreaName,
      initialB = 20, #this is necessary if including growthcurves > 0.9
      #TODO: actually there are still some cohorts that remain stuck at 20 with growthcurve 0.9 
      minCohortBiomass = 1 
    ), 
    Biomass_speciesParameters = list(
      standAgesForFitting = c(20, 125), 
      PSPdataTypes = c("BC", "NFI", "AB", "SK"),
      minDBH = 7 #minDBH is 9 in SK 
      #but I don't think the disparity matters much the way sppParams uses data 
      #and a lower DBH means more young plots are available 
    ), 
    Biomass_core = list(
      seedingAlgorithm = "noSeeding"
    )
  ), 
  outputs =  {
    outputs <- rbind(
      data.table(objectName = "pixelGroupMap", saveTime = c(seq(times$start, times$end, 10)), 
                 exts = ".tif", fun = "writeRaster", package = "terra"), 
      data.table(objectName = "cohortData", saveTime = c(seq(times$start, times$end, 10))), 
      data.table(objectName = "speciesEcoregion", saveTime = times$end), 
      data.table(objectName = "ecoregion", saveTime = times$end), 
      data.table(objectName = "species", saveTime = times$end),
      data.table(objectName = "ecoregionMap", saveTime = times$end, exts = ".tif", 
                 fun = "writeRaster", package = "terra"),
      data.table(objectName =  "simulatedBiomassMap", saveTime = times$start, 
                 exts = ".tif", fun = "writeRaster", package = "terra"),
      fill = TRUE
    )
    outputs[is.na(fun), c("exts", "fun", "package") := .("rds", "saveRDS", "base")]
    outputs <- as.data.frame(outputs)
    outputs$arguments <- list(overwrite = TRUE)
    return(outputs)
  }
)

focalSim <- simInitAndSpades2(simProject)

#run single

simProject$params$Biomass_speciesParameters$speciesFittingApproach <- "single"
simProject$paths$outputPath <- "outputs/single"
singleSim <- simInitAndSpades2(simProject)

# analyis is only performed on files in outputs, 
# so the simSingle and focalSim objects won't actually be used (the simList is not saved)
