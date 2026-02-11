

landscapeSummary <- function(outputPath) {

  pgMapFiles <- list.files(outputPath, full.names = TRUE, pattern = "pixelGroupMap")
  cdFiles <- list.files(outputPath, full.names = TRUE, pattern = "cohortData")
  years <- stringr::str_extract(cdFiles, "year\\d+") |>
    na.omit()|>
    sub(pattern = "year", replacement = "") |>
    as.numeric()

  yearSummary <- lapply(years, function(year){
    cd <- readRDS(grep(cdFiles, pattern = year, value = TRUE))
    pg <- terra::rast(grep(pgMapFiles, pattern = year, value = TRUE))
    #summarize ANPP, mortality, and B by species
    pixelTable <- LandR::addNoPixel2CohortData(cd, pixelGroupMap = pg)
    pixelTable <- pixelTable[, .(sumANPP = sum(aNPPAct), sumB = sum(B), sumMortality = sum(mortality)), .(pixelGroup, speciesCode, noPixels)]
    npix <- sum(pixelTable[, .N, .(pixelGroup, noPixels)]$noPixels)
    #sum across pixels, divide by n pixels
    pixelTable <- pixelTable[, .(meanANPP = sum(sumANPP)/npix, 
                   meanB = sum(sumB)/npix,
                   meanMort = sum(sumMortality)/npix),
               .(speciesCode)]
    pixelTable[, propB := meanB/sum(meanB)]
    pixelTable[, year := year]
    pixelTable
  })
}

#utility function to get cohortDataLongs from file path 
# (optionally subset by some pixelIndex, given the resulting size)
getCDLongFromOutput <- function(outputsPath, source = NULL, pixelIndex = NULL) {
  CDs <- list.files(path = outputsPath, full.names = TRUE) |>
  lapply(FUN = makeCDLongs, pixelIndex = pixelIndex) |>
    lapply(rbindlist) |>
    rbindlist()
  
  if (!is.null(source)){
    CDs[, source := source]
  }
  return(CDs)
}

makeCDLongs <- function(outputPath, pixelIndex = NULL) {
  pgMapFiles <- list.files(outputPath, full.names = TRUE, pattern = "pixelGroupMap")
  cdFiles <- list.files(outputPath, full.names = TRUE, pattern = "cohortData")
  years <- stringr::str_extract(cdFiles, "year\\d+") |>
    na.omit()|>
    sub(pattern = "year", replacement = "") |>
    as.numeric()
  
  CDLong <- lapply(years, function(year, pix = pixelIndex) {
    cd <- readRDS(grep(cdFiles, pattern = year, value = TRUE))
    pg <- terra::rast(grep(pgMapFiles, pattern = year, value = TRUE))
    cdL <- LandR::addPixels2CohortData(cd, pg)
    if (!is.null(pix)) {
      cdL <- cdL[pixelIndex %in% pix,]
    }
    cdL[, year := year]
    return(cdL)
  })
  return(CDLong)
}

#' Extract species biomass from pixels based on subset

#' @param PixUnderThresh vector of pixels to examine
#' @param focalOutputPath the outputs directory for focal
#' @param singleOutputPath the outputs directory for single
#' @param year the simulation year to retrieve for focal and single cohortData
#'
#' @returns
#' @export
#'
#' @examples
pixelSummaryByThreshold <- function(PixUnderThresh, year, focalOutputPath, singleOutputPath) {
  pgFilename <- paste0("pixelGroupMap_year", year, ".tif")
  cdFilename <- paste0("cohortData_year", year, ".rds")
  
  PG_focal <- terra::rast(file.path(focalOutputPath, pgFilename))
  CD_focal <- readRDS(file.path(focalOutputPath, cdFilename))
  PG_single <- terra::rast(file.path(singleOutputPath, pgFilename))
  CD_single <- readRDS(file.path(singleOutputPath, cdFilename))

  subsetPGs <- as.vector(PG_focal)[PixUnderThresh]
  subsetPGs <- data.table(pixelGroup = subsetPGs)
  subsetPGs <- subsetPGs[, .N, .(pixelGroup)]
  
  
  subsetPGs <- subsetPGs[, .N, .(pixelGroup)] 
  
  #this will contain duplicates (potentially)
  CD_focal <- CD_focal[subsetPGs, on = c("pixelGroup")]

  #multiply B y N when summing
  #track the area too
  CD_focal <- CD_focal[, .(B = B * N), .(speciesCode, pixelGroup, ecoregionGroup, N)]
  CD_focal <- CD_focal[, .(B = sum(B), N = sum(N)), .(speciesCode, ecoregionGroup)]
  CD_focal[, source := "focal"]
  
  ##### single ####
  subsetPGs <- as.vector(PG_single)[PixUnderThresh]
  subsetPGs <- data.table(pixelGroup = subsetPGs)
  subsetPGs <- subsetPGs[, .N, .(pixelGroup)]
  #this will contain duplicates (potentially)
  CD_single <- CD_single[subsetPGs, on = c("pixelGroup")]
  #multiply B y N when summing
  
  CD_single <- CD_single[, .(B = B * N), .(speciesCode, pixelGroup, ecoregionGroup, N)]
  CD_single <- CD_single[, .(B = sum(B), N = sum(N)), .(speciesCode, ecoregionGroup)]
  CD_single[, source := "single"]
  
  # return
  return(list(single = CD_single, focal = CD_focal))
}

#' Find pixels satisfying initial biomass and competition req.

#' @param prop proportion of speciesEcoregion biomass below which pixels are taken 
#' @param minSpp minimum number of species in a pixelGroup (ie competition)
#' @param speciesEcoregion species ecoregion map
#' @param initialCD cohortData in the initial year
#' @param initialPG  pixelGroupMap in the initial year
getPixSubset <- function(initialCD, initialPG, speciesEcoregion,
                         prop = 0.25, minSpp = 2) {
  #setup
  initialCD <- initialCD[, sumB := sum(B), .(pixelGroup, ecoregionGroup)]
  ecoMaxB <- speciesEcoregion[, .(maxB = max(maxB)), .(ecoregionGroup)]
  initialCD <- initialCD[ecoMaxB, on = "ecoregionGroup"]
  initialCD[, propB := sumB/maxB]
  PGunderQnt <- unique(initialCD[propB < prop,]$pixelGroup) #35K unique groups
  initialCD[, nSpp := length(unique(speciesCode)), .(pixelGroup)]
  PGoverMinSpp <- unique(initialCD[nSpp > minSpp, ]$pixelGroup) #900K
  
  PixUnderThresh <- which(as.vector(initialPG) %in% PGunderQnt &
                            as.vector(initialPG) %in% PGoverMinSpp)
  
  #### focal ####
}

#` Extract the cohortData tables from biomass_yieldTables and return a subset

#the biomassYield pixelGroups are re-made from cohortData and will be identical between the two groups
# as Biomass_yield dissolves landscape into age and B = 0, and then grows unique pixelGroups 
#' @param pixelGroupSubset optional subset of pixelGroups
#' used to subset. N is the number of unique species in a pixelGroup.  
#' @param outputPath  location of the B_yT cohortDatas
getBiomassYieldCD <- function(outputPath, pixelGroupSubset = NULL) {

  #the biomassYield pixelGroups are re-made from cohortData and will be identical between the two groups
  # as Biomass_yield dissolves landscape into age and B = 0, and then grows unique pixelGroups 
  cohortDatas <- list.files(file.path(outputPath, "cohortDataYield"), full.names = TRUE)
  cohortCDs <- lapply(cohortDatas, function(x) {
    x <- qs2::qs_read(x)
    if (!is.null(pixelGroupSubset)) {
      x <- x[pixelGroup %in% pixelGroupSubset,]
    }
    return(x)
  }) |>
    rbindlist(use.names = TRUE)
  return(cohortCDs)
}
