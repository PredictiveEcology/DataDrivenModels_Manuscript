

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

