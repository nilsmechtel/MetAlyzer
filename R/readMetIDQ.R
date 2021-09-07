## read metabolite names and classes, raw data, meta data and flag data ##
readMetIDQ <- function(file, verbose=TRUE) {
  if (verbose) {
    print("------------------------------")
    print("Loading full sheet...")
  }

  ## read the full Excel sheet ##
  fullSheet <- openxlsx::read.xlsx(file,
                                   sheet = 1,
                                   colNames = FALSE,
                                   skipEmptyRows=FALSE,
                                   skipEmptyCols = FALSE)
  fullSheet <- as.matrix(fullSheet)
  fullSheet[fullSheet == "NA"] <- NA # all entries are strings

  ## data range ##
  colClass <- ceiling(which(fullSheet == "Class") / nrow(fullSheet)) # column of header "Class"
  colsData <- (colClass+1):ncol(fullSheet) # columns

  colSampleType <- ceiling(which(fullSheet == "Sample Type") / nrow(fullSheet)) # column of "Sample Type"
  rowsData <- which(fullSheet[, colSampleType] == "Sample") # rows

  ## metabolites with classes ##
  rowClass <- which(fullSheet == "Class") %% nrow(fullSheet) # row of header "Class"
  metabolites <- fullSheet[rowClass-1, colsData] # metabolites are a row above classes
  classes <- fullSheet[rowClass, colsData]
  names(metabolites) <- classes

  if (verbose) {
    print(paste("Number of metabolites:", length(metabolites)))
    print(paste("Number of classes:", length(unique(names(metabolites)))))

    print("------------------------------")
    print("Reading data...")
  }

  ## raw data ##
  rawData <- as.data.frame(fullSheet[rowsData, colsData])
  colnames(rawData) <- metabolites
  rawData[] <- sapply(rawData[], as.numeric)

  ## meta data ##
  rowSampleType <- which(fullSheet == "Sample Type") %% nrow(fullSheet) # row of "Sample Type"
  metaHeader <- paste(fullSheet[rowSampleType, 1:colClass])
  metaData <- as.data.frame(fullSheet[rowsData, 1:colClass])
  colnames(metaData) <- metaHeader
  metaData <- metaData[colSums(is.na(metaData)) != nrow(metaData)] # remove full NA columns
  nas <- grep("NA", colnames(metaData))
  colnames(metaData)[nas] <- paste0("X", seq_along(nas)) # replace NAs in colnames with X and consecutive numbers

  if (verbose) {
    print(paste("Number of samples:", nrow(metaData)))
    print(paste("Columns in meta data:", ncol(metaData)))
  }

  ## background color ##
  if (verbose) {
    print("------------------------------")
    print("Reading background color...")
  }

  wb <- loadWorkbook(file)
  sheet1 <- getSheets(wb)[[1]]
  rows <- getRows(sheet1)
  cells <- getCells(rows)
  styles <- sapply(cells, getCellStyle) # get style of each cell
  bg <- sapply(styles, function(style) { # get background color of each cell
    fg  <- style$getFillForegroundXSSFColor()
    rgb <- tryCatch(fg$getRgb(), error = function(e) "")
    rgb <- paste(rgb, collapse = "")
    return(rgb)
  })
  rowIndex <- as.numeric(unlist(lapply(strsplit(names(bg),split = "\\."),
                                       "[", 1))) # row indexes are at first position
  rowIndex <- rowIndex - min(rowIndex) + 1 # start at 1
  colIndex <- as.numeric(unlist(lapply(strsplit(names(bg), split = "\\."),
                                       "[", 2))) # column indexes are at second position
  colIndex <- colIndex - min(colIndex) + 1 # start at 1
  matBG <- matrix("", ncol = max(colIndex), nrow = max(rowIndex))
  for (i in 1:length(bg)) { # fill background color matrix
    matBG[rowIndex[i], colIndex[i]] <- bg[i]
  }
  flagData <- as.data.frame(matBG)[rowsData, colsData]
  colnames(flagData) <- metabolites

  flagData[flagData == "6a5acd"] <- "LOD"
  flagData[flagData == "87ceeb"] <- "LOQ"
  flagData[flagData == "00cd66"] <- "Valid"
  flagData[flagData == "ffff33"] <- "Out of calibration range"
  flagData[is.na(rawData)] <- NA
  if (verbose) {
    sumLOD <- sum(flagData == "LOD", na.rm = TRUE)
    sumLOQ <- sum(flagData == "LOQ", na.rm = TRUE)
    sumValid <- sum(flagData == "Valid", na.rm = TRUE)
    sumCalRange <- sum(flagData == "Out of calibration range", na.rm = TRUE)
    nas <- sum(is.na(flagData))
    total <- nrow(flagData)*ncol(flagData)
    print(paste0("Number of LODs: ", sumLOD,
                 " (", round(sumLOD/total*100, 2), "%)"))
    print(paste0("Number of LOQs: ", sumLOQ,
                 " (", round(sumLOQ/total*100, 2),"%)"))
    print(paste0("Number of Valids: ", sumValid,
                 " (", round(sumValid/total*100, 2),"%)"))
    print(paste0("Number of calibration range passes: ", sumCalRange,
                 " (", round(sumCalRange/total*100, 2),"%)"))
    print(paste0("Number of NAs: ", nas,
                 " (", round(nas/total*100, 2),"%)"))
    ## checken!!! stimmt irgendwie nicht immer
  }

  outList <- list(rawData, flagData, metabolites, metaData)
  names(outList) <- c("raw data", "flag data", "metabolites", "meta data")
  return(outList)
}


filterIndicators <- function(inList, verbose=TRUE) {
  rawData <- inList[[1]]
  flagData <- inList[[2]]
  metabolites <- inList[[3]]
  metaData <- inList[[4]]

  ## filter matabolism indicators ##
  if ("Metabolism Indicators" %in% names(metabolites)) {
    notIndicators <- which(names(metabolites) != "Metabolism Indicators")
    metabolites <- metabolites[notIndicators]
    if (verbose) {
      print("------------------------------")
      print("Filtering metabolism indicators...")
      print(paste("Number of remaining metabolites:", length(metabolites)))
      print(paste("Number of remaining classes:",
                  length(unique(names(metabolites)))))
    }

    rawData <- rawData[, metabolites]
    flagData <- flagData[, metabolites]

  } else if (verbose) {
    print("------------------------------")
    print("Data doesn't contain metabolism indicators!")
  }

  outList <- list(rawData, flagData, metabolites, metaData)
  names(outList) <- c("raw data", "flag data", "metabolites", "meta data")
  return(outList)
}


# list("Group Column"="WT vs Glo1",
#      "Mutant"=c("WT", "Glo1 KO"),
#      "Treatment"=c("control", "overfed"))
annotateGroups <- function(inList, groupAnnotation, verbose=TRUE) {
  rawData <- inList[[1]]
  flagData <- inList[[2]]
  metabolites <- inList[[3]]


  metaData <- inList[[4]]

  ## annotate with groups
  if (verbose) {
    print("------------------------------")
    print("Annotating data...")
  }

  colGroup <- groupAnnotation[["Group Column"]]
  groupEntries <- metaData[, colGroup]

  groupNameOpt <- names(groupAnnotation)[2:length(groupAnnotation)]
  for (groupName in rev(groupNameOpt)) {
    options <- factor(groupAnnotation[[groupName]])
    if (grepl(options[1], paste(groupEntries, collapse = ""))) {
      pattern <- options[1]
      other <- options[2]
    } else {
      pattern <- options[2]
      other <- options[1]
    }
    groups <- sapply(groupEntries, function(GroupEntry) {
      if_else(grepl(pattern, GroupEntry, ignore.case = TRUE),
              pattern,
              other)
    })
    rawData <- mutate(rawData,
                      tmpName = groups,
                      .before = 1)
    colnames(rawData)[1] <- groupName
    flagData <- mutate(flagData,
                       tmpName = groups,
                       .before = 1)
    colnames(flagData)[1] <- groupName
    metaData <- mutate(metaData,
                       tmpName = groups,
                       .before = 1)
    colnames(metaData)[ncol(metaData)] <- groupName
  }

  if (verbose) {
    print(paste("Annotated colums:", paste(groupNameOpt, collapse = ", ")))
  }

  if (length(groupNameOpt) < 2) {
    groupNameOpt[2] <- "dummy"
    rawData$dummy <- ""
    flagData$dummy <- ""
    metaData$dummy <- ""
  }

  rawData <- arrange(rawData, !!sym(groupNameOpt[1]), !!sym(groupNameOpt[2]))
  flagData <- arrange(flagData, !!sym(groupNameOpt[1]), !!sym(groupNameOpt[2]))
  metaData <- arrange(metaData, !!sym(groupNameOpt[1]), !!sym(groupNameOpt[2]))

  rawData <- rawData[, which(colnames(rawData) != "dummy")]
  flagData <- flagData[, which(colnames(flagData) != "dummy")]
  metaData <- metaData[, which(colnames(metaData) != "dummy")]

  outList <- list(rawData, flagData, metabolites, metaData)
  names(outList) <- c("raw data", "flag data", "metabolites", "meta data")
  return(outList)
}
