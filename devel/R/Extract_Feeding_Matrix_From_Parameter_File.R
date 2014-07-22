# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Extract Feeding Matrix from Parameterfile "NorthSea_biol_fishing_final.prm"

# Author:    Alexander Keth
# Institute: Instutute for Hydrobiology and Fisheries Science (IHF)
# Date:      10.07.2014

# String to test the function!
# filename <- "NorthSea_biol_fishing_final.prm"

# Userinput: Modelpath!
model_path <- "Z://Atlantis//ATLANTIS NSmodel"

get.availability.matrix <- function(filename = "NorthSea_biol_fishing_final.prm")
{
  # Read in Parameterfile
  parfile <- readLines(file.path(model_path, filename), warn = F)
    
  # Extract search tags from the parameterfile! These are used to exract the data lateron and to create the pred/stanua columns in the final result!
  search <- parfile[grep(pattern = "pPREY", x = parfile)]
  search <- search[-c(1,2)]
  search <- str_split(search, pattern = "\t")
  search <- unlist(lapply(search, function(x)x[1]))
      
  # Extract positions in the Parameterfile. NOTE: The values are in the next row!
  pos <- unlist(lapply(search, grep, x = parfile))
  values <- parfile[pos + 1]
  values <- str_split(values, pattern = "\t")
  values <- lapply(values, as.numeric)

  # Transform to dataframe
  # NOTE: 2 Names are appened here namely "fishlarvae" and "import"! This is not needed when all names are provided in "FunctionalGroups.csv"
  result <- data.frame(do.call(rbind, values))
  names(result) <- c(as.vector(read.table(file.path(model_path, "functionalGroups.csv"), sep = ",", header = T)[,1]), "fishlarvae", "import")
  
  # Add Predator/Stanza columns to dataframe!
  result$code <- str_replace_all(search, pattern = "pPREY", replacement = "")  
  result$pred <- with(result, str_replace_all(code, pattern = c("[12]"), replacement = ""))  
  result$prey.cat <- ifelse(str_sub(result$code, start = -1) != "2", "adu", "juv")
  result$pred.cat <- ifelse(str_sub(result$code, start = 1, end = 1) != "2", "adu", "juv")

  # Transform dataframe! Conversion from wide to long! Not needed at the moment!
  # result <- melt(result, id.vars = c("pred", "pred.cat", "prey.cat"), measure.vars = names(result)[1:55], variable.name = "prey", value.name = "avail")
  return(result)
}






