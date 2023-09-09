#' Export an Rmd file from an R script and RData file
#'
#' This function generates an Rmd file with code chunks and setup.
#'
#' @param script_file Path to the R script file.
#' @param rdata_file  Path to the RData file to be loaded.
#' @param output_file Path to the output Rmd file.
#'
#' @export
exportrmd <- function(script_file, rdata_file, output_file) {
  # Read the content of the R script file
  script_content <- readLines(script_file, warn = FALSE)
  
  # Create a new Rmd file
  rmd_output <- file(output_file, "w")
  
  # Add the initial code chunk to set up knitr options and load RData
  cat("```{r , echo = FALSE, message = FALSE}\n", file = rmd_output)
  cat("knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.width = 15, fig.height = 10)\n", file = rmd_output)
  cat("load(\"", rdata_file, "\")\n", sep = "", file = rmd_output)
  cat("```\n\n", file = rmd_output)
  
  # Initialize a flag to determine if a code chunk is in progress
  in_chunk <- FALSE
  
  # Loop through the lines in the script file
  for (line in script_content) {
    # Check if the line starts a code chunk
    if (grepl("^## ---- chunk", line)) {
      # If a code chunk is in progress, close it
      if (in_chunk) {
        cat("```\n\n", file = rmd_output)
      }
      
      # Extract the chunk label
      chunk_label <- sub("^## ---- chunk", "", line)
      cat("```{r ", chunk_label, "}\n", file = rmd_output)
      in_chunk <- TRUE
    } else if (in_chunk) {
      # Write the code line into the current code chunk
      cat(line, "\n", file = rmd_output)
      
      # Check if the line ends the code chunk
      if (grepl("^## ---", line)) {
        cat("```\n\n", file = rmd_output)
        in_chunk <- FALSE
      }
    }
  }
  
  # Close the Rmd file
  close(rmd_output)
  
  cat("Rmd file saved as:", output_file, "\n")
  
}


