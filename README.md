# r_customfunctions
R functions to manipulate and create plots 
1- exportrmd - This is a custom function to export the code chunks from the R script, it detects the patterns of code chunks starts with ## ---- chunk and ends with ## ---. It will also load the RData to the rmd, from the parameters. Do add the necessary library to construct the plot.
example usage - exportrmd("/path/to/your/script.R", "/path/to/your/data.RData", "/path/to/output/output.Rmd")
