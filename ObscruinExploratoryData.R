exceldata <- read_excel("/Users/mill29ca 1/Desktop/Ig1213_analysis_simple.xlsx")
# had to set up excel file to have long list of the data 

domain1_linker_names <- c("Angle 60-90-91", "Angle 60-90-92", "Angle 60-90-93", "Angle 60-90-94", "Angle 60-90-95", "Angle 60-90-96", "Angle 60-90-97")
linker_domain2_names <- c("Angle 91-98-169", "Angle 92-98-169", "Angle 93-98-169", "Angle 94-98-169", "Angle 95-98-169", "Angle 96-98-169", "Angle 97-98-169")

df <- exceldata %>% mutate(Domain = ifelse(Location %in% domain1_linker_names, "Domain1", "Domain2"))

