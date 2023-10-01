#setwd("~/bahaha_assembly/synteny/genespace/example/");
library(data.table);
library(GENESPACE)
runwd <- file.path("/data/projects/rcui/bahaha_assembly/synteny/genespace/morespp/rundir")
#list.files(runwd, recursive = T, full.names = F)

# gpar <- init_genespace(
  # genomeIDs = c("BahahaTaipingensis", "LarimichthysCrocea", "CollichthysLucidus", "NibeaAlbiflora", "DatnioidesUndecimradiatus", "LutjanusErythropterus", "CheilinusUndulatus",  "MicropterusSalmoides",  "SanderLucioperca", "GasterosteusAculeatus"),
  # speciesIDs = c("Bahaha_taipingensis", "Larimichthys_crocea", "Collichthys_lucidus", "Nibea_albiflora", "Datnioides_undecimradiatus", "Lutjanus_erythropterus", "Cheilinus_undulatus",  "Micropterus_salmoides",  "Sander_lucioperca", "Gasterosteus_aculeatus"),
  # versionIDs = c("1.0", "JimeiU", "collichthys_lucidus","nibea_albiflora" , "1", "1", "1", "1", "1", "1"),
  # ploidy = rep(1,10),
  # diamondMode = "default",
  # orthofinderMethod = "default",
  # wd = runwd,
  # nCores = 120,
  # minPepLen = 30,
  # gffString = "gff",
  # pepString = "pep",
  # path2orthofinder = "/opt/miniconda3/bin/orthofinder",
  # path2diamond = "/opt/miniconda3/bin/diamond",
  # path2mcscanx = "/data/software/MCScanX/",
  # rawGenomeDir = file.path(runwd, "rawGenomes"))

gpar <- init_genespace(
  genomeIDs = c("BahahaTaipingensis", "LarimichthysCrocea", "NibeaAlbiflora","CollichthysLucidus",  "LutjanusErythropterus", "CheilinusUndulatus",  "MicropterusSalmoides",  "SanderLucioperca", "GasterosteusAculeatus"),
  speciesIDs = c("Bahaha_taipingensis", "Larimichthys_crocea", "Nibea_albiflora","Collichthys_lucidus",  "Lutjanus_erythropterus", "Cheilinus_undulatus",  "Micropterus_salmoides",  "Sander_lucioperca", "Gasterosteus_aculeatus"),
  versionIDs = c("1.0", "JimeiU", "nibea_albiflora" ,"collichthys_lucidus", "1", "1", "1", "1", "1"),
  ploidy = rep(1,9),
  diamondMode = "default",
  orthofinderMethod = "default",
  wd = runwd,
  nCores = 120,
  minPepLen = 30,
  gffString = "gff",
  pepString = "pep",
  path2orthofinder = "/opt/miniconda3/bin/orthofinder",
  path2diamond = "/opt/miniconda3/bin/diamond",
  path2mcscanx = "/data/software/MCScanX/",
  rawGenomeDir = file.path(runwd, "rawGenomes"))
  
#parse_annotations(
#  gsParam = gpar,
#  gffEntryType = "gene",
#  gffIdColumn = "locus",
#  gffStripText = "locus=",
#  headerEntryIndex = 1,
#  headerSep = " ",
#  headerStripText = "locus=")

#gpar <- run_orthofinder(
#  gsParam = gpar)

gpar <- synteny(gsParam = gpar)


arrCol <- rainbow(24);
arrLabelGenomes <- c('BahahaTaipingensis')

regs <- data.table(
  genome = c('BahahaTaipingensis'),
  chr = c('10'), start=c(0), end=(28742000))


datInvert <- data.frame(genome="LarimichthysCrocea", chr=toupper(c('cm014885.1', 'cm014896.1', 'cm014882.1', 'cm014893.1', 'cm014886.1', 'cm014881.1','cm014902.1','cm014900.1')) )
datInvert <- rbind(datInvert, data.frame(genome="CollichthysLucidus", chr=toupper(c('cm014081.1', 'cm014082.1', 'cm014083.1', 'cm014079.1', 'cm014095.1', 'cm014089.1', 'cm014094.1', 'cm014100.1', 'cm014101.1' ) ) ) );
datInvert <- rbind(datInvert, data.frame(genome="NibeaAlbiflora", chr=toupper(c('cm024795.1', 'cm024790.1', 'cm024800.1', 'cm024791.1', 'cm024794.1', 'cm024796.1', 'cm024807.1' ,'cm024806.1', 'cm024805.1', 'cm024802.1') ) ));
datInvert <- rbind(datInvert, data.frame(genome="LutjanusErythropterus", chr=toupper(c('cm034600.1','cm034606.1','cm034604.1','cm034607.1','cm034613.1') )));
datInvert <- rbind(datInvert, data.frame(genome="CheilinusUndulatus", chr=toupper(c('nc_054865.1', 'nc_054868.1', 'nc_054874.1', 'nc_054871.1', 'nc_054882.1', 'nc_054875.1', 'nc_054876.1', 'nc_054887.1', 'nc_054886.1', 'nc_054888.1', 'nc_054878.1'))));
datInvert <- rbind(datInvert, data.frame(genome="MicropterusSalmoides", chr=toupper(c('nw_024043372.1','nw_024044237.1','nw_024040374.1','nw_024040596.1','nw_024040152.1','nw_024044459.1','nw_024041039.1','nw_024041373.1','nw_024041484.1','nw_024040928.1'))));
datInvert <- rbind(datInvert, data.frame(genome="SanderLucioperca", chr=toupper(c('NC_050175.1','NC_050176.1','NC_050178.1','NC_050177.1','NC_050183.1','NC_050180.1','NC_050182.1','NC_050174.1','NC_050185.1','NC_050187.1','NC_050188.1','NC_050189.1','NC_050195.1'))));
datInvert <- rbind(datInvert, data.frame(genome="GasterosteusAculeatus", chr=toupper(c('nc_053214.1','nc_053218.1'))));

#datInvert <- NULL;

plot_riparian(
  gpar, 
  useOrder = FALSE, 
 onlyTheseRegions = regs, 
#  excludeChrOutOfRegion = TRUE, 
#  colByChrs = "white",
  verbose = FALSE)

#pg <- pangenome(gpar)
