### Instructions and sample code for installing the MEM package and running MEM analysis

##################################################################################################################
##################################################################################################################
###
### Part 0: Notes for Newbies / Refresher
###
##################################################################################################################
##################################################################################################################

##      Experienced users can skip Part 0, but don't forget to get the latest R/RStudio & set your path (Part 1).
##      To skip ahead, just press CTRL-ENTER here and you will go to Part 1.

##      You will use two programs for this script:
##      1) R, a stats programming environment that is primarily command line
##         (command line means typing instructions vs. a graphical user interface or GUI)
##      2) RStudio (likely what you are reading this in), a GUI for R.
##      You may also be interested in Shiny, a web app framework for R.  https://shiny.rstudio.com/

##      Most of the time people struggle with MEM, they need to re-install R and RStudio to get latest versions.
##      If you have an old version of R and/or RStudio, please download and install the latest (vs. "updating").
##      This will help ensure you're using the same version.  We've seen the updates be a little buggy.

##      To run a line in this script, hold "CONTROL" and press "ENTER" on a line to run that code.
##      Pressing CONTROL-ENTER on a line in this script window is a shortcut to type that line in the Console below.
##      When code is running for a while, you will see install bars or sometimes a stop sign in the Console below.

##      Lines here that start with the # symbol are "comments" -- they do not get executed.
##      You can take notes in the script if you "comment out" that line by putting a # at the start
##      You only need one # to comment out a line.  We're just being stylish with two...

##      To run RStudio as an administrator, right click the icon and choose "Run as Administrator"
##      You can also right click the icon, choose "Properties", then in the "Shortcut" tab, then
##      click the "Advanced" button and check "Run as Administrator" to always run as admin.

##################################################################################################################
### End of Part 0: Notes for Newbies / Refresher
##################################################################################################################

#
#
#
#
#

##################################################################################################################
##################################################################################################################
###
### Part 1: Set your MEM path and test1 path (do this every time)
###
##################################################################################################################
##################################################################################################################

#####################################################
##      ** Part 1 Notes **
##      The first thing you need to do is find the "path" to the MEM code on your computer.
##      The path is how to tell R/RStudio/MEM where to find your files.
##      Setting your path need to be done every time you use this script (other sections may be skipped).

##      If you put MEM into a folder on a PC Desktop called "2017 MEM Course" and your User name is "Jonathan"
##      then your path might be "C:/Users/Jonathan/Desktop/2017 MEM Course/MEM"
##      You need to edit the line below and run this line every time you do this script.

#####################################################
## Set MEM_path here and run this line each time
MEM_path = "C:/Users/Jonathan/Dropbox/2017 MEM Course/MEM"

#####################################################
## Make sure your MEM_path is working.  You should see a list of files that includes "MEM.Rproj" and "vignettes".
list.files(MEM_path)

##      If you are having trouble with MEM_path, uncomment the two lines of code below and try them.
##      This will only work if you launched this script from the folder "above" the MEM folder.
##      If this is the case, you should see a folder called MEM in the Files tab in the lower right RStudio window.
##      Running list.files(MEM_path) should show a list of files that includes "MEM.Rproj" and "vignettes".

# MEM_path = "./MEM"
# list.files(MEM_path)

#####################################################
## Set MEM_path_test1 here and run this line each time
MEM_path_test1 = paste(MEM_path, "/test1", sep="")

#####################################################
## Test that your MEM_path_test1 is working.  
list.files(MEM_path_test1)

##      If you get "Character(0)", this either means the test1 folder is empty or the path is wrong.
##      If you just made the folder, you should be all set (it actually is empty).  You use this in Part 4.
##      You may also see test files for the course like "CD4Tcells.fcs", or your own data.

#####################################################
## If you don't have a test1 folder in your MEM folder, make one and repeat the last two lines
shell(gsub("/", "\\\\", paste("explorer ", MEM_path, sep="")), intern=TRUE)

##      This line of code should open the MEM_path directory for Windows users.  Other users must open it manually.
##      If needed, make a new, empty folder here called "test1".  "Here" is your MEM_path (the "MEM" folder).

##################################################################################################################
### End of Part 1: Set your MEM path and test1 path (do this every time)
### Once your paths are set, if you have done Part 2 before, skip to Part 3, otherwise, continue on to Part 2.
##################################################################################################################

#
#
#
#
#

##################################################################################################################
##################################################################################################################
###
### Part 2: Install Packages (skip to Part 3 if you have already done this Part 2, recently)
###
##################################################################################################################
##################################################################################################################

#####################################################
##      ** Part 2 Notes **
##      Part 2 installs key packages (gplots, flowCore, and MEM) and only needs to be run once.
##      Your R should also needs four base packages: stats, tools, grDevices, and utils.  
##      These packages should come with R and if you don't have them, get R again (Part 0)
##      These take a little while and you should do them one at a time and wait for each to finish.

#####################################################
## Optional: Run RStudio as an administrator (see Part 0) & use the latest versions of MEM, R, and RStudio.

## This will open the page to download the latest MEM, once you get it here, you will be emailed when we update it
browseURL("https://mem.vueinnovations.com/", browser = getOption("browser"), encodeIfNeeded = FALSE)

## This will open the page to download the latest R
browseURL("http://cran.r-project.org/bin/windows/base/", browser = getOption("browser"), encodeIfNeeded = FALSE)

## This will open the page to download the latest RStudio
browseURL("http://www.rstudio.com/products/RStudio/", browser = getOption("browser"), encodeIfNeeded = FALSE)

#####################################################
## Install Step 1 - Install gplots
install.packages("gplots")

#####################################################
## Install Step 2 - Install flowCore, entering "a" to install all the packages, if asked
source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")

#####################################################
## Install Step 3 - Install the MEM package (now we're gettin' somewhere...)
install.packages(MEM_path, type="source", repos=NULL)

##################################################################################################################
### End of Part 2: Install Packages (skip to here if you have done Part 2 already, recently)
##################################################################################################################

#
#
#
#
#

##################################################################################################################
##################################################################################################################
###
### Part 3: Run an internal MEM dataset using peripheral blood mononuclear cells (PBMC) as a quick test
###
##################################################################################################################
##################################################################################################################

#####################################################
##      ** Part 3 Notes **
##      Part 3 loads an internal PBMC dataset.

##      For computational people:
##      PMBC = periperhal blood mononuclear cells (in this case from a healthy human blood donor)
##      This dataset is one big table where columns represent measured features (proteins) or information known.
##      Rows represent cell events.  The value in the cell is a mass cytometry "intensity" measurement.
##      The last column is a cluster indetifier for each cell's population ID.
##      Arcsinh data scale: https://my.vanderbilt.edu/irishlab/protocols/scales-and-transformation/

##      For flow cytometry people:
##      The population gating is indicated by a number for each cell in the last column of the data table. 
##      The other way to run MEM is in Part 4, where you have each population of cells in a different file.

##      This part is a quick check that MEM is working and can be skipped once it runs once.
##      This dataset is "internal", so you don't actually need the path working for this to load.
##      If Part 3 works and Part 4 does not, it usually means it is a problem with a path or your files.

#####################################################
## Test Step 1 - Load libraries
library(flowCore)
library(MEM)
library(gplots)

##     It's OK if you get a note about the lowess package being masked (this is normal and won't be a problem).
##     Errors usually can be fixed by updating the path (Part 1) or installing something (Part 2).

#####################################################
## Test Step 2 - Load MEM help to test that MEM is installed
?MEM

##      You should see the window in the lower right say "Marker Enrichment Modeling"
##      If this works: You have now gotten 95% of the way there!  Congratulations!  The rest is loading data and syntax.

#####################################################
## Test Step 3 - Run MEM on the internal PBMC dataset
MEM_values = MEM(PBMC, transform = TRUE, cofactor = 15, choose.markers = FALSE, choose.ref = FALSE, rename.markers = FALSE, file.is.clust = FALSE, add.fileID = FALSE, IQR_thresh = NULL)

#####################################################
## Test Step 4 - Load the help page for the build.heatmaps function
?build.heatmaps

##      You should see the window in the lower right say "Build heatmaps"

#####################################################
## Test Step 5 - Get your first MEM heatmaps with the internal dataset
build.heatmaps(MEM_values, cluster.MEM = "both", cluster.medians = "none", display.thresh = 0, newWindow.heatmaps=TRUE, output.files = FALSE )

##      This pops open two new windows: a yellow median heatmap and a clustered, blue and yellow MEM heatmap.
##      You should get 7 populations in rows with 25 protein features on the bottom.
##      If this works, you've just completed running MEM! Woot!
##      The population with CD8+10 is CD8 T cells.  The population with CD16+5 and no CD3 enrichment is NK cells.
##      Yes, the heatmap color scales are terrible.  We're working on it.
##      Yes, the yellow and black heatmap takes a while to load the color scale... sorry about that.

##################################################################################################################
### End of Part 3: Run an internal MEM dataset
##################################################################################################################

#
#
#
#
#

##################################################################################################################
##################################################################################################################
###
### Part 4: Analyzing a set of new files in MEM using FCS files 
###
##################################################################################################################
##################################################################################################################

#####################################################
##      ** Part 4 Notes **
##      Part 4 takes FCS files from the MEM paper and gets the MEM labels in Figure 1.
##      It will also make a PDF with the clustered MEM heatmap and MEM labels.

#####################################################
## Make sure your MEM_path and MEM_path_test1 variables are set
print(MEM_path)
print(MEM_path_test1)

#####################################################
## Make sure your libraries are loaded
library(flowCore)
library(MEM)
library(gplots)

##      You should see this print out in the Console the correct path, e.g.:
##      [1] "C:/Users/Jonathan/Dropbox/2017 MEM Course/MEM"
##      [1] "C:/Users/Jonathan/Dropbox/2017 MEM Course/MEM/test1"
##      If this fails or looks wrong, go to the first step and set your path

#####################################################
## Data Step 1 - Get flow cytometry files and put them in the test1 directory
##      We willl get the PBMC FCS files from Figure 1 of Diggins et al., Nature Methods 2017
##      and put these files into the "test1" folder that you made (see Part 0 as needed).

## This will open the page to download the Figure 1 test data
browseURL("https://flowrepository.org/experiments/1219/download_ziped_files", browser = getOption("browser"), encodeIfNeeded = FALSE)

##      FlowRepository ID FR-FCM-ZY63
##      Click "I'm not a Robot" and the button to "Zip and Download Files"
##      There should be 7 files named things like "CD4Tcells.fcs" totaling 8.1 in a ZIP
##      The downloaded file is called "FlowRepository_FR-FCM-ZYYJ_files.zip"
##      Find this file after you download it and get the files out of it and into your test1 directory.

#####################################################
## Data Step 2 - Put the unzipped files in the test1 directory
shell(gsub("/", "\\\\", paste("explorer ", MEM_path_test1, sep="")), intern=TRUE)

##      This line of code should open the directory for lazy Windows users.  Other users must open it themselves.

##      Unpack and put the flow cytometry data files into "./2017 MEM Course/MEM/test1"
##      Each file here represents a population of cells (and so we will set MEM's file.is.clust to TRUE).
##      There must be at least two FCS files in test1 for the following code to work.

##      If there are already files in test1, you should move or delete them.
##      Make sure there are not extra files with mis-matched channels in this directory or the code will fail.

##      The full FlowRepository page with the Figure 1 test data is here:
##      browseURL("https://flowrepository.org/id/FR-FCM-ZY63", browser = getOption("browser"), encodeIfNeeded = FALSE)

list.files(MEM_path_test1)

#####################################################
## Data Step 3 - Set MEM_path_test1 as your working director and read all the .FCS files into a data_set variable
setwd(MEM_path_test1)
data_set <- dir(pattern=".fcs")

##      MEM accepts FCS, CSV, and TXT file formats. 
##      Data format: cells in rows, features/markers in columns.
##      Normally you will have one file per population, in which case you are treating files as clusters
##      Alternatively, you can set up a new column in your file with the last column as an integer cluster id / index 
##      Change file extension pattern = "x" depending on file type (i.e. ".fcs", ".csv", ".txt")

##      You should see the Environment window in RStudio update to say data_set and list the files
##      This might be "CD4Tcells.fcs" and other file names.
##      If data_set in the Environment tab is "character (empty)", it didn't see any files in test1.
##      Note: MEM also accepts matrix or dataframe objects (if you have already loaded your data into R)

#####################################################
## Data Step 4 - Run MEM on the data_set
MEM_values = MEM(data_set, transform = TRUE, cofactor = 15, choose.markers = TRUE, choose.ref = FALSE, rename.markers = FALSE, file.is.clust = TRUE, add.fileID = FALSE, IQR_thresh = NULL)

##      For detailed explanation of the arguments, enter ?MEM in console.

##      If this works, you'll see a list of "column names" (the things measured for each event, e.g. proteins on cells)
##      You next need to enter which colun numbers to use -- go to Data Step 5.

#####################################################
## Data Step 5 - When it asks about channels in the Console, copy/paste or type the the following numbers:
##       12,14:17,19,20,22:23,26:29,31:33,35,36,38:40

##      This picks the channels / columns / proteins that have information and are on comparable scales.
##      You can use all the channels, but the scale for Time will "set" the +10/-10 scale for MEM.
##      Keep in mind that data with very different measurement scales will need some transformation before MEM will work well.

#####################################################
## Data Step 6 - Build the heatmaps, this time also making an output folder with PDFs in it
build.heatmaps(MEM_values, cluster.MEM = "both", cluster.medians = "none", display.thresh = 0, newWindow.heatmaps=TRUE, output.files = TRUE)

#####################################################
## Data Step 7 - Open the PDFs.  For lazy Windows users, you can use the following.  Sorry Mac+ users.
shell(gsub("/", "\\\\", paste("explorer ", paste(MEM_path_test1, "/output files", sep=""), sep="")), intern=TRUE)

#####################################################
## Data Step 8 - Compare to Figure 1 of Diggins et al., Nature Methods 2017, you should have similar MEM labels
browseURL("http://www.nature.com/nmeth/journal/v14/n3/fig_tab/nmeth.4149_F1.html", browser = getOption("browser"), encodeIfNeeded = FALSE)

##       This line opens the "output files" subdirectory within the test1 folder within your MEM_path
##       This will only work on Windows.  The gsub part is to replace the / with a \ for Windows
##       MEM values are be written to a text file in "output files".  The heatmaps are saved as PDFs.

##       Look for the "2017-04-15_094226 MEM heatmap.pdf" -- this has MEM labels and the heatmap
##       We're still working on getting beautifully formatted MEM labels and figures.  Let us know if you want to help!
##       You will also see a text file you can open in Excel with the MEM exponents for each channel vs. each population.
##       We're still working on automatically identifying each population... they will be numbered and you need to ID.

##       You can also replace the FCS files in test1 with your own files.  Have fun!

#####################################################
## Data Step 8 - Write Cell, Science, or Nature paper, referencing Diggins et al., Nature Methods 2017

##       Let us know if you encounter any issues!  
##       Official MEM code: https://mem.vueinnovations.com/

##################################################################################################################
### End of Part 4. Thank you for using MEM and for doing great science. 
##################################################################################################################