#A program that takes error-checked nest survival .csv files and converts them into a useable encounter history .inp file

# install/load necessary packages
my.packs <- c('MASS',"msm",'jagsUI')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)


rm(list=ls(all=TRUE)) #Clear computer's memory

#Read in each year of encounter history information
coei_data = read.csv('C:\\Users\\koonslt\\Documents\\Iles\\University\\PhD_Project\\Data\\COEI\\Database\\processed_files\\coei_processed_incage_CS_hatchdate.csv',header = T, sep = ",") #Processed datafile

#Years of study
study_years = c(  1976, 1977,
                  1984, 1985, 1986,
                  1991, 1992, 1993, 1994, 1995, 1996, 1997,
                  2000, 2002, 2003,
                  2009,
                  2010,
                  2011, 2012, 2013, 2014,2015)


#Years with actual encounter history data
study_years_data = c( 1976, 1977,
                      1984, 1985, 1986,
                      1991, 1992, 1993, 1994, 1995, 1996, 1997,
                      2000, 2002, 2003,
                      2009,
                      2010,
                      2011, 2012, 2013, 2014,2015)

#Years with simulated encounter history data
study_years_sim = c(1972, 1973, 1978, 1979, 1980)
coei_data = coei_data[coei_data$Year %in% study_years_data,] #Only use selected study years

#######
#Now that the full data matrix has been compiled, begin constructing the .inp file
#######

inp_colnames = c("nestIDyear", "FirstFound", "LastPresent", "LastChecked", "Fate", "Year", "Year.Scaled","NestAge")
inp = matrix(NA,nrow = length(unique(coei_data$nestIDyear)),ncol = length(inp_colnames))
colnames(inp) = inp_colnames
inp = as.data.frame(inp)
inp$nestIDyear = unique(coei_data$nestIDyear)

for (year in study_years){
    #For years with actual data (which appear in the compiled dataframe)
    if (year %in% study_years_data){

        year_data = coei_data[coei_data$Year == year,]

        nest_ids = unique(as.character(year_data$nestIDyear)) #Unique nests in the current year

        year_data$Int_Fate[is.na(year_data$Int_Fate)] = 1 #Since the first encounter currently has an int fate of "NA," replace it with a value of 1

        #Loop through the nests and generate the encounter histories
        for (nest in nest_ids){
            nest_data = year_data[year_data$nestIDyear == nest,]
            nest_dates = nest_data$Calendar
            if (length(nest_data$Year) == 1) next else { #Skip nests with fewer than 2 visits
                #Set all the values back to NA (to clear the computer's memory of the previous nest)
                first_found = NA
                last_active = NA
                last_checked = NA
                fail = NA

                #Extract relevant hatch dates
                hatch_date = nest_data$hatch.complete[1]
                NestAge = nest_data$Incubation_Age[1]

                if (is.na(hatch_date)) print(paste(nest,"has an error in hatch date"))
                NestAge = nest_data$Incubation_Age[1]

                first_found = as.numeric(min(nest_data$Calendar))  #First found

                fail_vector = nest_data$Fail #The fate of the nest (a value for each visit - will have 1 on a visit if it failed)

                #If the nest failed, 1 will be included in the fail vector (and 0 in the Int_Fate column will indicate when it was discovered as being failed)
                if (1 %in% fail_vector){
                    last_active = as.numeric(max(nest_data$Calendar[nest_data$Int_Fate == 1])) #Select the largest Julian Date for which the nest was active
                    last_checked = as.numeric(max(nest_data$Calendar)) #Last checked
                    if (last_checked > hatch_date) last_checked = hatch_date #Depredated nests could not have survived longer than the hatch date
                    fail = 1
                    if (last_active >= last_checked){
                        last_checked = last_active + 1
                        print (paste("Dates adjusted for nest",nest))}
                }

                #If the nest was successful, its last_active and last_checked dates should be the same (the projected date of hatch of the nest)
                if ((1 %in% fail_vector) == FALSE){
                    fail = 0
                    last_active = hatch_date
                    last_checked = hatch_date

                    #For successful nests, the hatch date cannot be greater than the last day the nest was physically checked
                    if (last_checked > max(nest_data$Calendar)){
                        last_active = max(nest_data$Calendar)
                        last_checked = max(nest_data$Calendar)}
                }

                #The basic information required for nest survival models
                inp$FirstFound[inp$nestIDyear == nest] = first_found
                inp$LastPresent[inp$nestIDyear == nest] = last_active
                inp$LastChecked[inp$nestIDyear == nest] = last_checked
                inp$Fate[inp$nestIDyear == nest] = fail

                #Additional nest-specific information included in these nest survival models
                #Year of study
                inp$Year[inp$nestIDyear == nest] = year
                inp$Year.Scaled[inp$nestIDyear == nest] = year - min(study_years)
                inp$NestAge[inp$nestIDyear == nest] = NestAge

            }
        }
    }
}

inp = inp[,-8]
inp = na.omit(inp)

write.csv(inp,file = "../../../Data/coei_inp.csv", row.names = F)
