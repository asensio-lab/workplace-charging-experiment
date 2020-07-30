# A Field Experiment on Workplace Norms and Electric Vehicle Charging Etiquette

This repository contains the replication code used in the "A Field Experiment on Workplace Norms and Electric Vehicle Charging Etiquette" paper as well as the anonymized version of the data upon which the paper is built. It is designed for easy replication of all figures and tables.

The repository contains the R code necessary to replicate all figures and tables from the paper. The anonymized data for replication is available at the following address:

https://doi.org/10.7910/DVN/NFPQLW

# Code Instructions & Information

All figures and tables are generated from the single R file saved in this folder, "**package.R**". Running the full file sequentially will automatically generate and save all figures and tables that are derived from data (including all but Figure 1, which is text-based).

## Dependencies

The following R packages are used in this analysis:
1. "rdd"
2. "rddensity"
3. "ggplot2"
4. "car"
5. "dplyr"
6. "lubridate"
7. "dummies"
8. "lpdensity"
9. "extrafont" (optional for replication)
10. "gridExtra"

## Navigation

Once the code has executed, each table and figure from the paper can be found it its own variable.

The naming conventions for these variables are designed to be intuitive and searchable, with the following naming conventions:

For figures -- "figure_#"

For tables -- "table_#"

In cases where a figure has multiple panes, the format is as follow: "figure_{figure #}{pane letter}". For example, Figure 2A appears as "figure_2a".

Figures may be generated non-sequentially in some cases and are organized by the raw data used to generate each rather than by the order in which they appear in the paper. 

# Data Description & Dictionary

This dataset contains information from 3,395 high resolution electric vehicle charging sessions. The data contains sessions from 85 EV drivers with repeat usage at 105 stations across 25 sites at a workplace charging program. The workplace locations include facilities such as research and innovation centers, manufacturing, testing facilities and office headquarters for a firm participating in the U.S. Department of Energy (DOE) workplace charging challenge. The data is in a human and machine readable .CSV format. The resolution of the data is to the nearest second, which is the same resolution as used in the analysis of the paper. It is directly importable into free software. All data is fully anonymized. 

## Main Data ("station_data.csv")

This file contains all data that provides the basis for our analysis, all from the firm within in our study. It is used in generating all figures and tables with the exception of Figure 1, which is not data dependent, and Figure S4b, which comes from an external sample of chargers. The data for generating Figure S4b may be made available upon request by contacting Omar Asensio, whose information can be found at this address: https://datasciencepolicy.gatech.edu/. 

### Data Dictionary

1. *sessionId*: a random number used to identify a specific electric vehcile (EV) charging session
2. *kwhTotal*: the total energy use of a given EV charging session, measured in kWh
3. *dollars*: the amount paid by the user for a given charging session, measured in U.S. dollars
4. *created*: the date and time a given session was initiated, expressed in the form "YYYY-MM-DD HOUR:MIN:SEC"
5. *ended*: the date and time a given session was terminated
6. *startTime*: the hour of day, from 1 to 24, during which the session was initiated
7. *startTime*: the hour of day, from 1 to 24, during which the session was terminated
8. *chargeTime*: the total duration of the session in hours
9. *weekday*: the day of the week on which the charging session took place
10. *platform*: the digital platform used by the EV driver to log the session
11. *distance*: the distance from a user's home to the charging location, expressed in miles except where user did not report address
12. *userId*: a random number used to uniquely identify a given user and his or her transactions
13. *stationId*: a random number unique to each specific charging station
14. *locationId*: a random number unique to a specific location owned by the firm, where chargers were installed
15. *managerVehicle*: a binary variable; 1 if the vehicle logging the transcation is of the type largely used by firm managers, 0 otherwise
16. *facilityType*: a categorical variable indicating the type of facility a station is installed at; manufacturing = 1, office = 2, research and development = 3, other = 4
17. *Mon, Tues, Wed, Thurs, Fri, Sat, Sun*: binary variables for day of week of a given session; 1 if the session occurred on that day of week, 0 otherwise
18. *reportedZip*: binary variable for if a user reported his or her zip code; 1 if zip code was reported, 0 otherwise
