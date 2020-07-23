This readme describes the structure of the data files. For more details please refer to the paper. 

i_f.csv contains the flood intensity index and the associated seasonality of the historical flood series. 
Description of columns in the file: 
Code: Code of series (see Extended Data Table 1)
Year: Calendar year 
Index: Three-scaled intensity index of flood events based on documentary evidence, comprising 'class1' (notable flood), 'class2' (great flood), 'class3' (extraordinary flood), 'no_info' (no information), 'prob_no_flood' (most probably no flood) and 'no_flood' (no flood). 
Season: Season of flood event

bias.csv contains a bias index reflecting the completeness of the source material in a historical context for a given flood series and year. 
Description of columns in the file: 
Year: Calendar year
Code of series: (103 columns) Bias index reflecting the completeness of the source material in a historical context for a given series and calendar year, comprising 1 (no data), 2 (periods with possibly missing data), 3 (average) and 4 (periods with overly dense data compared to the average of the series)

sites.csv contains meta information about the historical flood series including geographical coordinates. 
Description of columns in the file: 
code: Code of series (see Extended Data Table 1)
country: Country associated with flood series (see Extended Data Table 1)
river: River flood series (see Extended Data Table 1)
location: Location associated with flood series
lat1: Latitude (Degrees °)
lat2: Latitude (Minutes ')
lat3: Latitude (Seconds '') 
lon1: Longitude (Degrees °)
lon2: Longitude (Minutes ')
lon3: Longitude (Seconds '')
lon4: Longitude (East or West of Greenwich)
Repres: Representativeness index, reflecting the degree of data representativeness in a regional context, comprising 1 (low representativeness), 2 (average representativeness) and 3 (high representativeness)
Interp: Binary indicator: 1  means the series is used for the interpolation of the flood intensities, 0 means the series is supplementary and only used for the seasonality analysis 
regionEU: Regions associated with flood series, comprising n (Northern Europe), e (Eastern Europe), c (Central Europe), w (Western Europe) and s (Southern Europe). See Extended Data Figure 1

Delta_T.csv contains a 500-year Central European temperature reconstruction series. 
Description of columns in the file: 
Year: Calendar year
temp_anom: Air temperature deviations from the mean (1961-1990) in degree Celcius

500yr_literature table.pdf gives a description of the sources of the historical flood series. See Extended Data Table 1
Description of columns in the file:
Country: Country associated with flood series (see Extended Data Table 1)
River flood series: River associated with flood series (see Extended Data Table 1)
Series No.: Code of series (see Extended Data Table 1)
Details about the series sources and methodological data: Data sources including literature references and published data
