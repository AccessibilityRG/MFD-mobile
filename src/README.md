# Running the tool

You can run the mfd_interpolation.py tool with Python interpreter or with a terminal in a folder where the script file is located:
  
  ```
  $ cd /home/mfd-mobile/src
  $ python mfd_interpolation.py
  ``` 

## Parameters

**Before running the tool** you should adjust the necessary parameters and filepaths in the `main()` -function of the [script](mfd_interpolation.py). 
You need to specify parameters and filepaths as explained below. 

### Filepaths

| Filepath | Description | File format or data type used in the code | 
| ----------------------|----------------------------------------------------------------------------------------|------------------------------| 
| hat_fp      | File path to a file that contains the time-specific human activities (e.g. time-use survey data). | Excel file | 
| dps_fp      | File path to disaggregated physical surface layer. | Shapefile |
| cdr_fp      | File path to mobile phone data.  | Excel file |
| tz_fp       | File path to target zones layer (e.g. a 100 meter grid). | Shapefile |
| out_dir     | File path to a folder where the results will be saved. | Folder |
| out_prefix  | Prefix that will be placed in front of the result filename. | String |

### Column names 

You should adjust the column names with following parameters in a way that they are presented in your datafiles. 

#### Human activity data

| Parameter               | Description                                                        | Data type in the datafile | 
| ------------------------|--------------------------------------------------------------------|--------| 
| activity_function_type  | Column that contains information about human activity types        | String | 
| seasonal_factor_col     | File path to disaggregated physical surface layer.                 | Float  |
| spatial_unit_col        | Information about the spatial unit (building or area in our case)  | String |

#### Disaggregated physical surface layer

| Parameter                  | Description                                                        | Data type in the datafile                                   | 
|----------------------------|--------------------------------------------------------------------|-------------------------------------------------------------| 
| source_zone_col_dps        | Source zone ID (Base station ID in our case).                                                                | String or Int     | 
| target_zone_col            | Target zone ID (Grid Cell ID in our case).                                                                   | String or Integer |
| building_height_col        | Column with building height information.                                                                     | Integer or Float  |
| building_landuse_features  | Column in the disaggregated physical surface layer that has information about building and landuse types.    | String            | 

#### Mobile phone data

| Parameter                  | Description                                                        | Data type in the datafile  | 
|----------------------------|--------------------------------------------------------------------|----------------------------| 
| source_zone_col_cdr        | Source zone ID (Base station ID in our case).                      | String or Int              |

#### Target zone spatial layer

| Parameter                  | Description                                                        | Data type in the datafile     | 
|----------------------------|--------------------------------------------------------------------|-------------------------------| 
| target_zone_col_spatial    | Target zone ID (Grid Cell ID in our case).                         | String or Int                 |


### Time and projection parameters

| Parameter | Description | File format or data type used in the code |
|-----------|-------------|-------------------------------------------|
| start_h     | Start hour (e.g. number 9 for 9am). | Integer number |
| end_h       | End hour (e.g. number 18 for 6pm). | Integer number |
| epsg        | EPSG code for projection that will be used in the result Shapefile. | Integer number according the [EPSG codes](http://spatialreference.org/ref/epsg/) |