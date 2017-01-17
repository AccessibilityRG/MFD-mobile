# Running the tool

You can run the mfd_interpolation.py tool with Python interpreter or with a terminal in a folder where the script file is located:
  
  ```
  $ cd /home/mfd-mobile/src
  $ python mfd_interpolation.py
  ``` 

## Parameters

**Before running the tool** you should adjust the necessary parameters and filepaths in the `main()` -function of the [script](mfd_interpolation.py). 
You need to specify following parameters and filepaths:

| Parameter or filepath | Description | File format or data type used in the code | 
| ----------------------|----------------------------------------------------------------------------------------|------------------------------| 
| time_use_fp | File path to a file that contains the time-specific human activities (e.g. time-use survey data). | Excel file | 
| dps_fp      | File path to disaggregated physical surface layer. | Shapefile |
| cdr_fp      | File path to mobile phone data.  | Excel file |
| tz_fp       | File path to target zones layer (e.g. a 100 meter grid). | Shapefile |
| out_dir     | File path to a folder where the results will be saved. | Folder |
| out_prefix  | Prefix that will be placed in front of the result filename. | String |
| start_h     | Start hour (e.g. number 9 for 9am). | Integer number |
| end_h       | End hour (e.g. number 18 for 6pm). | Integer number |
| epsg        | EPSG code for projection that will be used in the result Shapefile. | Integer number according the [EPSG codes](http://spatialreference.org/ref/epsg/) |