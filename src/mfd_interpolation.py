"""
mfd_interpolation.py


Copyright (C) 2017.  Digital Geography Group / Accessibility Research Group, University of Helsinki (Järv, Tenkanen & Toivonen).
Programmed by: Henrikki Tenkanen, University of Helsinki, Finland.


This script is part of the following article:

  Järv Olle, Tenkanen Henrikki & Toivonen Tuuli (2017). Enhancing spatial accuracy of mobile
  phone data using multi-temporal dasymetric interpolation. Published in International Journal of 
  Geographical Information Science.


  PURPOSE:
  --------

  This script demonstrates how time-dependent dasymetric interpolation model can be used to interpolate population distributed in spatially varying
  coverage areas of mobile phone base stations into a desired layout of predefined statistical units using ancillary data sources. 

  REQUIREMENTS:
  -------------

  Python 3 with following packages and their dependencies: pandas, geopandas.
  
  DATA:
  -----

  This script is planned to work with data that is described in the main article, see Järv et al. (2017). In short, the script requires three datasets:
      
      1) a disaggregated physical surface layer, 
      2) time-dependent human activity data,
      3) mobile phone data (CDR), or similar. 

  CONTACT:
  --------
  
  You can contact Olle Järv or Henrikki Tenkanen if you have any questions related to the model or the Python implementation of the model:
      
     - olle.jarv (a) helsinki.fi
     - henrikki.tenkanen (a) helsinki.fi
     - http://www.helsinki.fi/science/accessibility
  

  LICENCE:
  --------

  mfd_interpolation.py by Accessibility Research Group (University of Helsinki) is licensed under
  a Creative Commons Attribution 4.0 International License.
  More information about license: http://creativecommons.org/licenses/by/4.0/

"""

import pandas as pd
import geopandas as gpd
from fiona.crs import from_epsg
import os

def main():    
    
    """ Main method that controls the Multi-temporal function-based dasymetric interpolation model (MFD interpolation). """
    
    # ------------------------------------------------------------------
    # 1. Read input data
    # -------------------------------------------------------------------
    
    # File paths and parameters
    # ............................
    
    # Time use survey
    time_use_fp = r"P:\h510\metropaccess\Aineistot\DYNAMO\Timeuse\TimeUse.xlsx"
    
    # Disaggregated physical surface layer ( landuse + buildings + coverage areas + predefined statistical units )
    dps_fp = r"P:\h510\metropaccess\Aineistot\DYNAMO\SpatialDomain\LandUseBuildingsGridMobileVoronoi_500m_Union.shp"
    
    # Mobile phone data (Call Detail Records)
    cdr_fp = r"P:\h510\metropaccess\Aineistot\DYNAMO\MPdata\MP_24h.xlsx"
    
    # Target zones (Predefined spatial units)
    tz_fp = r"P:\h510\metropaccess\Aineistot\DYNAMO\OD\Grid500m_popRegister_with_both_GRD_IDs.shp"
    
    # Output folder for the results
    out_dir = r"C:\HY-Data\HENTENKA\Testaus\Dasymetic"
    
    # Prefix for the output name (time info for the filename will be added automatically)
    out_prefix = "ZROP_test"
    
    # Start hour
    start_h = 9
    
    # End hour
    end_h = 16
    
    # EPSG code for desired output projection
    epsg = 3301
        
    # ----------------------------------------------------------
    
    print("Running Multi-temporal dasymetric interpolation ...")
    
    # Iterate over the desired hours of the day
    for xhour in range(start_h, end_h+1):
        print("Processing hour: %s" % xhour)
        
        # tu = time use
        # dps = disaggregated physical layer
        # cdr = mobile phone data
        # output = output spatial layer in statistical units
        tu, dps, cdr, output = readFiles(time_use_fp=time_use_fp, dps_fp=dps_fp, cdr_fp=cdr_fp, tz_fp=tz_fp)
    
        # Rename Site_ID to site_id
        dps = dps.rename(columns={'Site_ID': 'site_id'})
    
        # Use time window 4pm - 5 pm for the whole analysis (an example)
        # ..............................................................
        time_window = 'H%s' % xhour
                
        # --------------------------------------------------------------------
        # 2. Calculate Relative share of Mobile Phone users (RMP)
        # --------------------------------------------------------------------
            
        # Calculate RMP - i.e. normalize the Mobile Phone user counts to scale 0.0 - 1.0
        # Note: In the here this part is done earlier than in the manuscript (--> chapter 3.4) for practical reasons. 
        # ..........................................................................................
        
        # Name of the mobile phone user counts per site_id in the table
        twm = time_window + 'm'
        cdr = calculateRMP(cdr, time_window=time_window)
    
        # --------------------------------------------------------------------
        # 3. Reclassify Landuse layer ( based on Open Street Map information )
        # --------------------------------------------------------------------
    
        dps = reclassifyLanduse(dps_df=dps)
    
        # --------------------------------------------------------------------
        # 4. Join layers into same DataFrame
        # --------------------------------------------------------------------
    
        # Join necessary columns from <tu> table
        # .........................................
    
        dps = dps.merge(tu[[time_window+'t', 'ActivityLocation', 'Type', 'Type_weight']],
                            left_on=['AFT', 'SPUT', 'SF'],
                            right_on=['ActivityLocation', 'Type', 'Type_weight'])
    
    
        # Join necessary columns from <cdr> table
        # .........................................
    
        dps = dps.merge(cdr[[twm, 'RMP %s' % twm, 'site_id']], on='site_id')
    
        # ---------------------------------------------------------------------
        # 5. Calculate the Relative Floor Area (RFA)
        # ---------------------------------------------------------------------
    
        RFA = calculateRelativeFloorArea(dps)
    
        # ---------------------------------------------------------------------
        # 6. Calculate the Estimated Human Presences (EHP)
        # ---------------------------------------------------------------------
    
        EHP = calculateEHP(df=RFA, time_window=time_window)
        
        # ----------------------------------------------------------------------
        # 7. Calculate Relative Observed Population (ROP)
        # ----------------------------------------------------------------------
    
        ROP = calculateROP(df=EHP, time_window=time_window)
    
        # -----------------------------------------------------------------------
        # 8. Aggregate spatially to desired target zones (ZROP)
        # -----------------------------------------------------------------------
    
        ZROP = calculateZROP(df=ROP, time_window=time_window, grid_cell_id_col='GRD_INSPIR')
    
        # -----------------------------------------------------------------------
        # 9. Save result to disk in Shapefile format
        # -----------------------------------------------------------------------
        out_filename = "%s_%s.shp" % (out_prefix, time_window)
        out = os.path.join(out_dir, out_filename)
    
        # Save file to disk
        saveToShape(input_df=ZROP, grid_df=output, output_path=out, grid_cell_id='gridcode', epsg_code=epsg)
        

def readFiles(time_use_fp=None, dps_fp=None, cdr_fp=None, tz_fp=None):
    """ Read files into memory that are needed for Multi-temporal Dasymetric Interpolation """
    # Read input files
    time_use = pd.read_excel(time_use_fp)
    dps = gpd.read_file(dps_fp)
    cdr = pd.read_excel(cdr_fp, sheetname=0)
    tz = gpd.read_file(tz_fp)
    return time_use, dps, cdr, tz

def reclassifyLanduse(dps_df=None):
    """ Function reclassifies OSM landuse attributes based on criteria in Table S2."""

    # Create new columns for activity location categories
    # ----------------------------------------------------

    # AFT ==> Activity function type == 'ActivityLocation' in <time>
    dps_df['AFT'] = None

    # SPUT ==> Spatial Unit Type == 'Type' in <time>
    dps_df['SPUT'] = None

    # SF ==> Seasonal Factor M == 'Type_weight' in <time>
    dps_df['SF'] = None

    # Execute the classification row by row
    # --------------------------------------
    dps_df = dps_df.apply(classify, axis=1, incol='TYYPg', activityType='AFT', spatialUnit='SPUT', seasonalFactor='SF')  

    return dps_df

def classify(row, incol, activityType, spatialUnit, seasonalFactor):
    """ Reclassifies landuse features to following activity types (AT). See Table S2. """

    # RESIDENTIAL
    if row[incol] in ['Accommodation', 'Residential', 'ResiArea']:
      row[activityType] = 'Residential'

      if row[incol] == 'ResiArea':
        row[spatialUnit] = 'area'
        row[seasonalFactor] = 0.1
      else:
        row[spatialUnit] = 'building'
        row[seasonalFactor] = 0.9

    # WORK
    elif row[incol] in ['Work', 'WorkArea', 'Educational']:
      row[activityType] = 'Work'

      if row[incol] == 'WorkArea':
        row[spatialUnit] = 'area'
        row[seasonalFactor] = 0.1
      else:
        row[spatialUnit] = 'building'
        row[seasonalFactor] = 0.9

    # RETAIL & SERVICE
    elif row[incol] == 'Retail':
      row[activityType] = 'Retail'
      row[spatialUnit] = 'building'
      row[seasonalFactor] = 1.0

    # OTHER
    elif row[incol] in ['Eatery', 'Leisure', 'Public', 'GreenArea', 'Other_']:
      row[activityType] = 'Other'

      if row[incol] == 'GreenArea':
        row[spatialUnit] = 'area'
        row[seasonalFactor] = 0.1
      else:
        row[spatialUnit] = 'building'
        row[seasonalFactor] = 0.9

    # TRANSPORT (i.e. Road in the data)
    elif row[incol] == 'Road':
      row[activityType] = 'Road'
      row[spatialUnit] = 'area'
      row[seasonalFactor] = 1.0

    # RESTRICTED
    elif row[incol] == 'Restricted':
      row[activityType] = 'Restricted'
      row[spatialUnit] = 'all'
      row[seasonalFactor] = 0.0

    return row

def calculateRelativeFloorArea(df=None):
    """ 
    Calculate the relative floor area (RFA) for each subunit within a base station. 
    See chapter 3.2 in the article + chapter S2.2 and Table S2 in the supplementary materials.
    """

    # Floor ==> Calculate the number of floors in a building based on mean height
    df = calculateFloors(df=df, height_col='MEAN', type_col='Type4E', target_col='Floor')
    
    # FEA ==> Calculate the Feature Area (i.e. Area of UNION polygon)
    df['FEA'] = None
    df = calculateArea(df, geom_col='geometry', target_col='FEA')

    # FA ==> Floor Area (i.e. [Area of UNION polygon] * [Number of floors in building] )
    df['FA'] = df['Floor'] * df['FEA']

    # SSFA ==> Sum Site Floor Area (i.e. Sum 'FA' ("Floor Area") by 'Site_ID' of mobile phone cells)
    # --------
    df['SSFA'] = None

    # Group data by 'Site_ID'
    grouped = df.groupby('site_id')

    # Iterate over groups and sum the values
    for key, values in grouped:
      # Sum the 'Area Floor'
      ssaf = values['FA'].sum()

      # Get the indices of the values
      siteid_indices = values.index

      # Assign value to column 'SSFA'
      df.loc[siteid_indices, 'SSFA'] = ssaf

    # RFA ==> Relative Floor Area for each subunit within a base station (scale 0.0 - 1.0)
    df['RFA'] = None
    df['RFA'] = df['FA'] / df['SSFA']

    return df

def calculateFloors(df, height_col, type_col, target_col):
    """ 
    Calculate the number of floors in a building based on <height_col> into <target_col> based on building type <type_col>.
    See chapter 3.1 in the article + chapter S2.1 and Table S2 in the supplementary materials.
    
    """
    # Column for 'Floors'
    df[target_col] = None

    # Calculate the amount of Floors based on MEAN hight of the building
    df = df.apply(returnFloors, axis=1, height_col=height_col, type_col=type_col, target_col=target_col)

    # Set Floor value to 1 for values less than zeros
    df.ix[df[target_col]<1, target_col] = 1

    # Return DataFrame
    return df

def returnFloors(row, height_col, type_col, target_col):
    """ Return number of floors based on mean floor height of the building """

    # If building type is 'Residential' mean floor height is 3.5
    if row[type_col] == 'Residential':
      row[target_col] = row[height_col] / 3.5
    # If building type is 'Work' mean floor height is 4.5
    elif row[type_col] == 'Work':
      row[target_col] = row[height_col] / 4.5
    # In other cases mean floor height is 4.0
    else:
      row[target_col] = row[height_col] / 4.0
    return row

def calculateArea(df, geom_col, target_col):
    """ Calculate the area of input geometry into <target_col> """
    # Create empty column for Area
    df[target_col] = None

    # Calculate Area of the geometries
    df = df.apply(getArea, axis=1, geom_col=geom_col, target_col=target_col)

    # Return DataFrame
    return df

def getArea(row, geom_col, target_col):
    """ Return the area of geometry """
    row[target_col] = row[geom_col].area
    return row

def calculateEHP(df, time_window):
    """
    Calculate Estimated Human Presence (EHP).
    See chapter 3.3 in the article + chapter S2.3, Figure S2 and Table S3 in the supplementary materials.
    
    """
    # Time Window
    tw = time_window + 't'
    
    # Calculate (absolute) estimated human presence (aEHP) for selected time window  ==> [Relative Floor Area] * [Seasonal Factor Coefficient] * [Hour Factor H]
    df['aEHP %s' % tw] = df['RFA'] * df['SF'] * df[tw]

    # Create column for estimated human presence (EHP) that is a normalized aEHP (scale 0.0 - 1.0).
    df['EHP %s' % tw] = None

    # Group data by 'Site_ID'
    grouped = df.groupby('site_id')

    # Iterate over groups and normalize the values
    for key, values in grouped:

      # Get the indices of the values
      siteid_indices = values.index

      # Get the sum of 'aEHP' values within site
      sum_aEHP = values['aEHP %s' % tw].sum()
     
      # Normalize the aEHP values by 'Site_ID' for each time unit (scale 0.0 - 1.0) ==> RMP
      EHP = values['aEHP %s' % tw] / sum_aEHP

      # Assign values to 'RMP' columns
      df.loc[siteid_indices, 'EHP %s' % tw] = EHP.values

    return df
    

def calculateRMP(cdr, time_window):
    """
    Calculate Relative Mobile Phone data distribution (RMP) for each subunit within a given base station.
    See 3.4 in the article and S1 in the supplementary materials.
    
    Note: In the manuscript this part is done later than in this script for practical reasons. 
    
    """
    
    # Name of the mobile phone user counts per site_id in the table
    twm = time_window + 'm'
    
    # Normalize the Mobile Phone user counts to scale 0.0 - 1.0
    cdr['RMP %s' % twm] = cdr[twm] / cdr[twm].sum()
    return cdr
    

def calculateROP(df, time_window):
    """ 
    Calculate the Relative Observed Population (ROP) for each subunit within a base station. See Table E in S2. 
    See chapter 3.4 in the article + chapter S2.4 and Table S4 in the supplementary materials.
    
    """
    twt = time_window + 't'
    # Attribute name for normalized mobile phone users (RMP)
    twm = 'RMP %s' % time_window + 'm'
    # Calculate 'ROP' ==> 'EHP hh-hh' * 'RMP hh-hh' 
    df['ROP %s' % twt] = df['EHP %s' % twt] * df[twm]
    return df

def calculateZROP(df, time_window, grid_cell_id_col):
    """ 
    Sum Relative Observed Population (ROP) for each target zone, i.e. calculate ZROP. 
    See chapter 3.5 in the article + chapter S2.5 and Table S5 in the supplementary materials.
    """
    
    tw = time_window
    rop = 'ROP ' + tw + 't'
    
    # Group data by 'Grid cell id'
    grouped = df.groupby(grid_cell_id_col)

    # Create DataFrame for spatial units (e.g. a 100 m grid)
    ZROP_grid = pd.DataFrame()

    # Iterate over grid cells
    for key, values in grouped:

      # Sum all ROP features that belongs to the same 'Grid cell id'
      zrop = values[rop].sum()

      # Append to DataFrame
      ZROP_grid = ZROP_grid.append([[key, zrop]])

    # Set column names
    ZROP_grid.columns = [grid_cell_id_col, 'ZROP %s' % tw]

    # Change grid cell id to numeric if possible
    try:
        ZROP_grid[grid_cell_id_col] = ZROP_grid[grid_cell_id_col].astype(int)
    except ValueError:
        print("Warning: Could not convert the ZROP values to numeric.")
        pass

    return ZROP_grid

def saveToShape(input_df, grid_df, output_path, grid_cell_id, epsg_code):
    """ Save ZROP values in <input_df> as Shapefile defined in <grid_df> to <output_path> using projection in <epsg code> """
    
    # Join the data with grid GeoDataFrame
    geo = grid_df[[grid_cell_id, 'geometry']].merge(input_df, left_on=grid_cell_id, right_on='GRD_INSPIR', how='inner')
    
    # Re-project
    geo['geometry'] = geo['geometry'].to_crs(epsg=epsg_code)

    # Ensure that results is GeoDataFrame
    geo = gpd.GeoDataFrame(geo, geometry='geometry', crs=from_epsg(epsg_code))

    # Fill NaN values with 0
    geo = geo.fillna(value=0)

    # Save to disk
    geo.to_file(output_path)

    return geo

if __name__ == "__main__":
    geo = main()
