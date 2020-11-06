
"""
Making date and datesec arrays in netcdf format that we can
add to the ERA-interim metdata file so it can be used to nudge NorESM

Required modifications before you run it:
Insert information about the starting date, the length of the year and the length of february
Remember also to change the time.units so it matches the first date of your file.

Run it like this: python date_info_day.py 2001-01-01

Author: Inger Helene Hafsahl Karset
"""

import sys
from netCDF4 import Dataset
from numpy import dtype

def date_info_day(date_str, infile):
    """
    Making date and datesec arrays in netcdf format that we can
    add to the ERA-interim metdata file so it can be used to nudge NorESM

    Required modifications before you run it:
    Insert information about the starting date, the length of the year and the length of february
    Remember also to change the time.units so it matches the first date of your file.

    Author: Sara Blichner modified Inger Helene Hafsahl Karset's code

    :param date_str:
    :param infile:
    :param outfile:
    :return:
    """
    #date_str = str(sys.argv[1])
    #infile = './' + date_str + '.nc'

    # prepare date
    year,mon,day = date_str.split('-')
    year_num = int(float(year))
    mon_num = int(float(mon))
    day_num = int(float(day))


    datesec_calc = []
    val_pr_day = 4
    secstep = 86400/val_pr_day
    sec = [0, 1*secstep, 2*secstep, 3*secstep]
    for j in sec:
        datesec_calc.append(j)

    # Open a netCDF file for appending:
    ncfile = Dataset(infile,'a')
    #time_in = ncfile.variables['time'][:]
    #ncfile = Dataset('date_datesec' + date + '.nc','w')

    # Create the variable (4 byte integer in this case)
    # first argument is name of variable, second is datatype, third is
    # a tuple with the names of dimensions.
    date_str = ncfile.createVariable('date',dtype('int32').char,('time'))
    datesec = ncfile.createVariable('datesec',dtype('int32').char,('time'))

    # Write data to variable:
    date_str[:] = year_num*10000+mon_num*100+day_num
    datesec[:] = datesec_calc

    # Add attributes to the variables:
    date_str.long_name = 'current date (YYYYMMDD)'
    datesec.long_name = 'current seconds of current date'

    # close the file.
    ncfile.close()
    return
