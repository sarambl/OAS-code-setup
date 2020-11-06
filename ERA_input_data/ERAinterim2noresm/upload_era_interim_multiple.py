#!/usr/bin/env python
#
# (C) Copyright 2012-2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.
#
"""
Edited from ECMWF script for retrieving ERA-Interrim data.
Edits by Sara Blichner <s.m.blichner@geo.uio.no>
"""

import calendar
import sys
from calendar import monthrange

import numpy as np
from ecmwfapi import ECMWFDataServer

from ERA_input_data.ERAinterim2noresm.path_defs import input_folder

dic_month = dict((k, v.lower()) for k, v in enumerate(calendar.month_abbr))
# To run this example, you need an API key 
# available from https://api.ecmwf.int/v1/key/

server = ECMWFDataServer()


outfolder = input_folder


def import_years(start_y, end_y, fieldtype='3D'):
    """
    Download years
    :param start_y:
    :param end_y:
    :param fieldtype:
    :return:
    """
    import_monthly_files([start_y, 1], [end_y, 12], fieldtype=fieldtype)
    return


def import_monthly_files(start_yearmonth, end_yearmonth, fieldtype='3D'):
    """
    Import monthly files
    :param start_yearmonth:
    :param end_yearmonth:
    :param fieldtype: 3D is U,V,T etc, surf is PS etc
    :return:
    """
    start_y = start_yearmonth[0]
    start_m = start_yearmonth[1]
    end_y = end_yearmonth[0]
    end_m = end_yearmonth[1]
    for y in np.arange(start_y, end_y + 1):
        if y == start_y:
            if start_y == end_y:
                ran = np.arange(start_m, end_m + 1)
            else:
                ran = np.arange(start_m, 13)
        elif y == end_y:
            ran = np.arange(1, end_m + 1)
        else:
            ran = np.arange(1, 13)
        for mon in ran:
            get_monthly_file(y, mon, fieldtype=fieldtype)
    return


def get_monthly_file(year, month, fieldtype='3D'):
    """
    Download file for month
    :param year:
    :param month:
    :param fieldtype:
    :return:
    """
    start_date = int(year * 1e4 + month * 1e2 + 1)
    mon_len = monthrange(year, month)[1]
    end_date = int(year * 1e4 + month * 1e2 + mon_len)
    date = str(start_date) + '/to/' + str(end_date)

    # ofn = (dic_month[month]+str(year)+f'{fieldtype}.nc')
    ofn = (dic_month[month] + str(year) + f'{fieldtype}.grb')
    output_fname = outfolder / ofn  # 'from_%s_to_%s.grb'%(start_date, end_date)
    input_dic = get_input_dic(date, output_fname, fieldtype=fieldtype)
    print(input_dic)
    server_retrieval(input_dic)
    return


def get_input_dic(date, output_fname, fieldtype='3D'):
    """
    Makes request dictionary.
    :param date:
    :param output_fname:
    :param fieldtype:
    :return:
    """
    if fieldtype == '3D':
        levtype = "ml"
        param = "130/131/132/133"
    elif fieldtype == 'surf':
        levtype = 'sfc'
        param = "134"
    else:
        sys.exit(f'Fieldtype "{fieldtype}" not recognized. Aborting. ')
    input_dic = \
        {
            'dataset': "interim",
            'date': "%s" % date,
            'time': "00/06/12/18",
            'step': "0",
            'stream': "oper",
            'levtype': levtype,
            'levelist': "all",
            'type': "an",
            'class': "ei",
            'grid': "128",
            'param': param,
            'target': str(output_fname),
            # 'format': "netcdf"  # test
        }

    return input_dic


def server_retrieval(input_dic):
    server.retrieve(input_dic)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Correct usage: python upload_era_interim_multiple.py <star_year> <end_year>')
        sys.exit()
    if len(sys.argv) > 3:
        field_type = sys.argv[3]
    else:
        field_type = '3D'
    import_years(int(sys.argv[1]), int(sys.argv[2]), fieldtype=field_type)
