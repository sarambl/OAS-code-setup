"""
Author: Sara Blichner <s.m.blichner@geo.uio.no>
Based on conv_ERA-interim.sh:
# Author: Matthias Hummel <hummel@geo.uio.no>
# and Inger Helene H. Karset <i.h.h.karset@geo.uio.no>
# Date: 08.09.2016
# Modified by Moa Sporre 14.11.2017


"""
import os
from glob import glob
from pathlib import Path
from subprocess import run
import sys
from netCDF4 import Dataset
from numpy import dtype
# %%
import xarray as xr

from ERA_input_data.ERAinterim2noresm.path_defs import dic_month, data_folder, input_folder, res_file, res_file_T, vct_file, tmp_folder, out_folder

out_folder.mkdir(exist_ok=True)
tmp_folder.mkdir(exist_ok=True)


def get_daily_fl(mon_nr, year, prefix='', folder=tmp_folder):
    """
    Gets list of files for daily files
    :param mon_nr: month in format number
    :param year: year
    :param prefix: prefix before the filename
    :param folder: where to look.
    :return:
    """
    fn = fn_raw_nc(mon_nr, year, '3D')
    tmpfn = folder / f'{prefix}{fn}'
    fls = glob(f'{str(tmpfn)[:-3]}_day*.nc')
    return fls


def grb2nc(fn, output_folder=input_folder):
    """
    Compute netcdf file from grb files.
    :param fn: grb file to be converted.
    :param output_folder: where to put the nc files.
    :return:
    """
    fn_out = output_folder / f'{Path(fn).name[:-3]}nc'
    cdo_com = f'cdo -s -t ecmwf -f nc copy {fn} {fn_out}'
    if not fn_out.is_file():
        run(cdo_com, shell=True)
    else:
        print(f'File found, using existing file {fn_out}')
    return


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


def change_unit_time_and_save_final(year, mon_num, fn, day_num=None):
    """
    Unit of time needs to be changed and then final file is saved.
    Original file in "hours since", but we need "days since"
    :param year:
    :param mon_num:
    :param fn:
    :param day_num:
    :return:
    """
    ds = xr.open_dataset(fn)
    if 'units' in ds['time'].attrs.keys():
        if 'days' in str(ds['time'].units):
            print(f'unit already in days: {ds["time"].units}')
            return
    if day_num is None:
        day_num = fn[-5:-3]
    mon_nr = '%02d' % mon_num
    date = f'{year}-{mon_nr}-{day_num}'

    com1 = f'ncap2 -s "time@units=\\"days since {date} 00:00:00\\"" {fn}'
    print(com1)
    run(com1, shell=True)

    final_out_fn = out_folder / f'{date}.nc'
    com2 = f'ncap2 -s time=time/24  {fn} {final_out_fn}'
    print(com2)
    run(com2, shell=True)
    date_info_day(date, final_out_fn)

    return


def fn_raw(month_nr, year, fieldtype):
    """
    Filename raw input file
    :param month_nr:
    :param year:
    :param fieldtype:
    :return:
    """
    return f'{dic_month[month_nr]}{year}{fieldtype}.grb'


def fn_raw_nc(month_nr, year, fieldtype):
    """
    Filename raw in nc format
    :param month_nr:
    :param year:
    :param fieldtype:
    :return:
    """
    return f'{dic_month[month_nr]}{year}{fieldtype}.nc'


def get_fl_raw(fieldtype, year):
    """
    File list of raw files in year.
    :param fieldtype: "3D" og "surf"
    :param year:
    :return:
    """
    fl_3D = glob(str(input_folder) + f'/*{year}{fieldtype}.nc')
    fl_check = [fn_raw_nc(i, year, fieldtype) for i in range(1, 13)]
    fl_check2 = [Path(f).name for f in fl_3D]
    if set(fl_check) != set(fl_check2):
        ms_f = [f for f in fl_check if f not in fl_check2]
        # print(f'Lacking files for year {year}: \n {ms_f}. Want to proceed? (y/n)')
        if __name__ == '__main__':
            ans = input(f'Lacking files for year {year}: \n {ms_f}. Want to proceed? (y/n) ')
            if ans.strip() in ['n', 'N']:
                sys.exit()
    return fl_3D


def make_vct_file(res_fn=res_file, outfn=vct_file):
    """
    Make vct file for re-formatting lev.
    :param res_fn:
    :param outfn:
    :return:
    """

    res_f = xr.open_dataset(res_fn)
    hyai = 'hyai'
    hybi = 'hybi'
    hyai_da = res_f[hyai].values
    hybi_da = res_f[hybi].values
    outfile = open(str(outfn), 'w')

    info_line = '#   k         vct_a(k) [Pa]             vct_b(k) []\n'
    outfile.write(info_line)
    # hyai needs to be multiplied by 1e5
    for i in range(0, len(hyai_da)):
        outfile.write('%4d %25.17f %25.17f \n' % (i, 1e5 * hyai_da[i], hybi_da[i]))
    outfile.close()
    return outfn


def input_fn_from_ym(year, month):
    """
    filename from year/month input.
    :param year:
    :param month:
    :return:
    """
    return f'{dic_month[month]}{year}.nc'


# %%
def main(year):
    fl_3D = get_fl_raw('3D', year)
    print(fl_3D)
    # %%
    # Get resolution:
    rf = xr.open_dataset(res_file)
    xinc = float(rf['lon'][1] - rf['lon'][0])  # .values[0]
    yinc = float(rf['lat'][1] - rf['lat'][0])
    print(f'Resolution file found: dx = {xinc}, dy={yinc}')
    # %%
    # cp surface to new 3D files:
    daily_files = []
    for mon_nr in range(1, 13):
        fn_3D_grb = input_folder / fn_raw(mon_nr, year, '3D')
        fn_surf_grb = input_folder / fn_raw(mon_nr, year, 'surf')
        # if not fn
        if not fn_3D_grb.is_file():
            print(f'Could not find file {fn_3D_grb}')
            continue
        if not fn_surf_grb.is_file():
            print(f'Could not find file {fn_surf_grb}')
            continue
        print(f'converting {fn_3D_grb} to nc')
        grb2nc(fn_3D_grb)
        print(f'converting {fn_surf_grb} to nc')
        grb2nc(fn_surf_grb)
    print('Done converting to nc')
    print('Starting to merge files:')
    for mon_nr in range(1, 13):
        fn_3D = input_folder / fn_raw_nc(mon_nr, year, '3D')
        fn_surf = input_folder / fn_raw_nc(mon_nr, year, 'surf')
        if not fn_3D.is_file():
            print(f'Could not find file {fn_3D}')
            continue
        if not fn_surf.is_file():
            print(f'Could not find file {fn_surf}')
            continue
        tmpfn = tmp_folder / fn_3D.name

        fls = glob(f'{str(tmpfn)[:-3]}_day*.nc')
        if len(fls) >= 28:
            print('Daily files already computed')
            daily_files = daily_files + fls
            continue
        cdo_comm = f'cdo -s merge {str(fn_surf)} {str(fn_3D)} {str(tmpfn)}'
        print(cdo_comm)

        run(cdo_comm, shell=True)
        cdo_comm2 = f'cdo -s splitday {str(tmpfn)} {str(tmpfn)[:-3]}_day'
        print(cdo_comm2)
        run(cdo_comm2, shell=True)
        os.remove(tmpfn)
        fls = glob(f'{str(tmpfn)[:-3]}_day*.nc')
        daily_files = daily_files + fls

    # make vct file:
    vct_filen = make_vct_file(res_fn=res_file)
    # %%
    daily_files.sort()
    # %%
    print(f'Remapping level:')
    for fd in daily_files:
        pd = Path(fd)
        pout = pd.parent / f'vert_{pd.name}'
        if pout.is_file():
            continue
        comm_cdo = f'cdo -s --no_warnings remapeta,{vct_filen} -chname,SP,APS -selname,U,V,T,Q,SP {fd} {pout}'
        print(comm_cdo)
        run(comm_cdo, shell=True)
    # %%
    print('Done remapping eta/lev')
    print('Change horizontal resolution...')
    fd = daily_files[0]
    pd = Path(fd)
    vert_f = pd.parent / f'vert_{pd.name}'
    # change horizontal resolution
    print('creating weights file....')
    cdo_comm = f'cdo -s selname,T {res_file} {res_file_T}'
    p_weights = data_folder / 'weights.nc'
    cdo_comm2 = f'cdo -s genbil,{res_file_T} {vert_f} {p_weights}'
    print(cdo_comm)
    run(cdo_comm, shell=True)
    print(cdo_comm2)
    run(cdo_comm2, shell=True)

    for fd in daily_files:
        pd = Path(fd)
        v_fp = pd.parent / f'vert_{pd.name}'
        pout = v_fp.parent / f'horiz_{v_fp.name}'
        if pout.is_file():
            continue
        cdo_comm = f'cdo -s remap,{res_file_T},{p_weights} {v_fp} {pout}'
        print(cdo_comm)
        run(cdo_comm, shell=True)
        rn_comm = f'ncrename -v APS,PS {pout}'
        print(rn_comm)
        run(rn_comm, shell=True)

    for mon_nr in range(1, 13):
        fl = get_daily_fl(mon_nr, year, prefix='horiz_vert_', folder=tmp_folder)
        for fd in fl:
            change_unit_time_and_save_final(year, mon_nr, fd)

    return


# %%
if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.exit('Lacking input year \n Correct usage: '
                 'python conv_ERA_interim.py <year>')

    inyear = sys.argv[1]

    main(inyear)
