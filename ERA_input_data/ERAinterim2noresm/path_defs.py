import calendar
from pathlib import Path

dic_month = dict((k, v.lower()) for k, v in enumerate(calendar.month_abbr))
data_folder = Path('/cluster/work/users/sarambl/ERAinterim_data')
input_folder = data_folder / 'rawdata'
res_file = data_folder / 'resfile_f19_tn14_h0.nc'
res_file_T = data_folder / f'{res_file.name[:-3]}_T.nc'
vct_file = data_folder / 'vct_fixed.txt'
tmp_folder = data_folder / 'tmp_data'
out_folder = data_folder / 'outdata'
