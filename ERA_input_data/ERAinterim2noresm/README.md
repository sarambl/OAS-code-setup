# To make noresm input for nudging to ERA-Interim
- edit paths in [path_defs.py](path_defs.py), i.e. define where to download data etc. 
- run [download_and_compute.py](download_and_compute.py) as
```shell script
python download_and_compute.py <start_year> <end_year>

```
### OR
First run [upload_era_interim_multiple.py](upload_era_interim_multiple.py) and then [conv_ERA_interim.py](conv_ERA_interim.py)