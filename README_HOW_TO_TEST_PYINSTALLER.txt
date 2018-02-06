***

BUILD CODE
normal py36 environment, conda3

WORKING PART TO CREATE A UNIQUE FILE
pyinstaller --noupx --onefile sherpa.py
WORKING PART TO CREATE INSTALLATION DIR WITH LOT OF FILES
pyinstaller --noupx --onedir sherpa.py

but in runtime 'qt' is not working...fix needed:
https://github.com/pyinstaller/pyinstaller/commit/082078e30aff8f5b8f9a547191066d8b0f1dbb7e
https://github.com/pyinstaller/pyinstaller/commit/59a233013cf6cdc46a67f0d98a995ca65ba7613a

TESTING PART
pyinstaller --upx-dir=D:/WORK/projects/1_urbIam/1_CODE_MATLAB/SHERPA/SHERPA-GITHUB/upx --onedir sherpa.py
	
***

TEST CODE: PYTHON
python sherpa.py 1 input/base_emissions/BC_emi_PM25_Y.nc input/London_region.nc input/user_reduction_snap7.txt input/base_concentrations/BC_conc_PM25_Y.nc input/source_receptors/SR_PM25_Y.nc output/

python sherpa.py 1 input/base_emissions/BC_emi_PM10_Y.nc input/London_region.nc input/user_reduction_snap7.txt input/base_concentrations/BC_conc_PM10_Y.nc input/source_receptors/SR_PM10_Y.nc output/

python sherpa.py 1 input/base_emissions/BC_emi_NO2_Y.nc input/London_region.nc input/user_reduction_snap7.txt input/base_concentrations/BC_conc_NO2_NO2eq_Y_mgm3.nc input/source_receptors/SR_NO2eq_Y.nc output/

python sherpa.py 8 input/impacts/healthbl_nc.nc output/ config/sharedvariables.json input/base_concentrations/BC_conc_PM25_Y.nc

python sherpa.py 6 input/base_emissions/BC_emi_PM25_Y.nc input/EMI_RED_ATLAS_NUTS2.nc 45 10 input/user_reduction_snap7.txt input/base_concentrations/BC_conc_PM25_Y.nc input/source_receptors/SR_PM25_Y.nc output/

python sherpa.py 6 input/base_emissions/BC_emi_PM10_Y.nc input/EMI_RED_ATLAS_NUTS2.nc 45 10 input/user_reduction_snap7.txt input/base_concentrations/BC_conc_PM10_Y.nc input/source_receptors/SR_PM10_Y.nc output/

python sherpa.py 6 input/base_emissions/BC_emi_NO2_Y.nc input/EMI_RED_ATLAS_NUTS2.nc 45 10 input/user_reduction_snap7.txt input/base_concentrations/BC_conc_NO2_NO2eq_Y_mgm3.nc input/source_receptors/SR_NO2eq_Y.nc output/

python sherpa.py 7 input/base_emissions/BC_emi_PM25_Y.nc input/base_concentrations/BC_conc_PM25_Y.nc input/pDUST-pSALT/ input/source_receptors/SR_PM25_Y.nc input/selection/gridnew/fua/ input/selection/gridnew/nuts/ input/selection/gridnew/ input/AM_targets.txt output/ input/sherpa_icon_name_256x256.png fua PM25

python sherpa.py 31 input/base_emissions/BC_emi_PM25_Y.nc input/London_region.nc input/user_reduction_mod3b.txt input/base_concentrations/BC_conc_PM25_Y.nc input/source_receptors/SR_PM25_Y.nc output/

python sherpa.py 32 input/base_emissions/BC_emi_PM25_Y.nc input/London_region.nc input/user_reduction_mod3b.txt input/base_concentrations/BC_conc_PM25_Y.nc input/source_receptors/SR_PM25_Y.nc output/

python sherpa.py 31 input/base_emissions/BC_emi_PM10_Y.nc input/London_region.nc input/user_reduction_mod3b.txt input/base_concentrations/BC_conc_PM10_Y.nc input/source_receptors/SR_PM10_Y.nc output/

python sherpa.py 32 input/base_emissions/BC_emi_PM10_Y.nc input/London_region.nc input/user_reduction_mod3b.txt input/base_concentrations/BC_conc_PM10_Y.nc input/source_receptors/SR_PM10_Y.nc output/

python sherpa.py 31 input/base_emissions/BC_emi_NO2_Y.nc input/Essen_region.nc input/user_reduction_mod3b.txt input/base_concentrations/BC_conc_NO2_NO2eq_Y_mgm3.nc input/source_receptors/SR_NO2eq_Y.nc output/

python sherpa.py 32 input/base_emissions/BC_emi_NO2_Y.nc input/Essen_region.nc input/user_reduction_mod3b.txt input/base_concentrations/BC_conc_NO2_NO2eq_Y_mgm3.nc input/source_receptors/SR_NO2eq_Y.nc output/
***

TEST CODE: EXE
as previous, but remove 'python' and run 'sherpa.exe'

***
