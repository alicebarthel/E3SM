import subprocess
import re
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from netCDF4 import date2num, Dataset
import numpy as np
from datetime import datetime
import gzip, io, os, glob
import argparse, sys

verbose = False
savedata = False

def main():
    parser = argparse.ArgumentParser(description='Run a check on CPL vs AM sea ice budget(s).')
    parser.add_argument('-d','--dir', type=str, required=True, help='Path to simulations directory (where simulation subdir resides)')
    parser.add_argument('-r','--run', type=str, required=True, help='Simulation name (e.g. date.version.mesh.machine)')
    parser.add_argument('-p','--period', type=int, nargs=2, required=True, help='Period to analyse. Two arguments expected: yearstart yearend')
    parser.add_argument('-hem','--hemisphere', type=int, nargs='+', default=0,  help='Hemisphere(s) to run the check over: 0=GLOBAL; 1=NH; 2=SH')
    parser.add_argument('-b','--budget', type=str, nargs='+', help='Type(s) of budget check to run. Options: h=heat, fw=freshwater, s=salt, a=area')
    parser.add_argument('-rt','--reltolerance', type=float, default=1e-4, help='Relative tolerance for the error check')
    parser.add_argument('-at','--abstolerance', type=float, default=1e-10, help='Absolute tolerance for the error check')
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help='Verbose option to write details of the process')
    parser.add_argument('-csv', "--savedata", action="store_true", default=False, help='Save the data structures we produce for the check to local csv files')
    args = parser.parse_args()
    if args.verbose:
        print("verbose mode turned on")
        verbose = True
    if args.savedata:
        print("saving data mode on")
        savedata = True
        
    print(args)
    compare_budgets_seaice(args.dir, args.run, args.budget, args.period[0], args.period[1], args.reltolerance, args.abstolerance, args.hemisphere)




def compare_budgets_seaice(modeldir, runname, budgetlist, ystart, yend, rtol, atol, hemislist):
    '''Compares sea ice flux terms from a) the coupler log, b) the conservation check Analysis Member (assuming monthly values).'''
    
    checksfailed = 0
# Checks on input and paths
    _check_timeperiod(ystart, yend)
    directory = '{}/{}/run/'.format(modeldir, runname)
    assert(os.path.isdir(directory)), "Directory not found: %s"%directory
    _check_AMfiles_exist(modeldir, runname, ystart, yend)
    assert(all([hem in [0,1,2] for hem in hemislist])), "Expected options for hemisphere are 0 = global, 1 = Nh, 2 = SH"
    assert(all([b in ['h', 's', 'fw', 'a'] for b in budgetlist])), "Expected options for budget are h, s, fw, a"
    assert(rtol <= 1 and rtol >= 0), "Relative tolerance is expected to be between 0 and 1"
    assert(atol >= 0), "Absolute tolerance is expected to be positive"
    
# extracting data
    for b in budgetlist:
        _verboseprint('Processing budget = %s'%b)
        dictbudget = _return_dict_for_budget(b)

        for hemis in hemislist:
            _verboseprint('Processing hemis = %s'%hemis)
            CPLdf = _extract_cpldata_fromdir(directory, dictbudget, ystart, yend, hemis)
            AMdf = _build_df_AM_SI(modeldir,runname,dictbudget, ystart, yend, hemis)

            CPLdf, AMdf = _check_adjust_shape(CPLdf, AMdf, dictbudget)
            if (savedata):
                CPLdf.to_csv('CPLdf.csv', header=True, index=True)
                AMdf.to_csv('AMdf.csv', header=True, index=True)

            ### Actual check ###
            check = np.allclose(CPLdf, AMdf, rtol, atol, equal_nan=False)
            print ("Are the CPL and AM %s budget arrays for %s equal within the tolerance: \t"%(dictbudget['Fullname'], _return_str_forhemis(hemis)), check)

            if not(check):
                checksfailed = checksfailed + 1
                print('------------------------------------------------------' )
                print('%12s  %6s    %8s   %s' %('Variable', 'pass?', 'max diff', 'time index of fails' ))
                print('------------------------------------------------------' )
                for ii in np.arange(len(dictbudget['var_header_am'])):
                    varcheck = np.allclose(CPLdf[dictbudget['var_header_am'][ii]], AMdf[dictbudget['var_header_am'][ii]], rtol, atol,  equal_nan=False)
                    if varcheck:
                        print('%12s  %6s' %(dictbudget['var_header_am'][ii],varcheck))
                    else: 
                        print('%12s  %6s    %.2E  ' %(dictbudget['var_header_am'][ii],varcheck, max(abs(CPLdf[dictbudget['var_header_am'][ii]]-AMdf[dictbudget['var_header_am'][ii]]))), list(np.where(np.logical_not(np.isclose(CPLdf[dictbudget['var_header_am'][ii]], AMdf[dictbudget['var_header_am'][ii]], rtol, atol,  equal_nan=False)))[0]))
                _plot_fluxtimeseries(CPLdf, AMdf, dictbudget,runname, ystart, yend,  hemis)
    
    assert(checksfailed == 0), "Sea-ice budget conservation checks failed: %s"%checksfailed

    
# Helper Functions

def _extract_cpldata_fromdir(dir_name, dictbudget, ystart, yend, hemis):
    list_of_files = sorted( filter( os.path.isfile,
                            glob.glob(dir_name + 'cpl.log*') ) )
    if (len(list_of_files) == 0):
        raise ValueError('no cpl file matching expected pattern in the given directory')
        #os.exit(1)
    if (len(list_of_files) == 1):
        #print(list_of_files[0])
        df = _extract_fromcpl_budgetterms(list_of_files[0], dictbudget, ystart, hemis)
    if (len(list_of_files) > 1):
        dflist = [None] *  len(list_of_files) # I use this instead of append because it fails when appending several elements at a time
        for ii in np.arange(0, len(list_of_files)):
            _verboseprint("Extracting from cpl file: %s"%list_of_files[ii])
            dflist[ii] = _extract_fromcpl_budgetterms(list_of_files[ii], dictbudget, ystart, hemis)
            #print('v2 time range ', dflist[ii].index.min(), dflist[ii].index.max())
        df = pd.concat(dflist, verify_integrity=True) #replace with in-place?
        #print(df.index.min(), df.index.max())
        df = df[(df.index > ((ystart-1)*12)) & (df.index <= (yend)*12-1) ] # sub-select the time period we asked for

        if np.shape(df)[0] == 0:
            raise ValueError("the extracted CPL data does not overlap with the time period requested ")
        elif (np.shape(df)[0] < (ystart-yend+1)*12-1):
            raise ValueError("the extracted CPL data only partially covers the time period requested ") # this could a warning if we wanted to compare the overlap

        assert(np.shape(df)[0]==(yend-ystart+1)*12-1), "the extracted CPL data does not have the expected shape"
        _verboseprint("the extracted CPL data covers the requested time period: %s"%(np.shape(df)[0]==(yend-ystart+1)*12-1))
    
    return df


def _extract_fromcpl_budgetterms(filename, dictbudget, ystart, hemis): 
### Extract all sea-ice budget terms into hemisphereic dataframe (global, nh, sh)  -- format agnostic
    file_content = _open_cpl_file(filename).read()
    blocks = re.findall(
    dictbudget['regex_header'],
    file_content,
    re.DOTALL
    )
    
    mylist = []
    for block in blocks:
        time = re.findall("date\s\=\s\s*(-?\d+)", block)
        time = date2num(datetime.strptime('{:08d}'.format(eval(time[0])), '%Y%m%d'), units = 'months since 0001-01-01', calendar = '360_day')
        _verboseprint('t = %s'%time)
        for var, varc in zip(dictbudget['var_header_cpl'], dictbudget['budget_terms']):
            if (varc == 'SUM'):
                for row in re.findall(dictbudget['regex_SUM']%varc, block):
                    mylist.append([time] + [var] + [_extract_value_forhemis(row, hemis)])

            else:
                for row in re.findall(dictbudget['regex_var']%varc, block):
                    mylist.append([time] + [var] + [_extract_value_forhemis(row, hemis)])


    cpl_hterms_ice = pd.DataFrame(mylist, columns = ['time', 'var', 'ice %s'%_return_str_forhemis(hemis)])
    cpl_hterms_ice = cpl_hterms_ice.pivot(index="time", columns="var", values="ice %s"%_return_str_forhemis(hemis))

    return cpl_hterms_ice

    
def _build_df_AM_SI(modeldir,runname,dictbudget, ystart, yend, hemis):
# Extract Hemispheric SEA-ICE variables from the SEA-ICE conservation Analysis Member, into 1 hemispheric dataframe
    reftime = date2num(datetime.strptime('{:04d}'.format(ystart), '%Y'), units = 'months since 0001-01-01', calendar = '360_day')
    # Build the data_array that will make it into a dataframe
    var_header = dictbudget['var_header_am']
    var_consCheck = dictbudget['var_consCheck']
    lentime = (yend-ystart+1)*12 # assuming monthly data and full years
    time = np.arange(0, lentime) + reftime
    dataarray = np.zeros([lentime, len(var_consCheck)])
    for ii  in np.arange(len(var_consCheck)): 
        dataarray[:, ii] = _extract_fromConsAM(modeldir,runname, ystart, yend, var_consCheck[ii], hemis)

    # Make a dataframe with matching AM data
    AMdf = pd.DataFrame(data = dataarray, index=time, columns=var_header)
    
    return AMdf


def _extract_fromConsAM(modeldir,runname, ystart, yend, var, hemis):
### Extract the requested variable from the SEA-ICE conservation Analysis Member
# Assumes: monthly data; hard coded SI conservation check file names; skips time stamp zero (first year is 11 expected length)
# (Looping is necessary due to 0-11 structure of Time in AM.nc -- cannot automatically concatenate)
    fluxvar = np.zeros((yend-ystart+1)*12) # assuming monthly data and full years
    for year in np.arange(ystart, yend+1):
        modelfile = '{}/{}/run/{}.mpassi.hist.am.conservationCheck.{:04d}.nc'.format(modeldir, runname, runname,year)
        time = 12*(year-ystart)
        f = Dataset(modelfile, mode='r')
        _verboseprint('Extracting var: %s'%var +' hemis: %s'%hemis)
        _verboseprint('time: %s, %s'%(time, time+11))
        ttime=max(time,1)
        if var == 'netShortwaveFlux':
                fluxvar[ttime: time+12] = f.variables['energyConsAbsorbedShortwaveFlux'][:, hemis] + f.variables['energyConsOceanShortwaveFlux'][:, hemis] # 1 is for NH, add others
        else: 
                fluxvar[ttime: time+12] = f.variables[var][:,  hemis]
        
    return fluxvar


def _check_adjust_shape(CPLdf, AMdf, dictbudget):
    if (np.shape(CPLdf)==np.shape(AMdf)):
        _verboseprint('Coupler and Analysis Member data structures have the same shape')
    if (np.shape(CPLdf)[1] > np.shape(AMdf)[1]):
        CPLdf=CPLdf[dictbudget['var_header_am']] #subselecting variables (e.g. dropping runoff and h2otemp)

    if (np.shape(AMdf)[0] == np.shape(CPLdf)[0]+1 ):
        AMdf=AMdf[dictbudget['var_header_am']].drop(min(AMdf.index)) #dropping first index in AM to match CPL (which does not produce 1st time)
    else: 
        _verboseprint('Shape of data structures: ',np.shape(CPLdf),np.shape(AMdf))

    assert(np.shape(CPLdf) == np.shape(AMdf)), "Coupler and Analysis Member data structures have incompatible shapes"
    return CPLdf, AMdf


def _open_cpl_file(filename):
    root, extension = os.path.splitext(filename)
    #do I need a tar.gz case?
    if (extension == '.gz'):
        in_file = io.TextIOWrapper(gzip.open(filename, 'rb'), encoding='utf-8')
    else:
        in_file = open(filename, 'r')
    return in_file

def _verboseprint(string):
    if verbose:
        print(string)


def _return_dict_for_budget(budgetname):
    # Creating a Dictionary for each budget type
    if budgetname == 'h':
        Dict = {'Name': budgetname, 'Fullname': 'heat',
    'var_header_am': ['NetHeat','Frazil','Sensible','OceanHeat', 'LongwaveUp','LongwaveDown','NetSW','Latent','SnowHeat'],
    'var_consCheck' : ['netEnergyFlux', 'energyConsFreezingPotential', 'energyConsSensibleHeatFlux', 'energyConsOceanHeatFlux', 'energyConsLongwaveUp', 'energyConsLongwaveDown', 'netShortwaveFlux', 'energyConsLatentHeat', 'energyConsSnowfallHeat'],
    'var_header_cpl' : ['Frazil','OceanHeat','NetSW','LongwaveDown', 'LongwaveUp','Latent','SnowHeat','hiroff','Sensible','hh2otemp','NetHeat'],
    'budget_terms' : ['hfreeze', 'hmelt', 'hnetsw', 'hlwdn', 'hlwup', 'hlatvap', 'hlatfus', 'hiroff', 'hsen', 'hh2otemp', 'SUM'],
    'regex_header': r'NET\sHEAT\sBUDGET\s\(W\/m2\):\speriod\s=\s\smonthly(.*?)\n\s\s\n',
    'regex_SUM' : r"\s\*%s\*\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)",
    'regex_var' : r"\s%s\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)"
               }
    elif budgetname == 'fw':
        raise ValueError('Freshwater budget not implemented yet, this is a placeholder')
    elif budgetname == 's':
        raise ValueError('Salt budget not implemented yet, this is a placeholder')
    elif budgetname == 'a':
        raise ValueError('Area budget not implemented yet, this is a placeholder')        
    else: 
        raise ValueError('Budget name unkown; Working option is h, (future options are fw, s, a)')
    _verboseprint("\ncreated dictionary for the budget requested. ")
    
    return Dict

def _check_timeperiod(ystart, yend):
    for year in [ystart, yend]:
        if (type(year) != int):
            print('Warning: Year input is not an int, converted to int')
        val = int(year)
        if (val < 1): 
            raise ValueError("Year inputs must be a positive integer ( > 1)")
            
        if (val > 2100): 
            raise ValueError("Year inputs must be an integer below 2100")

    if int(yend) < int(ystart):
        raise ValueError('yend must be greater than or equal to ystart')


def _check_AMfiles_exist(modeldir, runname,ystart, yend):
    missingfiles=[]
    for year in np.arange(ystart, yend):
        modelfile = '{}/{}/run/{}.mpassi.hist.am.conservationCheck.{:04d}.nc'.format(modeldir, runname, runname,year)
        try:
            f = open(modelfile, "r")
            f.close()
        except IOError:
            missingfiles.append(year)
    if len(missingfiles) > 0:
        print("%s files missing in directory: "%len(missingfiles))
        print('{}/{}/run/{}.mpassi.hist.am.conservationCheck.YEAR.nc'.format(modeldir, runname, runname))
        print('for YEAR in: ')
        print(missingfiles)
        raise IOError("Files cannot be opened (or do not exist).")

def _plot_fluxtimeseries(CPLdf, AMdf, dictbudget,runname, ystart, yend, hemis):
    var_header = dictbudget['var_header_am']
    fig, ax = plt.subplots(len(var_header),1, figsize=(16,16))
    for ii in np.arange(len(var_header)):
        ax[ii].plot(CPLdf[var_header[ii]], 'k--')
        ax[ii].plot(AMdf[var_header[ii]], 'r:')
        ax[ii].set_title(var_header[ii], y=1.0, pad=-14)
    ax[-1].set_xlabel('Months since 0001-01-01')
    fig.suptitle('%s : %s budget %s  year %s to %s'%(runname, dictbudget['Fullname'], _return_str_forhemis(hemis), ystart, yend))
    fig.tight_layout()
    fig.savefig('seaice-ocean_fluxcomparison_cpl_am.%s_%sbudget_%s_%s-%s.png'%(runname, dictbudget['Name'], _return_str_forhemis(hemis), ystart, yend), dpi= 150, format='png')
    return 0


def _extract_value_forhemis(row, hemis):
    if hemis == 0:
        return (eval(row[4]) + eval(row[5]))
    elif hemis == 1:
        return eval(row[4])
    elif hemis == 2:
        return eval(row[5])
    else:
        raise ValueError("Hemis should be within 0,1,2 - this error should be caught earlier!")
        
def _return_str_forhemis(hemis):
    if hemis == 0:
        return "GL"
    elif hemis == 1:
        return "NH"
    elif hemis == 2:
        return "SH"
    else:
        raise ValueError("Hemis should be within 0,1,2 - this error should be caught earlier!")



if __name__ == "__main__":
    main()
