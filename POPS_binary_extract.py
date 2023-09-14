#!/Users/FarmerLab/anaconda3/envs/COMP/bin/python

# latest version duplicated from POPS_binary_extract-V5.py
# Ricardo L. Pena 9/14/23
# Colorado State University
# penar@colostate.edu
# Farmer Group

#import modules
import numpy as np
import pandas as pd
import struct
import time
from itertools import chain
import glob
from datetime import datetime
import os
import os.path

from argparse import ArgumentParser
from multiprocessing import Pool

### flatten list function
def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))

### extract binary function
def extract_binary(binary_filename):
    with open(binary_filename, mode='rb') as file: # b is important -> binary
        fileContent = file.read()
    lst = [fileContent]

    file_length = len(fileContent) # full length of file
    line_index_start = 0 # initialize the index to start at zero
    list_peakamplitude, list_timestamps = [], [] # empty the lists

    while file_length>line_index_start:
        line_index_end = line_index_start+12 # do this here to not have to do it again
        num_records, timestamp_val = struct.unpack('<Id',fileContent[line_index_start:line_index_end])
        num_elements = num_records*3
        line_record_end = line_index_end+(num_elements*4)

        # get all the data out of the binary record here and into a numpy array
        # might be a faster way to do this

        chunk_data = struct.unpack('I'*num_elements,fileContent[line_index_end:line_record_end])
        chunk_data = np.asarray(chunk_data).reshape(-1,3) # turn into an array and reshaped, could be done in a single step

        # make an empty array to do the timestamp shenanigans
        # dt = time between each record except for the first entry
        # first dt is the time since the timestamp_val, float time

        chunk_timestamp = np.empty(len(chunk_data),dtype='float64') # empty array with float precision
        chunk_timestamp = np.cumsum(chunk_data[:,2]/1e6) # cumulative sum of the dt values, divide by 1e6 is because dt is in microseconds
        chunk_timestamp += timestamp_val # adding the initial timestamp of the data record to get the true timestamp

        # list comprehension is faster than numpy
        chunk_timestamp = chunk_timestamp[:].tolist() # converting to list
        chunk_amplitude = chunk_data[:,0].tolist() # converting to a list

        # the two outputs PeakAmplitude and the Timestamp
        # these are lists of lists
        list_peakamplitude.append(chunk_amplitude)
        list_timestamps.append(chunk_timestamp)
        line_index_start = line_record_end # reset things

    # at the end of this the lists of lists need flattening
    t = time.time()
    PeakAmplitude = fast_flatten(list_peakamplitude)
    Timestamp = fast_flatten(list_timestamps)
    elapsed = time.time() - t
    del fileContent; del lst
    return PeakAmplitude,Timestamp

def input_arguments(): # arguments that are user defined
    # run script + -h to list current ports and see arguments.
    parser = ArgumentParser()
    parser.add_argument('--directory', type = str,
        help = 'Select the directory containing binary files or the directory containing subdirectories containing binaries.',
        default=".")
    parser.add_argument('--nbins', type = int,
        help = 'Select the desired number of bins (default = 16 bins).',
        default=16)
    parser.add_argument('--mie', type = str,
        help = 'Select the smoothed mie conversion table to interpolate bins from (default = Mie_scripps_1_41.csv).',
        default='./mie_conversion_tables/Mie_scripps_1_41.csv')
    parser.add_argument('--multiproc', action='store_true', help= 'call argument to enable binary processing with four processors.')
    parser.add_argument('--logmin', type = float,
        help = 'log of the amplitude corresponding to 120 nm in mie table (default = 0.871 for Scripps)',
        default=0.871)
    parser.add_argument('--logmax', type = float,
        help = 'log of the amplitude corresponding to 3000 nm in mie table (default = 4.067 for Scripps)',
        default=4.067)

    arguments = parser.parse_args()
    return arguments

### Mie Bin-creation functions
# Bin set-up and convertion based on Mie conversion table
def bin_headers(nbins, logmin, logmax):
    #replicating mie conversion excel sheet sent by Dr. Bryan Rainwater at Handix
    # The number of bins, logmin, and log max can change.
    upper_bin_edge = []
    lower_bin_edge = []
    for i in range(0,nbins):
        upper_n = i + 1
        lower_n = i

        l = logmin + ((logmax-logmin)/nbins)*lower_n
        u = logmin + ((logmax-logmin)/nbins)*upper_n
        # bin edges are reported as loggeed values; 10**bin_edge to get non-logged intensity/bin edges
        lower_bin_edge.append(l)
        upper_bin_edge.append(u)
    return lower_bin_edge, upper_bin_edge

def linear_interpolation(x_lower,y_lower,x_upper,y_upper,x_interp, logged):
    if logged == True: np.log10(x_interp)
    y_interp = y_lower + ((x_interp - x_lower)*((y_upper-y_lower)/(x_upper-x_lower)))
    return y_interp

def itensity_2_diameter(bin_range, Mie_map_amp, Mie_map_diam, Log_Intensity):
    raw_amplitude, interp_diameter = [], []
    for signal_amp in bin_range:
        try:
            IOR_raw_amp = np.log10(Mie_map_amp) # immediately take log to compare to logged bin edges(signal_amp)
            # bool = is amp value in mie table the log intensity of the amplitude or the raw amplitude?
            # True indicates that the mie table displays logged amplitudes or intensities.
            # False indicates that raw amplitudes/intensities are displayed in conversion table.
            if Log_Intensity == True: IOR_raw_amp = Mie_map_amp

            for i, amp in enumerate(IOR_raw_amp):
                    bottom = bool(signal_amp >= amp)
                    top = bool(signal_amp <=amp)

                    if ((bottom == False) and (top == True)):
                        upper_idx = i
                        lower_idx = i-1
                        break
            upper_amp, lower_amp =  Mie_map_amp[upper_idx], Mie_map_amp[lower_idx]
            upper_diam, lower_diam = Mie_map_diam[upper_idx], Mie_map_diam[lower_idx]

            # x is amp and y is diameter
            #if Log_Intensity == False: signal_amp = 10**signal
            signal = 10**signal_amp #the signal_amp is always logged, but log needs to be undone to correctly bin raw intensities which are not logged.
            signal_diam = linear_interpolation(lower_amp, lower_diam, upper_amp, upper_diam, signal, Log_Intensity)
            raw_amplitude.append(signal)
            interp_diameter.append(signal_diam)
            #print(signal, signal_diam)

        except Exception as e:
            print("outside of Mie amp range")
            print(e)
    return raw_amplitude, interp_diameter

# Calculate log intensity bin edges for desired bin number
def create_bins(nbins,logmin,logmax, mie_conv_table_intensity, mie_conv_table_diam, logged_amplitudes):
    l_bin, u_bin = bin_headers(nbins, logmin, logmax)

    upper_raw_amp, upper_interp_diam = itensity_2_diameter(u_bin,mie_conv_table_intensity, mie_conv_table_diam, logged_amplitudes)
    lower_raw_amp, lower_interp_diam = itensity_2_diameter(l_bin, mie_conv_table_intensity, mie_conv_table_diam,logged_amplitudes)
    mie_smooth_dict = {"Lower Amp": lower_raw_amp, "Upper Amp": upper_raw_amp,
                "Lower Bin Diameter": lower_interp_diam, "Upper Bin Diameter": upper_interp_diam}
    mie_smooth_df = pd.DataFrame(mie_smooth_dict)
    mie_smooth_df['Mean Diameter (nm)'] = mie_smooth_df[['Lower Bin Diameter','Upper Bin Diameter']].mean(axis=1)
    mie_smooth_df['Bin Header'] = round(mie_smooth_df['Lower Bin Diameter'],2).astype('str') + "_" + round(mie_smooth_df['Upper Bin Diameter'],2).astype('str')
    return mie_smooth_df

### Main working functions
def binary_2_csv(binary_mie_bins_info): # main working function
    # retrieve binary file to process and the bins and mie table to use
    binary_file = binary_mie_bins_info[0]
    mie_table = binary_mie_bins_info[1]
    bins = binary_mie_bins_info[2]

    file_csv_name = binary_file.split('.b')[0]+'_10Hz.csv'
    ### check if file already exists before extracting the binary again ###
    if (os.path.isfile(file_csv_name) == False):
        print('created {}'.format(file_csv_name))
        peaks, time_s = extract_binary(binary_file)

        # create a raw dataframe from binary data
        raw_peak_df = pd.DataFrame({'timestamp': time_s,'Peaks':peaks}) #.to_csv(file_csv_name)
        #raw_peak_df.to_csv(binary_file.split('.b')[0]+'_raw.csv') ### binning check!
        raw_peak_df['timestamp'] = raw_peak_df['timestamp'].apply(lambda x: datetime.utcfromtimestamp(x)) #recently updated
        timing_df = raw_peak_df.set_index('timestamp') # set timestamp as index

        # resampling the data to 10Hz
        try:
            df_10Hz = timing_df.groupby(pd.Grouper(freq='100ms'))['Peaks'].apply(list)
            #df_10Hz.to_csv(binary_file.split('.b')[0]+'_grouped-10Hz.csv') ### binning check!
            # bin and count the data now
            bin_counts = df_10Hz.apply(lambda x: pd.cut(x, bins=bins).value_counts())
            bin_counts.columns = mie_table['Bin Header']
            bin_counts.to_csv(file_csv_name)
            lst_p = [peaks, time_s, raw_peak_df, timing_df, df_10Hz, bin_counts]
            del peaks; del time_s; del raw_peak_df; del timing_df; del df_10Hz; del bin_counts
            del lst_p
        except: print("Fatal error in file: ",file_csv_name)
    else:
        print ('{} already exists!'.format(file_csv_name))
        pass

    lst = [binary_file, mie_table, bins]
    del binary_file; del mie_table; del bins
    del lst

def main(args):
    ### define Mie Table based on num bins:
    num_bins = args.nbins #16 # turn into arg
    # logmin and logmax change depending on the mie table used.
    # look at the amplitude for 120 nm and 3000 nm and take log10 of amplitude
    # logmin = 0.871 & logmax = 4.067 are values for mie scripps.
    # "chestnutridge
    # logmin = 1.526; 120 nm
    # logmax = 4.562; 3000nm"

    logmin = args.logmin
    logmax = args.logmax

    # also turn mie conversion table into arg
    # mie_conv_table = pd.read_csv('/Users/FarmerLab/Desktop/PENA/Projects/Mie_Scattering/Conversion_tables/Handix_PSL_MieTable.csv')
    mie_conv_table = pd.read_csv(args.mie)
    table_amp = mie_conv_table['amp_scale']
    table_diam = mie_conv_table['d_nm']
    mie_table = create_bins(num_bins,logmin,logmax,table_amp, table_diam, False)
    mie_table.to_csv(args.directory+"/resulting_mie-table.csv") ### binning check!
    bins = (list(mie_table['Lower Amp'])
    + list(mie_table['Upper Amp'][-1:]))

    ### Walk subdirectories for all binary files ###
    source_directory = args.directory
    # create list of binary files with complete path included
    binary_files_2_process = []
    for dirpath, dirnames, filenames in os.walk(source_directory):
        for filename in [f for f in filenames if f.endswith(".b")]:
            full_path = os.path.join(dirpath, filename)
            binary_files_2_process.append(full_path)
    binary_files_2_process.sort(key=os.path.getmtime)

    print('Extracting POPS binaries...')
    binary_mie_bins = [(f, mie_table, bins) for f in binary_files_2_process] # single arg for pool

    ### Perform extraction is pools of 4 (four processes at once) ###
    if args.multiproc == True:
        print('Multiprocessing enabled... 4 procs being used.')
        p = Pool(4)
        p.map(binary_2_csv,binary_mie_bins)
    else:
        for group in binary_mie_bins:
            try:
                binary_2_csv(group)
            except Exception as e:
                print("File Error!")
                print("{} failed".format(group))
                print("error:\n", e)
    print('done.')

if __name__ == '__main__':
    args = input_arguments()
    main(args)
