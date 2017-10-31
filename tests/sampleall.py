from downsample import downsample
from glob import glob
import time

'''
This script reads in all netCDF data from the VS0169_Full directory and downsamples it, 
saving the results to the VS0169_Sampled directory. 
'''

infilenames = glob('../Data/VS0169_Full/*.nc')
infilenames.extend(glob('../Data/VS0169_Full/*.volume'))

outfilenames = ['../Data/VS0169_Sampled/sampled_133_' + line.split('/')[-1] for line in infilenames]

for infn, outfn in zip(infilenames, outfilenames):
    start = time.time()
    fnshort = infn.split('/')[-1]
    print('\nNow processing:', fnshort)


    factors = (1, 3, 3)
    downsample(infn, outfn, factors)
    end = time.time()
    print('Finished processing: {}!\n {} seconds elapsed total \n'.format(fnshort, end-start))
