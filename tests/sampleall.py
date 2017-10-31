from downsample import downsample
from glob import glob
import time

infilenames = glob('../VS0169/*.nc')
infilenames.extend(glob('../VS0169/*.volume'))

outfilenames = ['sampled/sampled_133_' + line.split('/')[-1] for line in infilenames]

for infn, outfn in zip(infilenames, outfilenames):
    start = time.time()
    fnshort = infn.split('/')[-1]
    print('\nNow processing:', fnshort)


    factors = (1, 3, 3)
    downsample(infn, outfn, factors)
    end = time.time()
    print('Finished processing: {}!\n {} seconds elapsed total \n'.format(fnshort, end-start))
