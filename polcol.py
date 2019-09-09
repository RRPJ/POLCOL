from rht import xyt_name_factory, main
from POLCOL_tools import cutout, plot_dist_polangle, make_polarplot


# Variables for cutting out a square from the original fits file
x1, x2 = (500, 700)  # x pixel borders
y1, y2 = (500, 700)  # y pixel borders

position = ((x1 + x2) / 2., (y1 + y2) / 2.)  # x,y
size = (y2 - y1, x2 - x1)  # NOTE: y,x
image_filename = '/dir/file.fits'
catalogue = '/dir/file.csv'
# square_size = 7  # Must be odd, not implemented

# Calculation parameters for RHT (Standard 55, 0.70, 15, True)
wlen = 55
frac = 0.70
smr = 15
original = True

# Run the cutout script and save its resulting filename
# Cut a piece or whole image? True: cut a piece, false: whole image
cut = True
if(cut):
    cutout_filename = cutout(image_filename, position, size)
else:
    cutout_filename = image_filename  # For whole picture

# Run the RHT
main(cutout_filename, wlen=wlen, frac=frac, smr=smr, drht=not original)

# Get the RHT output filename
rht_output_filename = xyt_name_factory(
    cutout_filename, wlen, smr, frac, original)

# Plotting:
polar = False  # Polarplot? True if desired
distance = True  # Distance - polarisation angle plot? True if desired

if(polar):
    make_polarplot(rht_output_filename)

if(distance):
    plot_dist_polangle(catalogue, rht_output_filename)
