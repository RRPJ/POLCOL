# POLCOL
A POLarisation COrreLator for astronomers working with radio data.

This code works with the Rolling Hough Transform by Susan Clarke - https://github.com/seclark/RHT
Be sure to download this repository and add the rht.py and RHT_tools.py files to the same directory as this repository.

The current implementation has stellar catalogues ONLY. The script assumes the catalogue to have 'ra' and 'dec' axes. The error in distance is assumed to be lower and upper estimates. In other words: dist - dist_err_low is the lower distance, while dist_err_high - dist is the upper distance. Be careful when implementing other catalogues.

To use, simply open polcol.py in a text editor. There are a couple of things to change before running:

1. image_filename: input your image filename in string form
2. catalogue: input your stellar catalogue in string form. Future releases may also support other tracers.
3. Decide if you want to cut out a rectangle or take the whole image to run into the RHT. For the first, set "cut" to True, and change x1, x2 and y1, y2 to values desired. The script cuts out a rectangle of (x2 - x1) by (y2 - x1). 
4. Set the RHT parameters to your choice. Default values are given in comment. 
5. Open rht.py in a text editor, and set OUTPUT to a directory of your choice. This directory doesn't need to exist yet.
6. Set polar to True if you want a polarplot, and distance to True if you want a distance - polarisation angle plot. The script saves the images in current working directory.
7. Go to the function put_XYT and edit the variable header_keywords so that you have all the keywords in there that you want to keep. This is useful if you want to preserve certain axes you had in your input image, like RA and DEC.

Once that's done, save the files and simply run polcol.py by typing 'python polcol.py' in your terminal. 
