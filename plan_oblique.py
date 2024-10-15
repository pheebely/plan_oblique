import os
import math
import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()  # Enable GDAL exceptions

os.environ['GTIFF_SRS_SOURCE'] = 'EPSG'


def planOblique(inDEM, outRaster, angle):
    """
    Applies shearing to a terrain model along the vertical y axis.

    Inputs:
        inDEM -- Input DEM file path.
        outRaster -- Output sheared raster file path.
        angle -- Inclination angle for shearing.
    """
    try:
        # Open the input DEM file
        dem_ds = gdal.Open(inDEM)
        if dem_ds is None:
            raise FileNotFoundError(f"Error opening input DEM file: {inDEM}")

        # Get DEM raster properties
        geotransform = dem_ds.GetGeoTransform()
        projection = dem_ds.GetProjection()
        x_min, pixel_width, _, y_max, _, pixel_height = geotransform
        n_cols = dem_ds.RasterXSize
        n_rows = dem_ds.RasterYSize

        # Read DEM data as a NumPy array
        band = dem_ds.GetRasterBand(1)
        arrIn = band.ReadAsArray().astype(float)

        # Calculate shearing factor based on the inclination angle
        elevationScale = 1.0 / math.tan(math.radians(angle))

        # Get minimum and maximum elevations
        refElevation = np.nanmin(arrIn)
        maxElevation = np.nanmax(arrIn)

        # Compute maximum vertical displacement
        max_dy = (maxElevation - refElevation) * elevationScale

        # Calculate the additional rows needed due to shearing
        dRows = int(max_dy / abs(pixel_height)) + 1
        new_nRows = n_rows + dRows

        # Create enlarged array and initialize with no-data values
        noData = np.nan
        arr = np.full((new_nRows, n_cols), noData)

        # Copy original raster to the enlarged array
        arr[dRows:, :] = arrIn

        # Shear the grid
        arrOut = np.full((new_nRows, n_cols), noData)

        for col in range(n_cols):
            prevRow = new_nRows - 1
            prevZ = np.nan
            for row in reversed(range(prevRow + 1)):
                if np.isnan(arr[row, col]):
                    arrOut[row, col] = np.nan
                else:
                    prevRow = row
                    break

            prevShearedY = y_max - prevRow * abs(pixel_height)
            for row in reversed(range(prevRow + 1)):
                targetY = y_max - row * abs(pixel_height)
                interpolatedZ = np.nan
                for r in reversed(range(prevRow + 1)):
                    z = arr[r, col]
                    if np.isnan(z):
                        move = r - 1
                        while move >= 0 and np.isnan(arr[move, col]):
                            move -= 1
                        prevRow = move
                        interpolatedZ = noData
                        prevZ = noData
                        break

                    shearedY = y_max - r * \
                        abs(pixel_height) + (z - refElevation) * elevationScale
                    if shearedY > targetY:
                        w = (targetY - prevShearedY) / \
                            (shearedY - prevShearedY)
                        interpolatedZ = w * z + (1.0 - w) * prevZ
                        break

                    prevRow = r
                    if shearedY >= prevShearedY:
                        prevShearedY = shearedY
                        prevZ = z

                arrOut[row, col] = interpolatedZ

        # Remove empty rows and adjust the geotransform
        valid_rows = np.where(~np.isnan(arrOut).any(axis=1))[0]
        start = valid_rows[0]
        end = valid_rows[-1] + 1
        new_nRows_trimmed = end - start
        new_y_max = y_max - start * abs(pixel_height)

        # Before the Create call
        print(f"n_cols: {n_cols}, new_nRows_trimmed: {new_nRows_trimmed}")

        # Ensure new_nRows_trimmed is an integer
        new_nRows_trimmed = int(new_nRows_trimmed)

        # Create output raster
        try:
            driver = gdal.GetDriverByName('GTiff')
            out_ds = driver.Create(
                outRaster, n_cols, new_nRows_trimmed, 1, gdal.GDT_Float32)
            if out_ds is None:
                raise Exception(f"Failed to create output raster: {outRaster}")
        except Exception as e:
            print(f"Error while creating output raster: {e}")

        out_ds.SetGeoTransform(
            (x_min, pixel_width, 0, new_y_max, 0, pixel_height))
        out_ds.SetProjection(projection)

        # Write the sheared array to the output raster
        out_band = out_ds.GetRasterBand(1)
        out_band.WriteArray(arrOut[start:end, :])
        out_band.SetNoDataValue(noData)

        # Close datasets
        dem_ds = None
        out_ds = None
        print("Sheared raster successfully saved.")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == '__main__':
    inDEM = 'inDEM'
    outRaster = 'outRaster'
    angle = 25  # Example angle
    planOblique(inDEM, outRaster, angle)
