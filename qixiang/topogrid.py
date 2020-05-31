#!/usr/bin/python
# -*- coding:utf-8 -*-
#
# Neng Lu
# nengl@student.unimelb.edu.au
# ANU & Unimelb
# Canberra, Australia
# 
# Reference packages: 
# grid.py
# topopy: https://github.com/geolovic/topopy
# 
# Version: 1.0
# First version 24 May, 2020
# Last modified 26 May, 2020

import gdal
import numpy as np
from scipy import ndimage
from skimage.morphology import reconstruction

try: 
    import matplotlib.pyplot as plt
    from matplotlib import cm
    PLT = True
except Exception as e:
    PLT = False

class TopoGrid(object):
    """
    Class for holding and manipulating gridded data.
    ============
    Attributes
    
    _size : *Tuple*, size of the grid (XSize, YSize)
    _dims : *Tuple*, size of the internal array (nrow, ncol)
    _geot : *Tuple*, GeoTranstorm matrix of the grid *(ULx, Cx, Tx, ULy, Ty, Cy)*
            * ULx = Upper-Left X coordinate (upper-left corner of the pixel)
            * ULy = Upper-Left Y coordinate (upper-left corner of the pixel)
            * Cx = X Cellsize
            * Cy = Y Cellsize (negative value)
            * Tx = Rotation in X axis
            * Ty = Rotation in Y axis
    _cellsize : *Tuple*, grid cellsize [dx,dy]
    _proj   : *String*, the projection of the grid in WKT
    _ncells : *Int*, total number of cells of the Grid
    _nodata : *np.ndarray*, the nodata of the grid
    _array  : *np.ndarray*, the data of the grid
    _tipo   : *Str*, the dtype of _array
    ============
    Methods
    
    File I/O
    read_raster : Read a raster file and add the data
    save : Save the grid data array to tiff file
    plot : Imshow the grid data array 
    ---------
    Hydrologic
    identify_flats : Return two TopoGrid() (or numpy.ndarrays) with flats and sills 
    fill_sinks  : Fill sinks using scikit-image reconstruction algorithm
    fill_sinks2 : Fill sinks using fill depressions/sinks in floating point array
    ---------
    Data Processing
    get_***    : Get the attributes
    set_***    : Set the attributes 
    is_inside  : Check if one point is inside the raster
    copy_layout: Copy all the parameters from TopoGrid instance except grid data (and nodata).
    xy_2_cell  : Get row col indexes from XY coordinates
    cell_2_xy  : Get XY coordinates from (row, col) cell indexes
    ind_2_cell : Get row col indexes from cells linear indexes (row-major, C-style)
    cell_2_ind : Get cell linear indexes from row and column indexes
    read_array : Reads the internal array of the Grid instace
    find, max, min, mean : return the regading results of the grid data array
    ===========   
    References
    This class is reorganized the Class PRaster(), Grid(), DEM() in grid.py 
    from repo [geolovic/topopy](https://github.com/geolovic/topopy)
    
    """
    def __init__(self, path=""):
        """    
        ============
        Parameters 
        
        path : *Str*, the path of the raster file
        """      
        if path:
            self.read_raster(path)
               
        else:
            self._size = (1, 1)
            self._dims = (1, 1)
            self._geot = (0., 1., 0., 0., 0., -1.)
            self._cellsize = (self._geot[1], self._geot[5])
            self._proj = ""
            self._ncells = 1
            self._nodata = None
            self._array = np.array([[0]], dtype=np.float)
            self._tipo = str(self._array.dtype)
            
    def read_raster(self,path=""):
        if path:
            raster = gdal.Open(path)
            if not raster:
                raise FileNotFoundError
            banda = raster.GetRasterBand(1)
            self._size = (banda.XSize, banda.YSize)
            self._dims = (banda.YSize, banda.XSize)
            self._geot = raster.GetGeoTransform()
            self._cellsize = (self._geot[1], self._geot[5])
            self._proj = raster.GetProjection()
            self._ncells = banda.XSize * banda.YSize
            self._nodata = banda.GetNoDataValue()
            self._array = banda.ReadAsArray()
            self._tipo = str(self._array.dtype)
        else:
            raise FileNotFoundError
    
    def get_size(self):
        return self._size
    
    def get_dims(self):
        return self._dims
    
    def get_ncells(self):
        return self._ncells
    
    def get_projection(self):
        return self._proj
    
    def get_cellsize(self):
        return self._cellsize
    
    def get_geotransform(self):
        return self._geot
    
    def is_inside(self, x, y):
        """
        Check if one point is inside the raster
        """
        row, col = self.xy_2_cell(x, y)
        rowinside = np.logical_and(row >= 0, row < self._dims[0])
        colinside = np.logical_and(col >= 0, col < self._dims[1])
        inside = np.logical_and(rowinside, colinside)
        return np.all(inside)
    
    def get_extent(self):
        """
        Return a tuple (XMin, XMax, YMin, YMax) with the extension of the Grid
        """
        xmin = self._geot[0]
        xmax = self._geot[0] + self._size[0] * self._cellsize[0]
        ymin = self._geot[3] + self._size[1] * self._cellsize[1]
        ymax = self._geot[3]
        
        return (xmin, xmax, ymin, ymax)
    
    def copy_layout(self, grid):
        """
        Copy all the parameters from another TopoGrid instance except grid data (and nodata)
        """
        self._size = grid.get_size()
        self._dims = grid.get_dims()
        self._geot = grid.get_geotransform()
        self._cellsize = grid.get_cellsize()
        self._proj = grid.get_projection()
        self._ncells = grid.get_ncells()

    def xy_2_cell(self, x, y):
        """
        Get row col indexes from XY coordinates
        =========== 
        Parameters

        x : number, list, or numpy.ndarray, X coordinates (number, list, or numpy.ndarray)
        y : Y coordinates (number, list, or numpy.ndarray)
        ===========     
        Return

        *Tuple* : Tuple with (row, col) indices as np.ndarrays
        """
        x = np.array(x)
        y = np.array(y)       
        row = (self._geot[3] - y) / (self._geot[5] * -1)
        col = (x - self._geot[0]) / self._geot[1]
        return row.astype(np.int32), col.astype(np.int32)

    def cell_2_xy(self, row, col):
        """
        Get XY coordinates from (row, col) cell indexes
        =========== 
        Parameters
        
        row : row indexes (number, list, or numpy.ndarray)
        col : column indexes (number, list, or numpy.ndarray)
        ===========     
        Return
        
        **tuple** : Tuple with (x, y) coordinates as np.ndarrays
        """
        row = np.array(row)
        col = np.array(col)
        x = self._geot[0] + self._geot[1] * col + self._geot[1] / 2
        y = self._geot[3] + self._geot[5] * row + self._geot[5] / 2
        return x, y
    
    def ind_2_cell(self, ind):
        """
        Get row col indexes from cells linear indexes (row-major, C-style)
        =========== 
        Parameters
        
        ind : linear indexes (number, list, or numpy.ndarray)
        =========== 
        Return
        
        **tuple** : Tuple with (row, col) indices as numpy.ndarrays
        """
        return np.unravel_index(ind, self._dims) 
    
    def cell_2_ind(self, row, col):
        """
        Get cell linear indexes from row and column indexes
        =========== 
        Parameters
        
        row : row indexes (number, list, or numpy.ndarray)
        col : column indexes (number, list, or numpy.ndarray)
        ===========     
        Return
        
        **numpy.array** : Array with linear indexes (row-major, C-style)
        """
        return np.ravel_multi_index((row, col), self._dims)
    
    def set_array(self, array):
        """
        Set the data array for the current Grid object. 
        If the current Grid is an empty Grid [get_size( ) = (1, 1)], any input array is valid
        If The current Grid is not an empty Grid, the input array should match Grid dimensions
        =========== 
        Parameters
        
        array : numpy.ndarray, Numpy array with the data
        """
        # If the Grid is an empty Grid, any array is valid
        if self._size == (1, 1):       
            self._dims = array.shape
            self._size = (array.shape[1], array.shape[0])    
            self._array = np.copy(array)
            self._tipo = str(self._array.dtype)
        # If the Grid is not an empty Grid, input array shape must coincide with internal array
        elif array.shape == self._dims:    
            self._array = np.copy(array)
            self._tipo = str(self._array.dtype)
        else:
            return 0
        
    def read_array(self, ascopy=False):
        """
        Reads the internal array of the Grid instace
        =========== 
        Parameters
        
        ascopy : bool
          If True, the returned array is a memory view of the Grid original array.
        =========== 
        Return
        
        **numpy.ndarray** : Internal array of the current Grid object
        """
        if ascopy:
            return np.copy(self._array)
        else:
            return self._array
    
    def find(self):
        """
        Find the non-zero elements in the array. Return a tuple of arrays with
        row and col positions.
        """
        return np.where(self._array > 0)
    
    def max(self):
        """
        Return the maximun value of the Grid
        """
        datapos = np.where(self._array != self._nodata)
        return np.max(self._array[datapos])
    
    def min(self):
        """
        Return the minimun value of the Grid
        """
        datapos = np.where(self._array != self._nodata)
        return np.min(self._array[datapos])
    
    def mean(self):
        """
        Return the mean value of the Grid
        """
        datapos = np.where(self._array != self._nodata)
        return np.mean(self._array[datapos])
    
    def set_value(self, row, col, value):
        """
        Set the value for a cell of the grid at (row, col)
        =========== 
        Parameters
        
        row, col : int 
          Row and column indexes
        value : number
          Value for the cell (row, col)
        """
        self._array[row, col] = value
    
    def get_value(self, row, col):
        """
        Get the value for a cell/s of the grid at (row, col)
        =========== 
        Parameters
        
        row, col : ints or numpy.ndarrays
          Row and column indexes
        """
        return self._array[row, col]
    
    def get_nodata(self):
        """
        Return the value for NoData in the grid. This value could be None if NoData
        are not defined in the grid.
        """
        return self._nodata
    
    def set_nodata(self, value):
        """
        Sets the nodata value for the Grid
        """
        self._nodata = value
    
    def get_nodata_pos(self):
        """
        Return the position of the NoData values as a tuple of two arrays (rows, columns)
        """
        if self._nodata is None:
            return (np.array([], dtype=np.int), np.array([], dtype=np.int))
        else:
            return np.where(self._array == self._nodata)
        
    def nan_2_nodata(self):
        """
        Change nan values to NoData (if Grid nodata is defined). 
        """
        if self._nodata is None:
            return
        idx = np.isnan(self._array)
        self._array[idx] = self._nodata
    
    def values_2_nodata(self, value):
        """
        Change specific values to NoData (if Grid nodata is defined). 
        =========== 
        Parameters
        
        value : Value or sequence of values that will be changed to NoData
        """
        if self._nodata is None:
            return
        
        if type(value) == int or type(value)==float:
            ind = np.where(self._array==value)
            self._array[ind] = self._nodata
        else:
            for val in value:
                ind = np.where(self._array == val)
                self._array[ind] = self._nodata        
    
    def plot(self, ax=None):
        """
        Plots the grid in a new Axes or in a existing one
        =========== 
        Parameters
        
        ax : matplotlib.Axe
          If is not defined, the function will use plt.imshow()
        """
        if not PLT:
            return
        
        arr = self.read_array()
        if self._nodata:
            mask = self._array == self._nodata
            arr = np.ma.array (self._array, mask = mask)
        
        extent = self.get_extent()
        if ax:
            ax.imshow(arr,extent=extent,cmap=cm.rainbow)
        else:
            plt.imshow(arr,extent=extent,cmap=cm.rainbow)
    
    def save(self, path):
        """
        Saves the grid in the disk
        =========== 
        Parameters
        
        path : str 
          Path where new raster will be saved
        """
        # Check if the type of the internal array is compatible with gdal
        if str(self._tipo) not in NTYPES.keys():
            return 0
        else:
            tipo = NTYPES[str(self._tipo)]
        
        # Prepare driver to write raster
        driver = gdal.GetDriverByName("GTiff")
        raster = driver.Create(path, self._size[0], self._size[1], 1, tipo)
        if not raster:
            return
        raster.SetGeoTransform(self._geot)
        raster.SetProjection(self._proj)
        
        if not self._nodata is None:
            raster.GetRasterBand(1).SetNoDataValue(self._nodata)
        raster.GetRasterBand(1).WriteArray(self._array)
        
    def identify_flats(self, nodata=True, as_array=False):
        """
        This functions returns two topopy.Grids (or numpy.ndarrays) with flats and sills. 
        Flats are defined  as cells without downward neighboring cells. Sills are cells where 
        flat regions spill over into lower terrain. It the DEM has nodata, they will be maintained
        in output grids.
        =========== 
        Parameters
        
        nodata : *Boolean*, indicates if notada are kept (True) or replaced by zeros (False)
        as_array : *Boolean*, indicates if outputs are boolean numpy.ndarrays (True) or a Grid objects (False)
        =========== 
        Returns
        
        *Tuple* : Tuple (flats, sills) with two topopy.Grid (as_array=False) or numpy.ndarray (as_array=True)
        =========== 
        References
        
        This algoritm is adapted from identifyflats.m by Wolfgang Schwanghart 
        (version of 17. August, 2017) included in TopoToolbox matlab codes.
        
        Schwanghart, W., Kuhn, N.J., 2010. TopoToolbox: A set of Matlab functions 
        for topographic analysis. Environ. Model. Softw. 25, 770–781. 
        https://doi.org/10.1016/j.envsoft.2009.12.002
        
        Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
        MATLAB-based software for topographic analysis and modeling in Earth 
        surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
        """

        z_arr = np.copy(self._array)
        
        # Change nodata to very low values
        nodata_ids = self.get_nodata_pos()
        z_arr[nodata_ids] = -9999
        
        footprint = np.ones((3, 3), dtype=np.int8)
        # Identify flats throught a image binary erosion
        # Flats will be True where cells don't have lower neighbors
        flats = ndimage.morphology.grey_erosion(z_arr, footprint=footprint) == z_arr
        
        # Remove flats from the borders
        flats[:,0] = False
        flats[:,-1] = False
        flats[0,:] = False
        flats[-1,:] = False
        
        # Remove flats for nodata values and cells bordering them
        flats[nodata_ids] = False
        auxmat = np.zeros(flats.shape, dtype="bool")
        auxmat[nodata_ids] = True
        nodata_bord = ndimage.morphology.grey_dilation(auxmat, footprint=footprint)
        flats[nodata_bord] = False
        
        # Identify sills
        sills = np.empty(z_arr.shape)
        sills.fill(-9999.)
        sills[flats] = z_arr[flats]
        aux_dil = ndimage.morphology.grey_dilation(sills, footprint=footprint)
        sills = np.logical_and(aux_dil == z_arr, np.logical_not(flats))
        sills[nodata_ids] = False
        
        # Prepare outputs
        if as_array:
            return flats, sills
        else:
            res = []
            for arr in [flats, sills]:
                grid = TopoGrid()
                grid.copy_layout(self)
                grid.set_nodata(-1)
                grid._tipo = 'int8'
                grid.set_array(arr.astype(np.int8))
                
                if nodata:
                    grid._array[nodata_ids] =  -1
                else:
                    grid._array[nodata_ids] =  0
                    grid.set_nodata(None)
            
                res.append(grid)
    
            return res

    def fill_sinks(self, as_array=False):
        """
        Fill sinks in a DEM using scikit-image reconstruction algorithm
        ===========         
        Returns
        
        *TopoGrid* or *np.ndarray* : Filled DEM
        """
        # Get the seed to start the fill process
        seed = np.copy(self._array)
        seed[1:-1, 1:-1] = self._array.max()
        
        # Here we modified the inner array to save memory (instead of working with a copy)
        nodata_pos = self.get_nodata_pos()
        self._array[nodata_pos] = -9999.
        
        filled = reconstruction(seed, self._array, method='erosion')
        filled = filled.astype(self._array.dtype)
        # Put back nodata values and change type
        if self._nodata:
            self._array[nodata_pos] = self._nodata
            filled[nodata_pos] = self._nodata
        
        if as_array:
            # Return filled DEM as numpy.ndarray
            return filled
        else:
            # Return filled DEM as topopy.DEM
            filled_dem = TopoGrid()
            filled_dem.copy_layout(self)
            filled_dem.set_array(filled)
            filled_dem.set_nodata(self._nodata)
            return filled_dem
    
    def fill_sinks2(self, four_way=False):
        """
        Fill sinks method adapted from  fill depressions/sinks in floating point array
        =========== 
        Parameters
        
        input_array : *np.ndarray*, input array to be filled
        four_way : *Boolean*, searchs the 4 (True) or 8 (False) adjacent cells
        ===========  
        Returns
        
        *TopoGrid* or *np.ndarray* : Filled DEM
        ===========
        References
       
        This algorithm has been adapted (with minor modifications) from the 
        Charles Morton slow fill algorithm (with ndimage and python 3 was not slow
        at all). 
        
        Cite:
        Soile, P., Vogt, J., and Colombo, R., 2003. Carving and Adaptive
        Drainage Enforcement of Grid Digital Elevation Models.
        Water Resources Research, 39(12), 1366
        
        Soille, P., 1999. Morphological Image Analysis: Principles and
        Applications, Springer-Verlag, pp. 173-174
    
        """
        # Change nan values to a very low value
        copyarr = np.copy(self._array)
        nodata_pos = self.get_nodata_pos()
        copyarr[nodata_pos] = -9999.
        
        # Set h_max to a value larger than the array maximum to ensure
        #   that the while loop will terminate
        h_max = copyarr.max() + 100
    
        # Build mask of cells with data not on the edge of the image
        # Use 3x3 square Structuring element
        inside_mask = ndimage.morphology.binary_erosion(
            np.isfinite(copyarr),
            structure=np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]).astype(np.bool))
    
        # Initialize output array as max value test_array except edges
        output_array = np.copy(copyarr)
        output_array[inside_mask] = h_max
    
        # Array for storing previous iteration
        output_old_array = np.copy(copyarr)
        output_old_array[:] = 0
    
        # Cross structuring element
        if four_way:
            el = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]).astype(np.bool)
        else:
            el = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]).astype(np.bool)
    
        # Iterate until marker array doesn't change
        while not np.array_equal(output_old_array, output_array):
            output_old_array = np.copy(output_array)
            output_array = np.maximum(
                copyarr,
                ndimage.grey_erosion(output_array, size=(3, 3), footprint=el))

        # Put back nodata values and change type
        if self._nodata:
            output_array[nodata_pos] = self._nodata
        # Create output filled DEM
        filled_dem = TopoGrid()
        filled_dem.copy_layout(self)
        filled_dem.set_array(output_array)
        filled_dem.set_nodata(self._nodata)
        return filled_dem
