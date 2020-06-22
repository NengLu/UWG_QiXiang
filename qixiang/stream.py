#!/usr/bin/python
# -*- coding:utf-8 -*-
#
# Neng Lu
# nengl@student.unimelb.edu.au
# ANU & Unimelb
# Canberra, Australia
# 
# Reference packages: 
# network.py
# topopy: https://github.com/geolovic/topopy
# 
# Version: 1.0
# First version 28 May, 2020
# Last modified 31 May, 2020

import numpy as np
from scipy.sparse import csc_matrix
from . import TopoGrid, Flow

class Stream(Flow):
    """
    Class that define a stream object, defined by applying a threshold to 
    a flow accumulation raster derived from the *Flow* (Flow.get_flow_accumulation()).
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
    _proj   : *Str*, the projection of the grid in WKT
    _ncells : *Int*, total number of cells of the Grid
    _nodata_pos : *Numpy.ndarray*, the nodata postion of the grid
    _ix     : *Numpy.ndarray*, the id of the givers
    _ixc    : *Numpy.ndarray*, the id of the receiver
    
    _threshold : 
    _ax : *Numpy.ndarray*, area size of the channel cells 
    _zx : *Numpy.ndarray*, elevation of the channel cells 
    _dd : *Numpy.ndarray*, giver-receiver distance of the channel cells 
    _dx : *Numpy.ndarray*, distance to mouth of the channel cells 
    _x  : *Numpy.ndarray*, x-coordinate vector of the channel cells 
    _y  : *Numpy.ndarray*, y-coordinate vector of the channel cells 
    
    _chi : *Numpy.ndarray*, chi of the channel cells 
    _slp : *Numpy.ndarray*, slope of the channel cells 
    _ksn : *Numpy.ndarray*, ksn of the channel cells  
    ============
    Methods
    
    File I/O
    save : Save the Network instance to disk
    ---------
    Hydrologic
    calculate_chi : Calculate chi (_chi) for channel cells
    calculate_gradients : Calculate gradients (slope or ksn, _slp or _ksn) for all channel cells
    
    get_stream_poi : Find points of interest of the drainage network
    get_streams : Output a grid representation of the Network object
    get_stream_segments : Extract a drainage network by using a determined area threshold
                          and output the numerated stream segments.

    get_klargestconncomps : Return k largest connected components in an instance of *Stream*                 
    get_largerpart_lssorder : Get the part with strahler stream order (i-1) in the trunk
    get_trunk : Extract trunk stream of the connected component by 
                 Strahler stream order and the size of stream nodes number       
    =========== 
    References
    
    This class is reorganized the Class Network() in network.py 
    from repo [geolovic/topopy](https://github.com/geolovic/topopy) by José Vicente Pérez Peña,
    and added some functions which function like *klargestconncomps.m* and *trunk.m* 
    from TopoToolbox matlab codes by Wolfgang Schwanghart.
                   
    Cite:            
    Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
    MATLAB-based software for topographic analysis and modeling in Earth 
    surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
    """

    def __init__(self, dem=None, flow=None, threshold=0, verbose=False, thetaref=0.45, npoints=5):    
        """    
        ============
        Parameters   
        
        dem  : *TopoGrid*
        flow : *Flow* 
        threshold : *Int*, number of the cells to initiate a channel
        thetaref  : *Float*, m/n coeficient to calculate chi values in each channel cell 
        verbose   : *Boolean*, indicate if calculate chi, ksn, and slp
        npoints   : *Int*, number of points to calculate slope and ksn in each cell. 
                     Slope and ksn are calculated with a moving window of (npoints * 2 + 1) cells.
        """

        # Set raster properties
        self._size = flow.get_size()
        self._dims = flow.get_dims()
        self._geot = flow.get_geotransform()
        self._cellsize = flow.get_cellsize()
        self._proj = flow.get_projection()
        self._ncells = flow.get_ncells()
              
        # Get a threshold if not specified (Default 0.5% of the total number of cells)
        if threshold == 0:
            threshold = self._ncells * 0.005
        self._threshold = threshold        
        
        # Get sort Nodes for channel cells
        # In TopoToolbox:
        # FD = FLOWobj(DEMf);
        # A = flowacc(FD);
        # W = A>th_discharge;
        # S = STREAMobj(FD,W);
        
        fac = flow.get_flow_accumulation(nodata=False, asgrid=False)
        w = fac > threshold
        w = w.ravel()
        I   = w[flow._ix]
        self._ix  = flow._ix[I]
        self._ixc = flow._ixc[I]        
        
        # Get area for channel cells
        self._ax = fac.ravel()[self._ix] * self._cellsize[0] * self._cellsize[1] * -1 # Area in map units
        
        # Get elevation for channel cells
        self._zx = dem.read_array().ravel()[self._ix]
        
        # Get distances to mouth (self._dx) and giver-receiver distances (self._dd)
        di = np.zeros(self._ncells)
        self._dd = np.zeros(self._ix.shape) # Giver-Receiver distance
        for n in np.arange(self._ix.size)[::-1]:
            grow, gcol = self.ind_2_cell(self._ix[n])
            rrow, rcol = self.ind_2_cell(self._ixc[n])
            gx, gy = self.cell_2_xy(grow, gcol)
            rx, ry = self.cell_2_xy(rrow, rcol)
            d_gr = np.sqrt((gx - rx)**2 + (gy - ry)**2)
            self._dd[n] = d_gr
            di[self._ix[n]] = di[self._ixc[n]] + d_gr
        self._dx = di[self._ix]
                   
        if verbose == True:
            # Get chi values using the input thetaref
            self._thetaref = thetaref
            self.calculate_chi(thetaref)
        
            # Calculate slopes and ksn
            self.calculate_gradients(npoints, 'slp')
            self.calculate_gradients(npoints, 'ksn')
            
    def save(self, path):
        """
        Saves the Network instance to disk. It will be saved as a numpy array in text format.
        The first three rows will have the information of the raster
        ===========
        Parameters
        
        path : *Str*, path to save the network object, with *.net extension
        """
    
        path = os.path.splitext(path)[0]
            
        # Create *.net file with properties
        netfile = open(path + ".net", "w")
        params = [self._size, self._cellsize, self._ncells, self._ksn_npoints, self._slp_npoints, self._thetaref] 
        linea = ";".join([str(param) for param in params]) + "\n"
        netfile.write(linea)
        linea = ";".join([str(elem) for elem in self._geot]) + "\n"
        netfile.write(linea)
        netfile.write(str(self._proj))
        netfile.close()
        
        # Create data array
        data_arr = np.array((self._ix, self._ixc, self._ax, self._dx, self._zx,
                             self._chi, self._slp, self._ksn, self._r2slp, 
                             self._r2ksn, self._dd)).T
        
        # Save the network instance as numpy.ndarray in text format
        np.save(path + ".npy", data_arr)
            
    def calculate_chi(self, thetaref=0.45, a0=1.0):
        """
        Function that calculate chi for channel cells
        ============
        Parameters   
        
        thetaref : *Float*, m/n coeficient to calculate chi
        a0 : *Float*, reference area to avoid dimensionality (usually don't need to be changed)
        """
        chi = np.zeros(self._ncells)
        for n in np.arange(self._ix.size)[::-1]:
            chi[self._ix[n]] = chi[self._ixc[n]] + (a0 * self._dd[n]/self._ax[n]**thetaref)            
        self._chi = chi[self._ix]
        self._thetaref = thetaref
      
    def calculate_gradients(self, npoints, kind='slp'):
        """
        This function calculates gradients (slope or ksn) for all channel cells. 
        Gradients of each cell are calculated by linear regression using a number
        of points (npoints) up and downstream.
        ============
        Parameters 
        
        npoints : *Int*, window to analyze slopes. Slopes are calculated by linear regression using a window
                  of npoints * 2 + 1 pixel (using the central pixel)
        kind : *Str*, {'slp', 'ksn'}
        """
        if kind not in ['slp', 'ksn']:
            kind = 'slp'
        
        # Get arrays depending on type
        if kind == 'slp':
            x_arr = self._dx
            y_arr = self._zx
        elif kind == 'ksn':
            x_arr = self._chi
            y_arr = self._zx
            
        # Get ixcix auxiliar array
        ixcix = np.zeros(self._ncells, np.int)
        ixcix[self._ix] = np.arange(self._ix.size)
        
        # Get heads and sorted them by elevation
        heads = self.get_stream_poi("heads", "IND")
        elev = self._zx[ixcix[heads]]
        spos = np.argsort(-elev)
        heads = heads[spos]
        winlen = npoints * 2 + 1
        gi = np.zeros(self._ncells)
        r2 = np.zeros(self._ncells)
        
        # Taking sequentally all the heads and compute downstream flow
        for head in heads:
            processing = True
            lcell = head
            mcell = self._ixc[ixcix[head]]
            fcell = self._ixc[ixcix[mcell]]
        
            if ixcix[fcell] == 0 or ixcix[self._ixc[ixcix[fcell]]] == 0:
                continue
            
            # Obtenemos datos de elevacion y distancias
            win = [fcell, mcell, lcell]
            xi = x_arr[ixcix[win]]
            yi = y_arr[ixcix[win]]
            
            # Calculamos pendiente de celda central por regresion
            poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
            # To avoid issues with horizontal colinear points
            if yi.size * yi.var() == 0:
                R2 = 1
            else:
                R2 = float(1 - SCR/(yi.size * yi.var()))
            
            g = poli[0]
            gi[mcell] = g
            r2[mcell] = R2
            
            while processing:
                # Cogemos la siguiente celda del canal (next_cell)
                fcell = win[0]
                next_cell = self._ixc[ixcix[fcell]]
                
                # Comprobamos la siguiente celda del canal
                # Si la siguiente celda del canal no es cero
                if ixcix[next_cell] != 0:
                    # Añadimos siguiente celda
                    win.insert(0, next_cell)
                    fcell = win[0]
                    # Movemos celda central
                    mcell = self._ixc[ixcix[mcell]]
        
                    if len(win) < winlen:
                        # Si estamos al principio del canal, win crece
                        next_cell = self._ixc[ixcix[fcell]]
                        win.insert(0, next_cell)
                        fcell = win[0]
                    else:
                        # Si no estamos al principio, eliminamos celda final
                        win.pop()
                # Si la siguiente celda es cero, no añadimos celdas, sino que win decrece
                else:
                    mcell = self._ixc[ixcix[mcell]]
                    win.pop()
                    win.pop()
                    if len(win) == 3:
                        processing = False
                        gi[fcell] = 0.00001
                        r2[fcell] = 0.00001
                        
                # Obtenemos datos de elevacion y distancias
                xi = x_arr[ixcix[win]]
                yi = y_arr[ixcix[win]]
                
                # Calculamos pendiente de celda central por regresion
                poli, SCR = np.polyfit(xi, yi, deg = 1, full = True)[:2]
                
                # To avoid issues with horizontal colinear points
                if yi.size * yi.var() == 0:
                    R2 = 1
                else:
                    R2 = float(1 - SCR/(yi.size * yi.var()))
            
                g = poli[0]
                    
                if gi[mcell] == 0:
                    gi[mcell] = g
                    r2[mcell] = R2
                else:
                    processing = False
        
        if kind == 'slp':
            self._slp = gi[self._ix]
            self._r2slp = r2[self._ix]
            self._slp_npoints = npoints
        elif kind == 'ksn':
            self._ksn = gi[self._ix]
            self._r2ksn = r2[self._ix]
            self._ksn_npoints = npoints   

    def get_stream_poi(self, kind="heads", coords="CELL"):
        """
        This function finds points of interest of the drainage network. 
        These points of interest can be 'heads', 'confluences' or 'outlets'.
        =========== 
        Parameters
        
        kind : *Str*, {'heads', 'confluences', 'outlets'}, Kind of point of interest to return.
        coords : *Str*, {'CELL', 'XY', 'IND'}, Output coordinates for the stream point of interest. 
        ===========   
        Returns
        
        *Numpy.ndarray*, array with one (id) or two columns ([row, col] or [xi, yi] - depending on coords) 
        with the location of the points of interest 
        ===========   
        References
        
        The algoritms to extract the point of interest have been adapted to Python 
        from Topotoolbox matlab codes developed by Wolfgang Schwanghart (version of 17. 
        August, 2017). These smart algoritms use sparse arrays with giver-receiver indexes, to 
        derive point of interest in a really efficient way. 
        
        Cite:     
        Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
        MATLAB-based software for topographic analysis and modeling in Earth 
        surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
        """
        # Check input parameters
        if kind not in ['heads', 'confluences', 'outlets']:
            kind = 'heads'
        if coords not in ['CELL', 'XY', 'IND']:
            coords = 'CELL'

        # Get grid channel cells
        w = np.zeros(self._ncells, dtype=np.bool)
        w[self._ix] = True
        w[self._ixc] = True
            
        # Build a sparse array with giver-receivers cells
        aux_vals = np.ones(self._ix.shape, dtype=np.int8)
        sp_arr = csc_matrix((aux_vals, (self._ix, self._ixc)), shape=(self._ncells, self._ncells))
            
        # Get stream POI according the selected type
        if kind == 'heads':
            # Heads will be channel cells marked only as givers (ix) but not as receivers (ixc) 
            sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
            out_pos = (sum_arr == 0) & w
        elif kind == 'confluences':
            # Confluences will be channel cells with two or givers
            sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
            out_pos = sum_arr > 1
        elif kind == 'outlets':
            # Outlets will be channel cells marked only as receivers (ixc) but not as givers (ix) 
            sum_arr = np.asarray(np.sum(sp_arr, 1)).ravel()
            out_pos = np.logical_and((sum_arr == 0), w)  
                
        out_pos = out_pos.reshape(self._dims)
        row, col = np.where(out_pos)
            
        if coords=="CELL":
            return np.array((row, col)).T
        elif coords=="XY":
            xi, yi = self.cell_2_xy(row, col)
            return np.array((xi, yi)).T
        elif coords=="IND":
            return self.cell_2_ind(row, col)
        
    def get_streams(self, asgrid=True):
        """
        This function outputs a grid representation of the Network object
        ===========
        Parameters
        
        asgrid : *Boolean*, indicates if the network is returned as *TopoGrid* (True) or as a *numpy.array* (False)
        """
        # Get grid channel cells
        w = np.zeros(self._ncells, dtype=np.int8)
        w[self._ix] = 1
        w[self._ixc] = 1
        w = w.reshape(self._dims)
        # Return grid
        if asgrid:
            return self._create_output_grid(w, 0)
        else:
            return w

    def get_stream_segments(self, asgrid=True):
        """
        This function extract a drainage network by using a determined area threshold
        and output the numerated stream segments.
        ===========
        Parameters
        
        asgrid : *Boolean*, indicates if the network is returned as *TopoGrid* (True) or as a *numpy.array* (False)
        """
        # Get heads and confluences and merge them
        head_ind = self.get_stream_poi("heads", "IND")
        conf_ind = self.get_stream_poi("confluences", "IND")
        all_ind = np.append(head_ind, conf_ind)
        del conf_ind, head_ind # Clean up
        
        # We created a zeros arrays and put in confluences and heads their id
        # Those id will be consecutive numbers starting in one
        seg_arr = np.zeros(self._ncells, dtype=np.int32)
        for n, inds in enumerate(all_ind):
            seg_arr[inds] = n+1
        
        # Move throught channel list. If receiver is 0, give receiver the same id that giver.
        # If a receiver is not 0, that means that we are in a confluence. 
        for n in range(len(self._ix)):
            if seg_arr[self._ixc[n]] == 0:
                seg_arr[self._ixc[n]] = seg_arr[self._ix[n]]
        
        # Reshape and output
        seg_arr = seg_arr.reshape(self._dims)
        if asgrid:
            return self._create_output_grid(seg_arr, 0)
        else:
            return seg_arr        

    def get_stream_order(self, kind="strahler", asgrid=True):
        """
        This function extract streams orderded by strahler or shreeve. Cell values
        will have a value acording with the order of the segment they belong
        ===========
        Parameters
        
        kind : *Str*, {'strahler', 'shreeve'}
        asgrid : *Boolean*, indicates if the selfwork is returned as *TopoGrid* (True) or as a *np.array*
        """
        if kind not in ['strahler', 'shreeve']:
            return
        
        # Get grid channel cells
        str_ord = np.zeros(self._ncells, dtype=np.int8)
        str_ord[self._ix] = 1
        str_ord[self._ixc] = 1
        visited = np.zeros(self._ncells, dtype=np.int8)
    
        if kind == 'strahler':
            for n in range(len(self._ix)):
                if (str_ord[self._ixc[n]] == str_ord[self._ix[n]]) & visited[self._ixc[n]]:
                    str_ord[self._ixc[n]] = str_ord[self._ixc[n]] + 1 
                else:
                    str_ord[self._ixc[n]] = max(str_ord[self._ix[n]], str_ord[self._ixc[n]])
                    visited[self._ixc[n]] = True
        elif kind == 'shreeve':
            for n in range(len(self._ix)):
                if visited[self._ixc[n]]:
                    str_ord[self._ixc[n]] = str_ord[self._ixc[n]] + str_ord[self._ix[n]]
                else:
                    str_ord[self._ixc[n]] = max(str_ord[self._ix[n]], str_ord[self._ixc[n]])
                    visited[self._ixc[n]] = True
        str_ord = str_ord.reshape(self._dims)
        
        if asgrid:
            return self._create_output_grid(str_ord, nodata_value=0)
        else:
            return str_ord        
        
    def get_klargestconncomps(self, k=1, asgrid=True):
        """
        This function return k largest connected components in an instance of *Stream*
        =========== 
        Parameters
        
        k : *Int*, number of largest connected components 
        asgrid : *Boolean*, indicates if the network is returned as *TopoGrid* (True) or as a *np.ndarray* (False)
        ===========   
        Returns
        
        cc_arr : *TopoGrid* or numpy.ndarray with the k largest connected components
        cc_id  : *Numpy.array*, The id or ids of the connected components
        ===========   
        References
        
        This function functions similarly with the klargestconncomps.m @STREAMobj from Topotoolbox 
        matlab codes developed by Wolfgang Schwanghart (version of 17. August, 2017), but uses different method.
        """
        temp_ix = self._ix
        temp_ixc = self._ixc
        ncc = 0
        cc_arr = np.zeros(self._ncells, np.int)
        nix = len(temp_ix)
        for n in range(nix-1,-1,-1):
            # If receiver is zero, add a new cc
            if cc_arr[temp_ixc[n]] == 0:
                ncc += 1
                cc_arr[temp_ixc[n]] = ncc  
            # Mark cc giver with the id of the cc receiver
            cc_arr[temp_ix[n]] = cc_arr[temp_ixc[n]]
        del temp_ix,temp_ixc
        cc_k = cc_arr.max() 
        if k > cc_k:
            k = cc_k
            print('There are only %d connected components in the stream network' % k)
        cc_size = np.zeros((cc_k,2))    
        for i in range(0,cc_k):
            cc_id = i+1
            cc_size[i,0] = cc_id
            cc_size[i,1] = np.sum(cc_arr == cc_id)
        cc_size = cc_size[np.lexsort(-cc_size.T)]
        cc_size = cc_size[0:k]
        # del the smaller cc
        for i in range(0,cc_k):
            cc_id = i+1
            if cc_id not in cc_size[:,0]:
                cc_arr[np.where(cc_arr == cc_id)]=0
                
        # Reshape and return
        cc_arr = cc_arr.reshape(self._dims)  
        
        if asgrid:
            return self._create_output_grid(cc_arr, 0), cc_size[:,0].astype(int)
        else:
            return cc_arr, cc_size[:,0].astype(int)        

    def get_trunk(self,ccs_arr,ccs_id):
        """
        This function extract trunk stream (longest stream) of the connected component (cc),
        by Strahler stream order and the size of stream nodes number
        =========== 
        Parameters
        
        cc_arr : *Numpy.ndarray*, array indicate the k largest connected components with value = cc_id
        cc_id  : *Numpy.array*, The IDs of the connected components
        ===========   
        Returns

        trunk_arr : *Numpy.ndarray*, array indicate the trunk with value = k
        =========== 
        Example
    
        stream = Stream(dem=topo, flow=flow, threshold=threshold, verbose=False, thetaref=0.45, npoints=5)
        ccs_arr,ccs_id = stream.get_klargestconncomps(k=5,asgrid=False)
        trunk_arr = stream.get_trunk(ccs_arr,ccs_id)
        ===========     
        References
        
        This function functions similarly with the trunk.m @STREAMobj from Topotoolbox 
        matlab codes developed by Wolfgang Schwanghart (version of 17. August, 2017), but uses different method.
        Their algorithm identifies the trunk by sequently tracing the maximum downstream distance in upstream direction. 
        """ 
        ncells = self._ncells
        trunk_arr =  np.zeros(ncells, dtype=np.int8)
        for k in range(0,len(ccs_id)):
            cc_id = ccs_id[k]
        
            w = ccs_arr == cc_id
            w = w.ravel()
            I   = w[self._ix]
            ixc = self._ixc[I]     # ixc of this cc
            ix = self._ix[I]       # ixc of this cc
 
            # 1. calculate the strahler stream order of this cc
            str_ord = np.zeros(ncells, dtype=np.int8)
            str_ord[ix] = 1
            str_ord[ixc] = 1
            visited = np.zeros(ncells, dtype=np.int8)
    
            for n in range(len(ix)):
                if (str_ord[ixc[n]] == str_ord[ix[n]]) & visited[ixc[n]]:
                    str_ord[ixc[n]] = str_ord[ixc[n]] + 1 
                else:
                    str_ord[ixc[n]] = max(str_ord[ix[n]], str_ord[ixc[n]])
                    visited[ixc[n]] = True    
        
            norder =  str_ord.max()
            l_norder = []
            for n in range(len(ix)):
                if str_ord[ix[n]] == norder:
                    l_norder.append(n)  
        
            if norder == 1:
                l_trunk = l_norder
            else:
                l_trunk = l_norder
                for iorder in range(norder,1,-1):                
                    l_lssiorder = get_largerpart_lssorder(ncells, ix, ixc ,iorder,l_norder, str_ord)
                    l_norder = l_lssiorder
                    l_trunk = l_lssiorder + l_trunk
                
            trunk_arr[ix[l_trunk]] = k+1      
        trunk_arr = trunk_arr.reshape(self._dims)   
        return trunk_arr        
    def get_trunk_output(self,trunk_arr,path):
        """
        This function output the coordinates of the trunk stream 
        =========== 
        Parameters
        
        trunk_arr : *Numpy.ndarray*, array indicate the trunk with value = k
        path : *Str*, path for the output text file
        ===========   
        Returns

        text file: x, y, z, distance
        =========== 
        Example
    
        stream = Stream(dem=topo, flow=flow, threshold=threshold, verbose=False, thetaref=0.45, npoints=5)
        ccs_arr,ccs_id = stream.get_klargestconncomps(k=5,asgrid=False)
        trunk_arr = stream.get_trunk(ccs_arr,ccs_id)
        stream.get_trunk_output(trunk_arr)
        ===========     
        """ 
        cab = "x;y;z;distance"
        n_trunk = trunk_arr.max()
        for n in range(1,n_trunk+1):
            w = trunk_arr == n
            w = w.ravel()
            I = w[self._ix]
            ixc = self._ixc[I] # ixc of this trunk
            ix = self._ix[I]   # ixc of this trunk
            row,col = self.ind_2_cell(ix)
            x, y = self.cell_2_xy(row, col)
            out_arr = np.array((x, y, self._zx[I], self._dx[I])).T
            fname = path + 'river' + str(n) +'.txt'
            np.savetxt(fname, out_arr,fmt='%3.8f %3.8f %3.3f %3.3f', header=cab) 
               
def get_largerpart_lssorder(ncells, ix, ixc ,iorder,l_iorder, str_ord):
    """
    This function gets the part with strahler stream order (i-1) in the trunk,
    by comparing the nodes number of all less i-order streams in the connected component (cc)
    =========== 
    Parameters
        
    ncells : *Int*, number of cells in grid
    ix  : *Numpy.array*, ix of this cc
    ixc : *Numpy.array*, ixc of this cc
    iorder : *Int*, the strahler stream order number
    l_iorder : *List*, the index in ix of the part with order i in the trunk
    str_ord : *Numpy.ndarray*, stream order of the grid nodes
    ===========   
    Returns

    l_lssi : *List*, the index in ix of the part with order (i-1) in the trunk
    """     
    head_id = ix[l_iorder[0]]
    
    # 2. Find the indexs (in ix) of parts with all smaller orders 
    # and these parts' outlet is the head found above
    tem_arr =  np.zeros(ncells, np.int)
    tem_arr[head_id] = 1
    nix = len(ix)
    l_lssi_all = []
    for n in range(nix-1,-1,-1):
        # If the stream receiver is not zero and stream giver is zero
        if (tem_arr[ixc[n]] != 0) & (tem_arr[ix[n]] == 0): 
            # Mark giver with the id of the receiver
            tem_arr[ix[n]] = tem_arr[ixc[n]]
            l_lssi_all.append(n)
    l_lssi_all = l_lssi_all[::-1]
    
    # 3. Sort the parts with their size (nodes number) 
    part_arr = np.zeros(ncells, np.int)   
    npart = 0
    nl_all = len(l_lssi_all)
    for n in range(nl_all-1,-1,-1):
        # If receiver's order is i, add a new part
        id_i = l_lssi_all[n]
        if str_ord[ixc[id_i]] == iorder:
            npart += 1
            part_arr[ixc[id_i]] = npart                
         # Mark giver with the id of the receiver
        part_arr[ix[id_i]] = part_arr[ixc[id_i]]                    
 
    npart = int(npart)
    if npart > 1:
        part_size = np.zeros((npart,2),np.int)
        for i in range(0,npart):
            part_id = i+1
            part_size[i,0] = part_id 
            part_size[i,1] = np.sum(part_arr == part_id)
        part_size = part_size[np.lexsort(-part_size.T)] 
        keep_id = part_size[0,0]
    else:
        keep_id = 1
        
    # 4. Output the list of index in ix of the part with order (i-1)    
    l_lssi = []
    lss_iorder = iorder - 1
    for n in range(0,nix):
        id_arr = ix[n]
        if (str_ord[id_arr] == lss_iorder) & (part_arr[id_arr] == keep_id):
            l_lssi.append(n)
    return l_lssi
