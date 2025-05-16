# -------------------------------------------------------------------------
# Name:        Lakes and reservoirs module
# Purpose:
#
# Author:      PB
#
# Created:     01/08/2016
# Copyright:   (c) PB 2016
# -------------------------------------------------------------------------

from cwatm.management_modules.data_handling import *
from cwatm.hydrological_modules.routing_reservoirs.routing_sub import *

from cwatm.management_modules.globals import *
import importlib

class lakes_reservoirs(object):
    """
    LAKES AND RESERVOIRS

    Note:

        Calculate water retention in lakes and reservoirs 

        Using the **Modified Puls approach** to calculate retention of a lake
        See also: LISFLOOD manual Annex 3 (Burek et al. 2013)

        for Modified Puls Method the Q(inflow)1 has to be used. It is assumed that this is the same as Q(inflow)2 for the first timestep
        has to be checked if this works in forecasting mode!

        Lake Routine using Modified Puls Method (see Maniak, p.331ff)

        .. math::
             {Qin1 + Qin2 \\over{2}} - {Qout1 + Qout2 \\over{2}} = {S2 - S1 \\over{\\delta time}}

        changed into:

        .. math::
             {S2 \\over{time + Qout2/2}} = {S1 \\over{dtime + Qout1/2}} - Qout1 + {Qin1 + Qin2 \\over{2}}

        Outgoing discharge (Qout) are linked to storage (S) by elevation.

        Now some assumption to make life easier:

        1.) storage volume is increase proportional to elevation: S = A * H where: H: elevation, A: area of lake

        2.) :math:`Q_{\\mathrm{out}} = c * b * H^{2.0}` (c: weir constant, b: width)

             2.0 because it fits to a parabolic cross section see (Aigner 2008) (and it is much easier to calculate (that's the main reason)

        c: for a perfect weir with mu=0.577 and Poleni: :math:`{2 \\over{3}} \\mu * \\sqrt{2*g} = 1.7`

        c: for a parabolic weir: around 1.8 

        because it is a imperfect weir: :math:`C = c * 0.85 = 1.5`

        results in formular: :math:`Q = 1.5 * b * H^2 = a*H^2 -> H = \\sqrt{Q/a}`

        Solving the equation:

        :math:`{S2 \\over{dtime + Qout2/2}} = {S1 \\over{dtime + Qout1/2}} - Qout1 + {Qin1 + Qin2 \\over{2}}`

        :math:`SI = {S2 \\over{dtime}} + {Qout2 \\over{2}} = {A*H \\over{DtRouting}} + {Q \\over{2}} = {A \\over{DtRouting*\\sqrt{a}* \\sqrt{Q}}} + {Q \\over{2}}`

        -> replacement: :math:`{A \\over{DtSec * \\sqrt{a}}} = Lakefactor, Y = \\sqrt{Q}`

        :math:`Y^2 + 2 * Lakefactor *Y - 2 * SI=0`

        solution of this quadratic equation:

        :math:`Q = (-LakeFactor + \\sqrt{LakeFactor^2+2*SI})^2`


    **Global variables**

    =====================================  ======================================================================  =====
    Variable [self.var]                    Description                                                             Unit 
    =====================================  ======================================================================  =====
    modflow                                Flag: True if modflow_coupling = True in settings file                  --   
    load_initial                           Settings initLoad holds initial conditions for variables                input
    wastewater_to_reservoirs                                                                                       --   
    saveInit                               Flag: if true initial conditions are saved                              --   
    waterBodyID                            lakes/reservoirs map with a single ID for each lake/reservoir           --   
    waterBodyOut                           biggest outlet (biggest accumulation of ldd network) of a waterbody     --   
    dirUp                                  river network in upstream direction                                     --   
    ldd_LR                                 change river network (put pits in where lakes are)                      --   
    lddCompress_LR                         compressed river network lakes/reservoirs (without missing values)      --   
    dirUp_LR                               river network direction upstream lake/reservoirs                        --   
    dirupLen_LR                            number of bifurcation upstream lake/reservoir                           --   
    dirupID_LR                             index river upstream lake/reservoir                                     --   
    downstruct_LR                          river network downstream lake/reservoir                                 --   
    catchment_LR                           catchments lake/reservoir                                               --   
    dirDown_LR                             river network direktion downstream lake/reservoir                       --   
    lendirDown_LR                          number of river network connections lake/reservoir                      --   
    compress_LR                            boolean map as mask map for compressing lake/reservoir                  --   
    decompress_LR                          boolean map as mask map for decompressing lake/reservoir                --   
    waterBodyOutC                          compressed map biggest outlet of each lake/reservoir                    --   
    waterBodyID_C                                                                                                  --   
    resYear                                Settings waterBodyYear, with first operating year of reservoirs         map  
    resYearC                               Compressed map of resYear                                               --   
    waterBodyTyp                           Settings, waterBodyTyp, with waterbody type 1-4                         map  
    waterBodyTyp_unchanged                                                                                         --   
    includeType4                           True if there is a reservoir of waterbody type 4 in waterBodyTyp map    bool 
    waterBodyTypC                          water body types 3 reservoirs and lakes (used as reservoirs but before  --   
    resVolumeC                             compressed map of reservoir volume                                      Milli
    resId_restricted                                                                                               --   
    waterBodyBuffer                                                                                                --   
    waterBodyBuffer_wwt                                                                                            --   
    lakeArea                               area of each lake/reservoir                                             m2   
    lakeAreaC                              compressed map of the area of each lake/reservoir                       m2   
    lakeDis0                               compressed map average discharge at the outlet of a lake/reservoir      m3/s 
    lakeDis0C                              average discharge at the outlet of a lake/reservoir                     m3/s 
    lakeAC                                 compressed map of parameter of channel width, gravity and weir coeffic  --
    reservoir_transfers_net_M3             net reservoir transfers, after exports, transfers, and imports          m3
    reservoir_transfers_in_M3              water received into reservoirs                                          m3
    reservoir_transfers_out_M3             water given from reservoirs                                             m3
    resVolumeOnlyReservoirs                                                                                        --
    resVolumeOnlyReservoirsC                                                                                       --   
    resVolume                                                                                                      --   
    lakeEvaFactorC                         compressed map of a factor which increases evaporation from lake becau  --   
    reslakeoutflow                                                                                                 --   
    lakeVolume                             volume of lakes                                                         m3   
    outLake                                outflow from lakes                                                      m    
    lakeInflow                                                                                                     --   
    lakeOutflow                                                                                                    --   
    reservoirStorage                                                                                               --   
    MtoM3C                                 conversion factor from m to m3 (compressed map)                         --   
    EvapWaterBodyMOutlet                   Evaporation from lakes and reservoirs summed up at outlet               m
    EvapWaterBodyM                         Evaporation from lakes and reservoirs                                   m
    lakeResInflowM                                                                                                 --   
    lakeResOutflowM                                                                                                --   
    lakedaycorrect                                                                                                 --   
    lakeFactor                             factor for the Modified Puls approach to calculate retention of the la  --   
    lakeFactorSqr                          square root factor for the Modified Puls approach to calculate retenti  --   
    lakeInflowOldC                         inflow to the lake from previous days                                   m/3  
    lakeOutflowC                           compressed map of lake outflow                                          m3/s 
    lakeLevelC                             compressed map of lake level                                            m    
    conLimitC                                                                                                      --   
    normLimitC                                                                                                     --   
    floodLimitC                                                                                                    --   
    adjust_Normal_FloodC                                                                                           --   
    norm_floodLimitC                                                                                               --   
    minQC                                                                                                          --   
    normQC                                                                                                         m3/s
    nondmgQC                                                                                                       m3/s
    deltaO                                                                                                         --   
    deltaLN                                                                                                        --   
    deltaLF                                                                                                        --   
    deltaNFL                                                                                                       --   
    reservoirFillC                                                                                                 --   
    waterBodyTypCTemp                                                                                              --   
    waterBodyTypTemp                                                                                               --   
    sumEvapWaterBodyC                                                                                              --   
    sumlakeResInflow                                                                                               --   
    sumlakeResOutflow                                                                                              --   
    lakeResStorage_release_ratio                                                                                   --   
    lakeResStorage_release_ratioC                                                                                  --   
    lakeIn                                                                                                         --   
    lakeEvapWaterBodyC                                                                                             --   
    resEvapWaterBodyC                                                                                              --   
    downstruct                                                                                                     --
    lakeStorage                                                                                                    --   
    resStorage                                                                                                     --   
    cellArea                               Area of cell                                                            m2   
    DtSec                                  number of seconds per timestep (default = 86400)                        s    
    MtoM3                                  Coefficient to change units                                             --   
    InvDtSec                                                                                                       --   
    UpArea1                                upstream area of a grid cell                                            m2   
    lddCompress                            compressed river network (without missing values)                       --   
    lakeEvaFactor                          a factor which increases evaporation from lake because of wind          --   
    dtRouting                              number of seconds per routing timestep                                  s    
    evapWaterBodyC                         Compressed version of EvapWaterBodyM                                    m    
    sumLakeEvapWaterBodyC                                                                                          --   
    noRoutingSteps                                                                                                 --   
    sumResEvapWaterBodyC                                                                                           --   
    discharge                              Channel discharge                                                       m3/s 
    inflowDt                                                                                                       --   
    prelakeResStorage                                                                                              --   
    runoff                                                                                                         --   
    includeWastewater                                                                                              --   
    lakeVolumeM3C                          compressed map of lake volume                                           m3
    lakeStorageC                                                                                                   --   
    reservoirStorageM3C                                                                                            --   
    lakeResStorageC                                                                                                --   
    lakeResStorage                                                                                                 --   
    =====================================  ======================================================================  =====

    **Functions**
    """

    def __init__(self, model):
        self.var = model.var
        self.model = model

    def reservoir_releases(self, xl_settings_file_path):
        pd = importlib.import_module("pandas", package=None)
        df = pd.read_excel(xl_settings_file_path, sheet_name='Reservoirs_downstream')
        waterBodyID_C_tolist = self.var.waterBodyID_C.tolist()

        reservoir_release = [[-1 for i in self.var.waterBodyID_C] for i in range(366)]
        for res in list(df)[2:]:
            if res in waterBodyID_C_tolist:
                res_index = waterBodyID_C_tolist.index(int(float(res)))

                for day in range(366):
                    reservoir_release[day][res_index] = df[res][day]

        reservoir_supply = [[-1 for i in self.var.waterBodyID_C] for i in range(366)]
        #reservoir_release.copy()
        if 'Reservoirs_supply' in pd.read_excel(xl_settings_file_path, None).keys():
            df2 = pd.read_excel(xl_settings_file_path, sheet_name='Reservoirs_supply')
            for res in list(df2)[2:]:
                if res in waterBodyID_C_tolist:
                    res_index = waterBodyID_C_tolist.index(int(float(res)))

                    for day in range(366):
                        reservoir_supply[day][res_index] = df2[res][day]
        else:
            reservoir_supply = reservoir_release.copy()
        
        return reservoir_release, reservoir_supply


    def initWaterbodies(self):
        """
        Initialize water bodies
        Read parameters from maps e.g
        area, location, initial average discharge, type 9reservoir or lake) etc.

        Compress numpy array from mask map to the size of lakes+reservoirs
        (marked as capital C at the end of the variable name)
        """

        def buffer_waterbody(rec, waterBody):
            """
            Puts a buffer of a rectangular rec around the lakes and reservoirs
            parameter rec = size of rectangular
            output buffer = compressed buffer
            """

            # add waterBody as input - allow to create a buffer on customized maps
            # waterBody = decompress(self.var.waterBodyID)
            rows, cols = waterBody.shape
            buffer = np.full((rows, cols), 1.0e15)
            for y in range(rows):
                for x in range(cols):
                    id = waterBody[y, x]
                    if id > 0:
                        for j in range(1, rec + 1):
                            addj = j // 2
                            if j % 2: addj = -addj
                            for i in range(1, rec + 1):
                                addi = i // 2
                                if i % 2: addi = -addi
                                yy = y + addj
                                xx = x + addi
                                if yy >= 0 and yy < rows and xx >= 0 and xx < cols:
                                    if id < buffer[yy, xx]:
                                        buffer[yy, xx] = id
            buffer[buffer == 1.0e15] = 0.

            # In the case of overlapping buffers, we ensure that each waterbody
            # is at least inside its own command area, and not the CA of another waterbody
            buffer = np.where(waterBody > 0, waterBody, buffer)

            return compressArray(buffer).astype(np.int64)



        if checkOption('includeWaterBodies'):

            # load lakes/reservoirs map with a single ID for each lake/reservoir
            self.var.waterBodyID = loadmap('waterBodyID').astype(np.int64)

            # look into the Excel to find if there are new reservoirs
            if checkOption('reservoir_add_info_in_Excel',True):
                resnew = decompress(globals.inZero.copy())
                resnew1 = []
                remove = []

                for i in range(len(self.var.reservoir_info)):

                    if self.var.reservoir_info[i][1]:
                        # create new reservoir/lake
                        col = int((float(self.var.reservoir_info[i][3]) - maskmapAttr['x']) * maskmapAttr['invcell'])
                        row = int((maskmapAttr['y'] - float(self.var.reservoir_info[i][2])) * maskmapAttr['invcell'])
                        if (col<0) or (row<0) or (col> (resnew.shape[1]-1)) or  (row> (resnew.shape[0]-1)):

                            msg = "---- New lake/reservoir No: " + str(int(self.var.reservoir_info[i][0])) + " with coordinates:\n"
                            msg +=  "lat or y: "+ str(self.var.reservoir_info[i][2]) + "  lon or x: " + str(self.var.reservoir_info[i][3]) + "\n"
                            msg += "is not in Mask map (result in row/col " + str(row) + " " + str(col) + ")\n"
                            if Flags['loud']: print (msg)
                            remove.append(self.var.reservoir_info[i])

                        else:
                            resnew[row,col] = int(self.var.reservoir_info[i][0])
                            resnew1.append(int(self.var.reservoir_info[i][0]))

                if len(resnew1) > 0:
                    resnewC = compressArray(resnew).astype(np.int64)
                    # check if lakes/res is in map
                    remove1 = []
                    for i in resnew1:
                        if not(i in resnewC):
                            msg = "--- New lake/reservoir No: "+ str(i) + " is not in the Mask Map\n"
                            if Flags['loud']: print (msg)

                    # check if lake/res already exists
                    index = np.where(resnewC>0)[0]
                    for i in index:
                        if self.var.waterBodyID[i] > 0:
                            no = resnewC[i]
                            msg = "New lake/reservoir No: "+ str(resnewC[i]) + \
                                  " is at the place of an existing one No: "+ str(self.var.waterBodyID[i]) + "\n"
                            raise CWATMError(msg)
                    self.var.waterBodyID = np.where(self.var.waterBodyID == 0, resnewC, self.var.waterBodyID)

                for i in remove:
                    self.var.reservoir_info.remove(i)


            self.var.includeWastewater = False
            # if checkoption 2nd argument = True then check first for argument ,if False print an error if not in setttings
            self.var.includeWastewater = checkOption('includeWastewater',True)
                
            # calculate biggest outlet = biggest accumulation of ldd network
            lakeResmax = npareamaximum(self.var.UpArea1, self.var.waterBodyID)
            self.var.waterBodyOut = np.where(self.var.UpArea1 == lakeResmax, self.var.waterBodyID, 0)

            # dismiss water bodies that are not a subcatchment of an outlet
            sub = subcatchment1(self.var.dirUp, self.var.waterBodyOut, self.var.UpArea1)
            self.var.waterBodyID = np.where(self.var.waterBodyID == sub, sub, 0)

            

            # and again calculate outlets, because ID might have changed due to the operation before
            lakeResmax = npareamaximum(self.var.UpArea1, self.var.waterBodyID)
            self.var.waterBodyOut = np.where(self.var.UpArea1 == lakeResmax, self.var.waterBodyID, 0)

            # change ldd: put pits in where lakes are:
            self.var.ldd_LR = np.where(self.var.waterBodyID > 0, 5, self.var.lddCompress)

            # create new ldd without lakes reservoirs
            self.var.lddCompress_LR, dirshort_LR, self.var.dirUp_LR, self.var.dirupLen_LR, self.var.dirupID_LR, \
            self.var.downstruct_LR, self.var.catchment_LR, self.var.dirDown_LR, self.var.lendirDown_LR = defLdd2(
                self.var.ldd_LR)

            # boolean map as mask map for compressing and decompressing
            self.var.compress_LR = self.var.waterBodyOut > 0
            self.var.decompress_LR = np.nonzero(self.var.waterBodyOut)[0]
            self.var.waterBodyOutC = np.compress(self.var.compress_LR, self.var.waterBodyOut)
            self.var.waterBodyID_C = np.compress(self.var.compress_LR, self.var.waterBodyID)

            # First year that the reservoir is operating
            self.var.resYear = loadmap('waterBodyYear')
            self.var.resYearC = np.compress(self.var.compress_LR, self.var.resYear)
            self.var.resVolumeC = np.compress(self.var.compress_LR, loadmap('waterBodyVolRes')) * 1000000
            self.var.waterBodyTyp = loadmap('waterBodyTyp').astype(np.int64)
            self.var.waterBodyTypC = np.compress(self.var.compress_LR, self.var.waterBodyTyp)

            # ================================
            # Lakes

            # Surface area of each lake [m2]
            self.var.lakeArea = loadmap('waterBodyArea') * 1000 * 1000  # mult with 1000000 to convert from km2 to m2
            self.var.lakeAreaC = np.compress(self.var.compress_LR, self.var.lakeArea)

            # lake discharge at outlet to calculate alpha: parameter of channel width, gravity and weir coefficient
            # Lake parameter A (suggested  value equal to outflow width in [m])
            # multiplied with the calibration parameter LakeMultiplier
            self.var.lakeDis0 = np.maximum(loadmap('waterBodyDis'), 0.1)
            self.var.lakeDis0C = np.compress(self.var.compress_LR, self.var.lakeDis0)
            chanwidth = 7.1 * np.power(self.var.lakeDis0C, 0.539)
            lakeAFactor = globals.inZero + loadmap('lakeAFactor')
            lakeAFactorC = np.compress(self.var.compress_LR, lakeAFactor)
            self.var.lakeAC = lakeAFactorC * 0.612 * 2 / 3 * chanwidth * (2 * 9.81) ** 0.5

            # ================================
            # Reservoirs

            self.var.reservoir_transfers_net_M3C = np.compress(self.var.compress_LR, globals.inZero.copy())
            self.var.reservoir_transfers_net_M3 = globals.inZero.copy()
            self.var.reservoir_transfers_in_M3C = np.compress(self.var.compress_LR, globals.inZero.copy())
            self.var.reservoir_transfers_in_M3 = globals.inZero.copy()
            self.var.reservoir_transfers_out_M3C = np.compress(self.var.compress_LR, globals.inZero.copy())
            self.var.reservoir_transfers_out_M3 = globals.inZero.copy()

            # a factor which increases evaporation from lake because of wind
            self.var.lakeEvaFactor = globals.inZero + loadmap('lakeEvaFactor')
            self.var.lakeEvaFactorC = np.compress(self.var.compress_LR, self.var.lakeEvaFactor)

            # initial only the big arrays are set to 0,
            # the  initial values are loaded inside the subroutines of lakes and reservoirs
            self.var.reslakeoutflow = globals.inZero.copy()
            self.var.lakeVolume = globals.inZero.copy()
            self.var.outLake = self.var.load_initial("outLake")

            self.var.lakeStorage = globals.inZero.copy()
            self.var.lakeInflow = globals.inZero.copy()
            self.var.lakeOutflow = globals.inZero.copy()
            self.var.reservoirStorage = globals.inZero.copy()

            self.var.MtoM3C = np.compress(self.var.compress_LR, self.var.MtoM3)

            # init water balance [m]
            self.var.EvapWaterBodyMOutlet = globals.inZero.copy()
            self.var.lakeResInflowM = globals.inZero.copy()
            self.var.lakeResOutflowM = globals.inZero.copy()

            if checkOption('calcWaterBalance'):
                self.var.lakedaycorrect = globals.inZero.copy()

            if 'reservoir_add_info_in_Excel' in option:
                if checkOption('reservoir_add_info_in_Excel'):
                    for i in range(len(self.var.reservoir_info)):
                        resint = int(self.var.reservoir_info[i][0])
                        resindex = np.where(self.var.waterBodyID_C == resint)
                        # test if reservoir is found -> later on maqke a new one
                        if resindex[0].size > 0:
                            resindex = resindex[0].tolist()[0]
                            if int(self.var.reservoir_info[i][4]) >0: self.var.waterBodyTypC[resindex] = int(self.var.reservoir_info[i][4])

                            if float(self.var.reservoir_info[i][6]) >0: self.var.lakeAreaC[resindex] = float(self.var.reservoir_info[i][6]) * 1000 * 1000
                            if float(self.var.reservoir_info[i][7]) > 0: self.var.lakeDis0C[resindex] = float(self.var.reservoir_info[i][7])
                            if float(self.var.reservoir_info[i][8]) > 0: self.var.resVolumeC[resindex] = float(self.var.reservoir_info[i][8]) * 1000000
                            if float(self.var.reservoir_info[i][9]) > 0:
                                lakeAFactorC = float(self.var.reservoir_info[i][9])
                                chanwidth = 7.1 * np.power(self.var.lakeDis0C[resindex], 0.539)
                                self.var.lakeAC[resindex] = lakeAFactorC * 0.612 * 2 / 3 * chanwidth * (2 * 9.81) ** 0.5

                            if float(self.var.reservoir_info[i][10]) >0: self.var.lakeEvaFactorC[resindex] = float(self.var.reservoir_info[i][10])
                            if float(self.var.reservoir_info[i][11]) >0: self.var.resYearC[resindex] = int(self.var.reservoir_info[i][11])

            # back to lakeArea , because it is used in routing_kinematic
            np.put(self.var.lakeArea, self.var.decompress_LR, self.var.lakeAreaC)
            # correcting reservoir volume for lakes, just to run them all as reservoirs
            self.var.resVolumeC = np.where(self.var.resVolumeC > 0, self.var.resVolumeC, self.var.lakeAreaC * 10)
            self.var.resVolume = globals.inZero.copy()
            np.put(self.var.resVolume, self.var.decompress_LR, self.var.resVolumeC)

            # Flag if res type-4 are used
            self.var.includeType4 = False
            if (self.var.waterBodyTypC == 4).any():
                self.var.includeType4 = True

            np.put(self.var.waterBodyTyp, self.var.decompress_LR, self.var.waterBodyTypC)
            self.var.waterBodyTyp_unchanged = self.var.waterBodyTyp.copy()
            self.var.waterBodyTypC = np.where(self.var.waterBodyOutC > 0, self.var.waterBodyTypC.astype(np.int64), 0)

            # changing reservoirs type 2 and 3 to lakes if volumes are zero:
            self.var.waterBodyTypC = np.where(self.var.resVolumeC > 0., self.var.waterBodyTypC,
                                              np.where(self.var.waterBodyTypC == 2, 1, self.var.waterBodyTypC))
            self.var.waterBodyTypC = np.where(self.var.resVolumeC > 0., self.var.waterBodyTypC,
                                              np.where(self.var.waterBodyTypC == 3, 1, self.var.waterBodyTypC))

            np.put(self.var.waterBodyTyp, self.var.decompress_LR, self.var.waterBodyTypC)

            # needs to be changed - maybe update all reservoir to spreadsheet
            self.var.resId_restricted = globals.inZero.copy()
            if self.var.includeWastewater:
                resIdList = np.unique(list(self.var.wastewater_to_reservoirs.values()))
                for rid in resIdList:
                    self.var.resId_restricted += np.where(self.var.waterBodyID == rid, rid, 0)

            # Create a buffer around water bodies as command areas for lakes and reservoirs
            if checkOption('includeWaterDemand'):
                waterBody_UnRestricted = self.var.waterBodyID.copy()

                if self.var.includeWastewater:
                    waterBody_UnRestricted = np.where(np.in1d(waterBody_UnRestricted, self.var.resId_restricted), 0, waterBody_UnRestricted)
                rectangular = 1
                if "buffer_waterbodies" in binding:
                    rectangular = int(loadmap('buffer_waterbodies'))

                self.var.waterBodyBuffer = buffer_waterbody(rectangular, decompress(waterBody_UnRestricted))
                if self.var.includeWastewater:
                    self.var.waterBodyBuffer_wwt = buffer_waterbody(rectangular, decompress(self.var.resId_restricted))

                # rectangular = 1
                # if "buffer_waterbodies" in binding:
                #    rectangular = int(loadmap('buffer_waterbodies'))
                # self.var.waterBodyBuffer = buffer_waterbody(rectangular)
            else:
                self.var.waterBodyBuffer = self.var.waterBodyID.copy()




    def initial_lakes(self):
        """
        Initial part of the lakes module
        Using the **Modified Puls approach** to calculate retention of a lake
        """

        # self_.var.lakeInflowOldC = np.bincount(self_.var.downstruct, weights=self.var.ChanQ)[self.var.LakeIndex]

        # self.var.lakeInflowOldC = np.compress(self_.var.compress_LR, self_.var.chanQKin)
        # for Modified Puls Method the Q(inflow)1 has to be used.
        # It is assumed that this is the same as Q(inflow)2 for the first timestep
        # Does this work in forecasting mode?

        self.var.lakeFactor = self.var.lakeAreaC / (self.var.dtRouting * np.sqrt(self.var.lakeAC))
        self.var.lakeFactorSqr = np.square(self.var.lakeFactor)
        # for faster calculation inside dynamic section

        lakeInflowIni = self.var.load_initial("lakeInflow")  # inflow in m3/s estimate
        if not (isinstance(lakeInflowIni, np.ndarray)):
            self.var.lakeInflowOldC = self.var.lakeDis0C.copy()
        else:
            self.var.lakeInflowOldC = np.compress(self.var.compress_LR, lakeInflowIni)
        # if a new lake is not initialized use assumption calculation
        self.var.lakeInflowOldC = np.where(self.var.lakeInflowOldC > 0,self.var.lakeInflowOldC,self.var.lakeDis0C)

        lakeVolumeIni = self.var.load_initial("lakeStorage")
        if not (isinstance(lakeVolumeIni, np.ndarray)):
            self.var.lakeVolumeM3C = self.var.lakeAreaC * np.sqrt(self.var.lakeInflowOldC / self.var.lakeAC)
        else:
            self.var.lakeVolumeM3C = np.compress(self.var.compress_LR, lakeVolumeIni)

        # if a new lake ist not initialized use assumption
        self.var.lakeVolumeM3C = np.where(self.var.lakeVolumeM3C > 0, self.var.lakeVolumeM3C,self.var.lakeAreaC * np.sqrt(self.var.lakeInflowOldC / self.var.lakeAC))
        self.var.lakeStorageC = self.var.lakeVolumeM3C.copy()

        lakeOutflowIni = self.var.load_initial("lakeOutflow")
        lakeStorageIndicator = np.maximum(0.0, self.var.lakeVolumeM3C / self.var.dtRouting + self.var.lakeInflowOldC / 2)
        # SI = S/dt + Q/2
        lakeOutflowC1 = np.square(-self.var.lakeFactor + np.sqrt(self.var.lakeFactorSqr + 2 * lakeStorageIndicator))
        # solution of quadratic equation
        #  it is as easy as this because:
        # 1. storage volume is increase proportional to elevation
        #  2. Q= a *H **2.0  (if you choose Q= a *H **1.5 you have to solve the formula of Cardano)
        if not (isinstance(lakeOutflowIni, np.ndarray)):
            self.var.lakeOutflowC = lakeOutflowC1.copy()
        else:
            self.var.lakeOutflowC = np.compress(self.var.compress_LR, lakeOutflowIni)
        # lake storage ini
        self.var.lakeOutflowC = np.where(self.var.lakeOutflowC>0,self.var.lakeOutflowC,lakeOutflowC1)

        self.var.lakeLevelC = self.var.lakeVolumeM3C / self.var.lakeAreaC
        ii =1

    def initial_reservoirs(self):
        """
        Initial part of the reservoir module
        Using the appraoch of LISFLOOD

        See Also:
            LISFLOOD manual Annex 1: (Burek et al. 2013)
        """
        #  Conservative, normal and flood storage limit (fraction of total storage, [-])
        self.var.conLimitC = np.compress(self.var.compress_LR, loadmap('conservativeStorageLimit') + globals.inZero)
        self.var.normLimitC = np.compress(self.var.compress_LR, loadmap('normalStorageLimit') + globals.inZero)
        self.var.floodLimitC = np.compress(self.var.compress_LR, loadmap('floodStorageLimit') + globals.inZero)

        # Minimum, Normal and Non-damaging reservoir outflow  (fraction of average discharge, [-])
        # multiplied with the given discharge at the outlet from Hydrolakes database
        self.var.minQC = np.compress(self.var.compress_LR, loadmap('MinOutflowQ') * self.var.lakeDis0)
        self.var.normQC = np.compress(self.var.compress_LR, loadmap('NormalOutflowQ') * self.var.lakeDis0)
        self.var.nondmgQC = np.compress(self.var.compress_LR, loadmap('NonDamagingOutflowQ') * self.var.lakeDis0)
        self.var.adjust_Normal_FloodC = np.compress(self.var.compress_LR,loadmap('adjust_Normal_Flood') + globals.inZero)

        if 'reservoir_add_info_in_Excel' in option:
            if checkOption('reservoir_add_info_in_Excel'):
               for i in range(len(self.var.reservoir_info)):
                    resint = int(self.var.reservoir_info[i][0])
                    resindex = np.where(self.var.waterBodyID_C == resint)
                    # test if reservoir is found -> later on maqke a new one
                    if resindex[0].size > 0:
                        resindex = resindex[0].tolist()[0]
                        if float(self.var.reservoir_info[i][12]) > 0: self.var.conLimitC[resindex] = float(self.var.reservoir_info[i][12])
                        if float(self.var.reservoir_info[i][13]) > 0: self.var.normLimitC[resindex] = float(self.var.reservoir_info[i][13])
                        if float(self.var.reservoir_info[i][14]) > 0: self.var.floodLimitC[resindex] = float(self.var.reservoir_info[i][14])
                        if float(self.var.reservoir_info[i][15]) > 0: self.var.adjust_Normal_FloodC[resindex] = float(self.var.reservoir_info[i][15])

                        if float(self.var.reservoir_info[i][16]) > 0: self.var.minQC[resindex] = float(self.var.reservoir_info[i][16]) * self.var.lakeDis0C[resindex]
                        if float(self.var.reservoir_info[i][17]) > 0: self.var.normQC[resindex] = float(self.var.reservoir_info[i][17]) * self.var.lakeDis0C[resindex]
                        if float(self.var.reservoir_info[i][18]) > 0: self.var.nondmgQC[resindex] = float(self.var.reservoir_info[i][18]) * self.var.lakeDis0C[resindex]

        self.var.norm_floodLimitC = self.var.normLimitC + self.var.adjust_Normal_FloodC * (
                    self.var.floodLimitC - self.var.normLimitC)

        # Repeatedly used expressions in reservoirs routine
        self.var.deltaO = self.var.normQC - self.var.minQC
        self.var.deltaLN = self.var.normLimitC - 2 * self.var.conLimitC
        self.var.deltaLF = self.var.floodLimitC - self.var.normLimitC
        self.var.deltaNFL = self.var.floodLimitC - self.var.norm_floodLimitC

        reservoirStorageIni = self.var.load_initial("reservoirStorage")

        self.var.reservoirFillC = self.var.normLimitC.copy()
        # Initial reservoir fill (fraction of total storage, [-])
        self.var.reservoirStorageM3C = self.var.reservoirFillC * self.var.resVolumeC
        # set initial fill of distribution reservoir to zero (res. type-4)
        loadres = self.var.reservoirStorageM3C * 0.

        if isinstance(reservoirStorageIni, np.ndarray):
            loadres = np.compress(self.var.compress_LR, reservoirStorageIni)

        self.var.reservoirStorageM3C = np.where(loadres == 0., self.var.reservoirStorageM3C, loadres)
        # for waterbodytyp 4 and 5
        self.var.reservoirStorageM3C = np.where(self.var.waterBodyTypC > 3, 0., self.var.reservoirStorageM3C)
        self.var.reservoirFillC = self.var.reservoirStorageM3C / self.var.resVolumeC

        # water balance
        self.var.lakeResStorageC = np.where(self.var.waterBodyTypC == 0, 0.,
                                            np.where(self.var.waterBodyTypC == 1, self.var.lakeStorageC,
                                                     self.var.reservoirStorageM3C))
        lakeStorageC = np.where(self.var.waterBodyTypC == 1, self.var.lakeStorageC, 0.)
        resStorageC = np.where(self.var.waterBodyTypC > 1, self.var.reservoirStorageM3C, 0.)
        self.var.lakeResStorage = globals.inZero.copy()
        self.var.lakeStorage = globals.inZero.copy()
        self.var.resStorage = globals.inZero.copy()
        np.put(self.var.lakeResStorage, self.var.decompress_LR, self.var.lakeResStorageC)
        np.put(self.var.lakeStorage, self.var.decompress_LR, lakeStorageC)
        np.put(self.var.resStorage, self.var.decompress_LR, resStorageC)

        self.var.reservoir_releases_excel_option = False
        if checkOption('reservoir_releases_in_Excel_settings', True):
            if 'Excel_settings_file' in binding:
                self.var.reservoir_releases_excel_option = True
                xl_settings_file_path = cbinding('Excel_settings_file')
                self.var.reservoir_releases, self.var.reservoir_supply = self.reservoir_releases(xl_settings_file_path)

                self.var.reservoir_releases = np.array(self.var.reservoir_releases)
                self.var.reservoir_supply = np.array(self.var.reservoir_supply)

        #  check if transfer has a valid receiver or giver
        if checkOption('reservoir_transfers', True):
            remove = []
            for trans in self.var.reservoir_transfers:
                if not(trans[1] in self.var.waterBodyID):
                    msg = "Reservoir transfer: giving reservoir is missing in Excel: " + str(trans[1])
                    if Flags['loud']: print (msg)
                    remove.append(trans)

            for trans in self.var.reservoir_transfers:
                if not(trans[2] in self.var.waterBodyID):
                    msg = "Reservoir transfer: receiving reservoir is missing in Excel: " + str(trans[1])
                    if Flags['loud']: print (msg)
                    remove.append(trans)

            for i in remove:
                try:
                    self.var.reservoir_transfers.remove(i)
                except:
                    ii =0  # transfer seems to be in both receiving and giving missing
            ii =1


    # ------------------ End init ------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------

    def dynamic(self):
        """
        Dynamic part set lakes and reservoirs for each year
        """
        if checkOption('includeWaterBodies'):
            # check years
            if dateVar['newStart'] or dateVar['newYear']:
                year = dateVar['currDate'].year

                # water body types:
                # - 5 = Fake reservoir: outflow = inflow, no evaporation -> used for water transfer
                # - 4 = water distribution reservoirs, input = output; only allows input from specified sources,
                #       e.g., inter-basin transfers, wastewater treatment, desalination, etc.
                # - 3 = reservoirs and lakes (used as reservoirs but before the year of construction as lakes)
                # - 2 = reservoirs (regulated discharge)
                # - 1 = lakes (weirFormula)
                # - 0 = non lakes or reservoirs (e.g. wetland)

                if returnBool('useResAndLakes'):
                    if returnBool('dynamicLakesRes'):
                        year = dateVar['currDate'].year
                    else:
                        year = loadmap('fixLakesResYear')

                    self.var.waterBodyTypCTemp = np.where((self.var.resYearC > year) & (self.var.waterBodyTypC == 2), 0,
                                                          self.var.waterBodyTypC)
                    self.var.waterBodyTypCTemp = np.where((self.var.resYearC > year) & (self.var.waterBodyTypC > 3), 0,
                                                          self.var.waterBodyTypCTemp)
                    self.var.waterBodyTypCTemp = np.where((self.var.resYearC > year) & (self.var.waterBodyTypC == 3), 1,
                                                          self.var.waterBodyTypCTemp)
                    # we also need the uncompressed version to compute leakage
                    if self.var.modflow or self.var.includeType4:
                        self.var.waterBodyTypTemp = np.where((self.var.resYear > year) & (self.var.waterBodyTyp == 2),
                                                             0, self.var.waterBodyTyp)
                        self.var.waterBodyTypTemp = np.where((self.var.resYear > year) & (self.var.waterBodyTyp > 3),
                                                             0, self.var.waterBodyTypTemp)
                        self.var.waterBodyTypTemp = np.where((self.var.resYear > year) & (self.var.waterBodyTyp == 3),
                                                             1, self.var.waterBodyTypTemp)
                else:
                    self.var.waterBodyTypCTemp = np.where(self.var.waterBodyTypC == 2, 0, self.var.waterBodyTypC)
                    self.var.waterBodyTypCTemp = np.where(self.var.waterBodyTypC > 3, 0, self.var.waterBodyTypCTemp)
                    self.var.waterBodyTypCTemp = np.where(self.var.waterBodyTypC == 3, 1, self.var.waterBodyTypCTemp)

                    if self.var.modflow or self.var.includeType4:
                        self.var.waterBodyTypTemp = np.where(self.var.waterBodyTyp == 2, 0, self.var.waterBodyTyp)
                        self.var.waterBodyTypTemp = np.where(self.var.waterBodyTyp > 3, 0, self.var.waterBodyTypTemp)
                        self.var.waterBodyTypTemp = np.where(self.var.waterBodyTyp == 3, 1, self.var.waterBodyTypTemp)

            self.var.sumEvapWaterBodyC = 0
            self.var.sumlakeResInflow = 0
            self.var.sumlakeResOutflow = 0

            # Reservoir_releases holds variables for different reservoir operations, each with 366 timesteps.
            # The value of the variable at the reservoir is the maximum fraction of available storage to be
            # released for the associated operation.
            # Downstream release is the water released downstream into the river.
            # This is overridden only in flooding conditions.

            if 'Reservoir_releases' in binding:
                day_of_year = dateVar['currDate'].timetuple().tm_yday
                self.var.lakeResStorage_release_ratio = readnetcdf2(
                    'Reservoir_releases', day_of_year,
                    useDaily='DOY', value='Downstream release')

                self.var.lakeResStorage_release_ratioC = np.compress(self.var.compress_LR,
                                                                     self.var.lakeResStorage_release_ratio)

            elif self.var.reservoir_releases_excel_option:
                self.var.lakeResStorage_release_ratioC = self.var.reservoir_releases[dateVar['doy']-1]

    def dynamic_inloop(self, NoRoutingExecuted):
        """
        Dynamic part to calculate outflow from lakes and reservoirs

        * lakes with modified Puls approach
        * reservoirs with special filling levels

        :param NoRoutingExecuted: actual number of routing substep
        :return: outLdd: outflow in m3 to the network

        Note:
            outflow to adjected lakes and reservoirs is calculated separately
        """

        def dynamic_inloop_lakes(inflowC, NoRoutingExecuted):
            """
            Lake routine to calculate lake outflow

            :param inflowC: inflow to lakes and reservoirs [m3]
            :param NoRoutingExecuted: actual number of routing substep
            :return: QLakeOutM3DtC - lake outflow in [m3] per subtime step
            """

            # ************************************************************
            # ***** LAKE
            # ************************************************************

            if checkOption('calcWaterBalance'):
                #    ii = 3
                oldlake = self.var.lakeStorageC.copy()

            # if (dateVar['curr'] == 3):

            # Lake inflow in [m3/s]
            lakeInflowC = inflowC / self.var.dtRouting

            # just for day to day waterbalance -> get X as difference
            # lakeIn = in + X ->  (in + old) * 0.5 = in + X  ->   in + old = 2in + 2X -> in - 2in +old = 2x
            # -> (old - in) * 0.5 = X
            lakedaycorrectC = 0.5 * (
                        inflowC / self.var.dtRouting - self.var.lakeInflowOldC) * self.var.dtRouting  # [m3]

            self.var.lakeIn = (lakeInflowC + self.var.lakeInflowOldC) * 0.5
            # for Modified Puls Method: (S2/dtime + Qout2/2) = (S1/dtime + Qout1/2) - Qout1 + (Qin1 + Qin2)/2
            #  here: (Qin1 + Qin2)/2

            self.var.lakeEvapWaterBodyC = np.where((self.var.lakeVolumeM3C - self.var.evapWaterBodyC) > 0.,
                                                   self.var.evapWaterBodyC, self.var.lakeVolumeM3C)
            self.var.sumLakeEvapWaterBodyC += self.var.lakeEvapWaterBodyC
            self.var.lakeVolumeM3C = self.var.lakeVolumeM3C - self.var.lakeEvapWaterBodyC
            # lakestorage - evaporation from lakes

            self.var.lakeInflowOldC = lakeInflowC.copy()
            # Qin2 becomes Qin1 for the next time step [m3/s]

            lakeStorageIndicator = np.maximum(0.0,
                                              self.var.lakeVolumeM3C / self.var.dtRouting - 0.5 * self.var.lakeOutflowC + self.var.lakeIn)
            # here S1/dtime - Qout1/2 + LakeIn , so that is the right part of the equation above

            self.var.lakeOutflowC = np.square(
                -self.var.lakeFactor + np.sqrt(self.var.lakeFactorSqr + 2 * lakeStorageIndicator))

            QLakeOutM3DtC = self.var.lakeOutflowC * self.var.dtRouting
            # Outflow in [m3] per timestep

            # New lake storage [m3] (assuming lake surface area equals bottom area)
            self.var.lakeVolumeM3C = (lakeStorageIndicator - self.var.lakeOutflowC * 0.5) * self.var.dtRouting
            # Lake storage

            self.var.lakeStorageC += self.var.lakeIn * self.var.dtRouting - QLakeOutM3DtC - self.var.lakeEvapWaterBodyC

            if self.var.noRoutingSteps == (NoRoutingExecuted + 1):
                self.var.lakeLevelC = self.var.lakeVolumeM3C / self.var.lakeAreaC

            # expanding the size
            # self.var.QLakeOutM3Dt = globals.inZero.copy()
            # np.put(self.var.QLakeOutM3Dt,self.var.LakeIndex,QLakeOutM3DtC)
            # if  (self.var.noRoutingSteps == (NoRoutingExecuted + 1)):
            if self.var.saveInit and (self.var.noRoutingSteps == (NoRoutingExecuted + 1)):
                np.put(self.var.lakeVolume, self.var.decompress_LR, self.var.lakeVolumeM3C)
                np.put(self.var.lakeInflow, self.var.decompress_LR, self.var.lakeInflowOldC)
                np.put(self.var.lakeOutflow, self.var.decompress_LR, self.var.lakeOutflowC)

            # Water balance
            if self.var.noRoutingSteps == (NoRoutingExecuted + 1):
                np.put(self.var.lakeStorage, self.var.decompress_LR, self.var.lakeStorageC)

            if checkOption('calcWaterBalance'):
                self.model.waterbalance_module.waterBalanceCheck(
                    [self.var.lakeIn],  # In [m3/s]
                    [self.var.lakeOutflowC, self.var.lakeEvapWaterBodyC / self.var.dtRouting],
                    # Out  self.var.evapWaterBodyC
                    [oldlake / self.var.dtRouting],  # prev storage
                    [self.var.lakeStorageC / self.var.dtRouting],
                    "lake", False)

            if checkOption('calcWaterBalance'):
                np.put(self.var.lakedaycorrect, self.var.decompress_LR, lakedaycorrectC)
                self.model.waterbalance_module.waterBalanceCheck(
                    [inflowC / self.var.dtRouting],  # In [m3/s]
                    [self.var.lakeOutflowC, self.var.lakeEvapWaterBodyC / self.var.dtRouting,
                     lakedaycorrectC / self.var.dtRouting],  # Out  self.var.evapWaterBodyC
                    [oldlake / self.var.dtRouting],  # prev storage
                    [self.var.lakeStorageC / self.var.dtRouting],
                    "lake2", False)

            if checkOption('calcWaterBalance'):
                self.model.waterbalance_module.waterBalanceCheck(
                    [inflowC],  # In [m3/s]
                    [QLakeOutM3DtC, self.var.lakeEvapWaterBodyC, lakedaycorrectC],  # Out  self.var.evapWaterBodyC
                    [oldlake],  # prev storage
                    [self.var.lakeStorageC],
                    "lake3", False)

            return QLakeOutM3DtC

        # ---------------------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------------------

        def dynamic_inloop_reservoirs(inflowC, NoRoutingExecuted):
            """
            Reservoir outflow

            :param inflowC: inflow to reservoirs [m3] per subtime step
            :param NoRoutingExecuted: actual number of routing substep
            :return: qResOutM3DtC - reservoir outflow in [m3] per subtime step
            """

            # ************************************************************
            # ***** Reservoirs
            # ************************************************************

            if checkOption('calcWaterBalance'):
                oldres = self.var.reservoirStorageM3C.copy()

            # QResInM3Dt = inflowC
            # Reservoir inflow in [m3] per timestep
            #self.var.reservoirStorageM3C += inflowC
            # New reservoir storage [m3] = plus inflow for this sub step
            # if reservoirtype = 4 -> no inflow
            self.var.reservoirStorageM3C = np.where(self.var.waterBodyTypC == 4, self.var.reservoirStorageM3C, self.var.reservoirStorageM3C + inflowC)

            # check that reservoir storage - evaporation > 0
            self.var.resEvapWaterBodyC = np.where(self.var.reservoirStorageM3C - self.var.evapWaterBodyC > 0.,
                                                  self.var.evapWaterBodyC, self.var.reservoirStorageM3C)

            # Reservoir type 5 is only in for water transfer without a reservoir retention without ET
            self.var.resEvapWaterBodyC = np.where(self.var.waterBodyTypC == 5, 0., self.var.resEvapWaterBodyC)
            self.var.sumResEvapWaterBodyC += self.var.resEvapWaterBodyC

            self.var.reservoirStorageM3C = self.var.reservoirStorageM3C - self.var.resEvapWaterBodyC

            self.var.reservoirFillC = self.var.reservoirStorageM3C / self.var.resVolumeC
            # New reservoir fill [-]

            # MODIFIED DOR FRIDMAN 
            # limit res. outflows to reservoir with water level > limit_to_resOutflows
            # limit_to_resOutflows is defined relative to res. Volume
            if "limit_to_resOutflows" in binding:
                reservoirOutflowLimit = loadmap('limit_to_resOutflows')
                # load map of outflow limit np.compress(self.var.compress_LR, self.var.waterBodyOut)
                reservoirOutflowLimitC = np.compress(self.var.compress_LR,
                                                     globals.inZero.copy() + reservoirOutflowLimit)
                # compress to lake&res data
                reservoirOutflowLimitMask = np.where(reservoirOutflowLimitC < self.var.reservoirFillC, 1, 0)
                # Create mask for limit out flow: 1 allow out flow (in case fill > limit) 
            else:
                reservoirOutflowLimitMask = np.compress(self.var.compress_LR, globals.inZero.copy() + 1)

            reservoirOutflow1 = np.minimum(self.var.minQC, self.var.reservoirStorageM3C * self.var.InvDtSec)
            # Reservoir outflow [m3/s] if ReservoirFill is nearing absolute minimum. 

            reservoirOutflow2 = self.var.minQC + self.var.deltaO * (
                        self.var.reservoirFillC - 2 * self.var.conLimitC) / self.var.deltaLN
            # 2*ConservativeStorageLimit
            # Reservoir outflow [m3/s] if NormalStorageLimit <= ReservoirFill > 2*ConservativeStorageLimit

            reservoirOutflow3 = self.var.normQC + (
                        (self.var.reservoirFillC - self.var.norm_floodLimitC) / self.var.deltaNFL) * (
                                            self.var.nondmgQC - self.var.normQC)
            # Reservoir outflow [m3/s] if FloodStorageLimit le ReservoirFill gt NormalStorageLimit
            inflowCQ = inflowC * self.var.InvDtSec * self.var.noRoutingSteps

            temp = np.minimum(self.var.nondmgQC,np.maximum(inflowCQ, self.var.normQC))
            reservoirOutflow4 = np.maximum(
                (self.var.reservoirFillC - self.var.floodLimitC - 0.01) * self.var.resVolumeC * self.var.InvDtSec, temp)

            # Reservoir outflow [m3/s] if ReservoirFill gt FloodStorageLimit
            # Depending on ReservoirFill the reservoir outflow equals ReservoirOutflow1, ReservoirOutflow2,
            # ReservoirOutflow3 or ReservoirOutflow4

            reservoirOutflow = reservoirOutflow1.copy()

            reservoirOutflow = np.where(self.var.reservoirFillC > 2 * self.var.conLimitC, reservoirOutflow2,
                                        reservoirOutflow)
            reservoirOutflow = np.where(self.var.reservoirFillC > self.var.normLimitC, self.var.normQC,
                                        reservoirOutflow)
            reservoirOutflow = np.where(self.var.reservoirFillC > self.var.norm_floodLimitC, reservoirOutflow3,
                                        reservoirOutflow)
            reservoirOutflow = np.where(self.var.reservoirFillC > self.var.floodLimitC, reservoirOutflow4,
                                        reservoirOutflow)

            temp = np.minimum(reservoirOutflow, np.maximum(inflowCQ, self.var.normQC))

            reservoirOutflow = np.where((reservoirOutflow > 1.2 * inflowCQ) &
                                        (reservoirOutflow > self.var.normQC) &
                                        (self.var.reservoirFillC < self.var.floodLimitC), temp, reservoirOutflow)

            # Reservoir_releases holds variables for different reservoir operations, each with 366 timesteps.
            # The value of the variable at the reservoir is the maximum fraction of available storage to be
            # released for the associated operation.
            # Downstream release is the water released downstream into the river.
            # This is overridden only in flooding conditions.

            if 'Reservoir_releases' in binding:
                reservoirOutflow = np.where(self.var.reservoirFillC > self.var.floodLimitC, reservoirOutflow,
                                            np.where(self.var.lakeResStorage_release_ratioC > 0,
                                            self.var.lakeResStorage_release_ratioC * self.var.reservoirStorageM3C *
                                            self.var.InvDtSec, 0))

            if self.var.reservoir_releases_excel_option:
                reservoirOutflow = np.where(self.var.reservoirFillC > self.var.floodLimitC, reservoirOutflow,
                                            np.where(self.var.lakeResStorage_release_ratioC > -1,
                                            self.var.lakeResStorage_release_ratioC * self.var.reservoirStorageM3C *
                                            self.var.InvDtSec, reservoirOutflow))

            # Apply reservoirOutflowLimitMask to limit resOutflows based on storagre relative to threshold volume
            # This may override the reservoir_releases behaviour - ATTENTION FROM MIKHAIL/DOR
            reservoirOutflow = reservoirOutflow * reservoirOutflowLimitMask

            qResOutM3DtC = reservoirOutflow * self.var.dtRouting

            # Reservoir outflow in [m3] per sub step
            qResOutM3DtC = np.where(self.var.reservoirStorageM3C - qResOutM3DtC > 0., qResOutM3DtC,
                                    self.var.reservoirStorageM3C)
            # check if storage would go < 0 if outflow is used
            qResOutM3DtC = np.maximum(qResOutM3DtC, self.var.reservoirStorageM3C - self.var.resVolumeC)
            # Check to prevent reservoir storage from exceeding total capacity
            # for watertyp 4: In = Out
            qResOutM3DtC = np.where(self.var.waterBodyTypC == 4, inflowC, qResOutM3DtC)
            # for watertyp 5: In = Out. Outflow = inflow + water that is transferred
            qResOutM3DtC = np.where(self.var.waterBodyTypC == 5, self.var.reservoirStorageM3C, qResOutM3DtC)

            self.var.reservoirStorageM3C -= qResOutM3DtC

            # Transfer water between reservoirs
            # Send storage between reservoirs using the Excel sheet reservoir_transfers within cwatm_settings.xlsx
            # Using the waterBodyIDs defined in the settings, designate
            # the Giver, the Receiver, and the daily fraction of live storage the Giver sends to the Receiver.
            # If the Receiver is already at capacity, the Giver does not send any storage.
            # Reservoirs can only send to one reservoir. Reservoirs can receive from several reservoirs.

            # for in and out rewservoirs the old storage volume is stored. The difference between new and old is added to outflow
            #self.var.oldReservoirStM3C = self.var.reservoirStorageM3C.copy()

            if checkOption('reservoir_transfers',True):

                inZero_C = np.compress(self.var.compress_LR, globals.inZero.copy())
                for transfer in self.var.reservoir_transfers:

                    if returnBool('dynamicLakesRes'):
                        year = dateVar['currDate'].year
                    else:
                        year = loadmap('fixLakesResYear')

                    if transfer[1] > 0:
                        # using Giving reservoir ID from Excel - if this is 0 it is from outside
                        giver = np.where(self.var.waterBodyID_C == transfer[1])[0][0]
                        giver_already_constructed = self.var.resYearC[giver] <= year
                    else:
                        giver_already_constructed = True

                    if transfer[2] > 0:
                        # using the Receiving reservoir from Excel
                        receiver = np.where(self.var.waterBodyID_C == transfer[2])[0][0]
                        receiver_already_constructed = self.var.resYearC[receiver] <= year
                    else:
                        receiver_already_constructed = True


                    if receiver_already_constructed and giver_already_constructed:
                        # if giving and receiving station already exist (is build before the year which is modelled)
                        if (transfer[2] == 0) or (self.var.waterBodyTypC[receiver] > 3):
                            # if receiver is outside OR the receiving is a reservoir type > 3
                            reservoir_unused_receiver = 10e12
                        else:
                            # check if there is space for the water in the receiving reservoir
                            reservoir_unused_receiver = self.var.resVolumeC[receiver] - self.var.reservoirStorageM3C[receiver]

                        if transfer[1] > 0:
                            reservoir_storage_giver = self.var.reservoirStorageM3C[giver]
                        else:
                            # In this case, the fraction refers to the fraction of the receiver,
                            # as the giver is infinite
                            reservoir_storage_giver = self.var.resVolumeC[receiver]



                        # --- Rulesets -------------------
                        # use different rulesets:
                        limit =  transfer[3]
                        # Rule 1: Fraction of live storage (values  ≤ 1) (default)
                        if transfer[0] == 1:
                            reservoir_transfer_actual = reservoir_storage_giver * transfer[4][dateVar['doy'] - 1] / self.var.noRoutingSteps

                        # Rule 2: Fraction of live storage (values  ≤ 1) e, LIMIT: miminum reservoir volume which should be preserved
                        if transfer[0] == 2:
                            reservoir_transfer_actual = reservoir_storage_giver * transfer[4][dateVar['doy'] - 1] / self.var.noRoutingSteps
                            if (self.var.reservoirStorageM3C[giver] - reservoir_transfer_actual) < limit:
                                reservoir_transfer_actual = np.maximum(0, self.var.reservoirStorageM3C[giver] - limit)

                        # Rule 3: m3 of storage volume, LIMIT: miminum reservoir volume which should be preserved
                        if transfer[0] == 3:
                            reservoir_transfer_actual = transfer[4][dateVar['doy'] - 1] / self.var.noRoutingSteps
                            if (self.var.reservoirStorageM3C[giver] - reservoir_transfer_actual) < limit:
                                reservoir_transfer_actual = np.maximum(0, self.var.reservoirStorageM3C[giver] - limit)


                        # if rule is based on volume:
                        if transfer[0] < 4:
                            reservoir_transfer_actual = np.minimum(reservoir_unused_receiver * 0.95,reservoir_transfer_actual)

                        # --- Outflow -------------------------------
                        # Outflow based rules
                        # Rule 4: Fraction of outflow (values  ≤ 1) Limit: minimum discharge
                        if transfer[0] == 4:
                            limit = limit * self.var.dtRouting
                            if transfer[1] == 0:  # giver is outside
                                limit = 0
                            # to get from m3/s to m3 of the substep routing poeriod
                            qnew = qResOutM3DtC[giver] * (1 - transfer[4][dateVar['doy'] - 1])  # [m3]
                            if qnew < limit:
                                qnew = np.minimum(qResOutM3DtC[giver], limit)
                            reservoir_transfer_actual = qResOutM3DtC[giver] - qnew
                            qResOutM3DtC[giver] = qnew

                        # Rule 5: Fraction of outflow (v alues  ≤ 1)
                        # Limit: minimum discharge, maximum discharge to receiving reservoir [m3/s]
                        if transfer[0] == 5:
                            limit= transfer[3].split(",")
                            limit1 = float(limit[0]) * self.var.dtRouting  # m3/s to m3 in the routing substep period
                            limit2 = float(limit[1]) * self.var.dtRouting
                            if transfer[1] == 0:  # giver is outside
                                limit1 = 0
                            # to get to m3 from m3/s but including routing steps
                            #mult = self.var.DtSec / self.var.noRoutingSteps
                            qnew = qResOutM3DtC[giver]  * (1 - transfer[4][dateVar['doy'] - 1]) # in m3
                            if qnew < limit1:
                                qnew = np.minimum(qResOutM3DtC[giver] , limit1)
                            if (qResOutM3DtC[giver] - qnew) < limit2:
                                reservoir_transfer_actual = qResOutM3DtC[giver] - qnew
                            else:
                                reservoir_transfer_actual = limit2
                                qnew = qResOutM3DtC[giver] - reservoir_transfer_actual
                            qResOutM3DtC[giver] = qnew


                        # Rule 6: Fix of outflow, [m3/s]
                        # Limit: minimum discharge [m3/s]
                        if transfer[0] == 6:
                            limit = limit * self.var.dtRouting
                            if transfer[1] == 0:  # giver is outside
                                limit = 0
                            reservoir_transfer_actual = transfer[4][dateVar['doy'] - 1]  * self.var.dtRouting
                            # from m3/s to m3
                            qnew = qResOutM3DtC[giver] - reservoir_transfer_actual
                            if qnew < limit:
                                qnew = np.minimum(qResOutM3DtC[giver], limit)
                            reservoir_transfer_actual = qResOutM3DtC[giver] - qnew
                            qResOutM3DtC[giver] = qnew


                        # -----------------------------------------

                        if transfer[1] > 0:
                            # There is a giver, not the ocean
                            inZero_C[giver] = -reservoir_transfer_actual
                            self.var.reservoir_transfers_out_M3C[giver] += reservoir_transfer_actual

                            if transfer[0] < 4:
                                # if given from storage than (not if split from outflow)
                                self.var.reservoirStorageM3C[giver] = self.var.reservoirStorageM3C[giver] - reservoir_transfer_actual

                        if transfer[2] > 0:  # There is a receiver, not the ocean
                            inZero_C[receiver] = reservoir_transfer_actual  # receiver
                            self.var.reservoirStorageM3C[receiver] = self.var.reservoirStorageM3C[receiver] + reservoir_transfer_actual
                            self.var.reservoir_transfers_in_M3C[receiver] += reservoir_transfer_actual

                    #print(transfer[0], 'donated', reservoir_transfer_actual, 'm3 to', transfer[1])

                self.var.reservoir_transfers_net_M3C += inZero_C

            #--------------------


            # New reservoir storage [m3]
            self.var.reservoirStorageM3C = np.maximum(0.0, self.var.reservoirStorageM3C)
            # New reservoir fill
            self.var.reservoirFillC = self.var.reservoirStorageM3C / self.var.resVolumeC

            # if  (self.var.noRoutingSteps == (NoRoutingExecuted + 1)):
            if self.var.noRoutingSteps == (NoRoutingExecuted + 1):
                np.put(self.var.reservoirStorage, self.var.decompress_LR, self.var.reservoirStorageM3C)

            if checkOption('calcWaterBalance'):
                self.model.waterbalance_module.waterBalanceCheck(
                    [inflowC / self.var.dtRouting],  # In
                    [qResOutM3DtC / self.var.dtRouting, self.var.resEvapWaterBodyC / self.var.dtRouting],
                    # Out  self.var.evapWaterBodyC
                    [oldres / self.var.dtRouting],  # prev storage
                    [self.var.reservoirStorageM3C / self.var.dtRouting],
                    "res1", False)

            return qResOutM3DtC

        # ---------------------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------------------
        # lake and reservoirs

        if checkOption('calcWaterBalance'):
            prereslake = self.var.lakeResStorageC.copy()
            prelake = self.var.lakeStorageC.copy()

        # ----------
        # inflow lakes
        # 1.  dis = upstream1(self_.var.downstruct_LR, self_.var.discharge)   # from river upstream
        # 2.  runoff = npareatotal(self.var_.waterBodyID, self_.var.waterBodyID)  # from cell itself
        # 3.                  # outflow from upstream lakes

        # ----------
        # outflow lakes res -> inflow ldd_LR
        # 1. out = upstream1(self_.var.downstruct, self_.var.outflow)

        # collect discharge from above waterbodies
        dis_LR = upstream1(self.var.downstruct_LR, self.var.discharge)
        # only where lakes are and unit convered to [m]
        dis_LR = np.where(self.var.waterBodyID > 0, dis_LR, 0.) * self.var.DtSec

        # sum up runoff and discharge on the lake
        inflow = npareatotal(dis_LR + self.var.runoff * self.var.cellArea, self.var.waterBodyID)

        # only once at the outlet
        inflow = np.where(self.var.waterBodyOut > 0, inflow, 0.) / self.var.noRoutingSteps + self.var.outLake

        if checkOption('inflow'):
            # if inflow ( from module inflow) goes to a lake this is not counted, because lakes,reservoirs are dislinked from the network
            inflow2basin = npareatotal(self.var.inflowDt, self.var.waterBodyID)
            inflow2basin = np.where(self.var.waterBodyOut > 0, inflow2basin, 0.)
            inflow = inflow + inflow2basin

        # calculate total inflow into lakes and compress it to waterbodie outflow point
        # inflow to lake is discharge from upstream network + runoff directly into lake + outflow from upstream lakes
        inflowC = np.compress(self.var.compress_LR, inflow)

        # ------------------------------------------------------------
        self.var.resEvapWaterBodyC = globals.inZero.copy()
        outflowLakesC = dynamic_inloop_lakes(inflowC, NoRoutingExecuted)
        outflowResC = dynamic_inloop_reservoirs(inflowC, NoRoutingExecuted)
        outflow0C = inflowC.copy()  # no retention
        outflowC = np.where(self.var.waterBodyTypCTemp == 0, outflow0C,
                            np.where(self.var.waterBodyTypCTemp == 1, outflowLakesC, outflowResC))

        # outflowC =  outflowLakesC        # only lakes
        # outflowC = outflowResC
        # outflowC = inflowC.copy() - self.var.evapWaterBodyC
        # outflowC = inflowC.copy()

        # waterbalance
        inflowCorrC = np.where(self.var.waterBodyTypCTemp == 1, self.var.lakeIn * self.var.dtRouting, inflowC)
        # EvapWaterBodyC = np.where( self.var.waterBodyTypCTemp == 0, 0. , np.where( self.var.waterBodyTypCTemp == 1, self.var.sumLakeEvapWaterBodyC, self.var.sumResEvapWaterBodyC))
        EvapWaterBodyC = np.where(self.var.waterBodyTypCTemp == 0, 0.,
                                  np.where(self.var.waterBodyTypCTemp == 1, self.var.lakeEvapWaterBodyC,
                                           self.var.resEvapWaterBodyC))

        self.var.lakeResStorageC = np.where(self.var.waterBodyTypCTemp == 0, 0.,
                                            np.where(self.var.waterBodyTypCTemp == 1, self.var.lakeStorageC,
                                                     self.var.reservoirStorageM3C))
        lakeStorageC = np.where(self.var.waterBodyTypCTemp == 1, self.var.lakeStorageC, 0.)
        resStorageC = np.where(self.var.waterBodyTypCTemp > 1, self.var.reservoirStorageM3C, 0.)

        self.var.sumEvapWaterBodyC += EvapWaterBodyC  # in [m3]
        self.var.sumlakeResInflow += inflowCorrC
        self.var.sumlakeResOutflow += outflowC
        # self.var.sumlakeResOutflow = self.var.sumlakeResOutflow  + self.var.lakeOutflowC * self.var.dtRouting

        # decompress to normal maskarea size waterbalance
        if self.var.noRoutingSteps == (NoRoutingExecuted + 1):
            np.put(self.var.EvapWaterBodyMOutlet, self.var.decompress_LR, self.var.sumEvapWaterBodyC)
            self.var.EvapWaterBodyMOutlet = self.var.EvapWaterBodyMOutlet / self.var.cellArea
            np.put(self.var.lakeResInflowM, self.var.decompress_LR, self.var.sumlakeResInflow)
            self.var.lakeResInflowM = self.var.lakeResInflowM / self.var.cellArea
            np.put(self.var.lakeResOutflowM, self.var.decompress_LR, self.var.sumlakeResOutflow)
            self.var.lakeResOutflowM = self.var.lakeResOutflowM / self.var.cellArea

            # --------------- correction PB May 2024
            # calculate evaporation for each cell of the lake
            # each lake cell fraction of the total lake area
            # a lake cell has at minimum 5% water
            fracwatermin = np.where(self.var.waterBodyID > 0, np.maximum(self.var.fracVegCover[5],0.05),0)
            wlakefracsum = npareatotal(fracwatermin, self.var.waterBodyID)
            # -> part of each cell of the total lake -> sum for each lake = 1
            wlakefrac = divideValues(self.var.fracVegCover[5], wlakefracsum)
            # all lake id cells get the evaporation of the outlet cell
            ebody = npareatotal(self.var.EvapWaterBodyMOutlet, self.var.waterBodyID)
            # 3) step evoporation is distributed by the water frac of each lake
            self.var.EvapWaterBodyM = ebody * wlakefrac
            self.var.EvapWaterBodyM[np.isnan(self.var.EvapWaterBodyM)] = 0.

            np.put(self.var.lakeResStorage, self.var.decompress_LR, self.var.lakeResStorageC)
            np.put(self.var.lakeStorage, self.var.decompress_LR, lakeStorageC)
            np.put(self.var.resStorage, self.var.decompress_LR, resStorageC)

            #water transfer
            if checkOption('reservoir_transfers', True):
                np.put(self.var.reservoir_transfers_net_M3, self.var.decompress_LR, self.var.reservoir_transfers_net_M3C)
                np.put(self.var.reservoir_transfers_out_M3, self.var.decompress_LR, self.var.reservoir_transfers_out_M3C)
                np.put(self.var.reservoir_transfers_in_M3, self.var.decompress_LR, self.var.reservoir_transfers_in_M3C)
                self.var.reservoir_transfers_net_M3C = np.compress(self.var.compress_LR,globals.inZero.copy())
                self.var.reservoir_transfers_out_M3C = np.compress(self.var.compress_LR, globals.inZero.copy())
                self.var.reservoir_transfers_in_M3C = np.compress(self.var.compress_LR, globals.inZero.copy())


        # ------------------------------------------------------------

        np.put(self.var.reslakeoutflow, self.var.decompress_LR, outflowC)
        lakeResOutflowDis = npareatotal(self.var.reslakeoutflow, self.var.waterBodyID) / (
                    self.var.DtSec / self.var.noRoutingSteps)
        # shift outflow 1 cell downstream
        out1 = upstream1(self.var.downstruct, self.var.reslakeoutflow)
        # everything with is not going to another lake is output to river network
        outLdd = np.where(self.var.waterBodyID > 0, 0, out1)

        # everything what is not going to the network is going to another lake
        outLake1 = np.where(self.var.waterBodyID > 0, out1, 0)
        # sum up all inflow from other lakes
        outLakein = npareatotal(outLake1, self.var.waterBodyID)
        # use only the value of the outflow point
        self.var.outLake = np.where(self.var.waterBodyOut > 0, outLakein, 0.)

        # Puts the value of lakeResStorage into all cells covered by the waterbody
        self.var.lakeResStorage_filled = npareamaximum(self.var.lakeResStorage, self.var.waterBodyID)
        self.var.lakeResStorage_buffer = npareamaximum(self.var.lakeResStorage, self.var.waterBodyBuffer)

        if checkOption('calcWaterBalance'):
            self.model.waterbalance_module.waterBalanceCheck(
                [inflowCorrC],  # In
                [outflowC, EvapWaterBodyC],  # Out  EvapWaterBodyC
                [prereslake],  # prev storage
                [self.var.lakeResStorageC],
                "lake1", False)

        if checkOption('calcWaterBalance'):
            self.model.waterbalance_module.waterBalanceCheck(
                [self.var.sumlakeResInflow],  # In
                [self.var.sumlakeResOutflow, self.var.sumEvapWaterBodyC],  # Out  self.var.evapWaterBodyC
                [np.compress(self.var.compress_LR, self.var.prelakeResStorage)],  # prev storage
                [self.var.lakeStorageC],
                "lake2", False)

        if checkOption('calcWaterBalance'):
            self.model.waterbalance_module.waterBalanceCheck(
                [self.var.lakeResInflowM],  # In
                [self.var.lakeResOutflowM, self.var.EvapWaterBodyMOutlet],  # Out  self.var.evapWaterBodyC
                [self.var.prelakeResStorage / self.var.cellArea],  # prev storage
                [self.var.lakeResStorage / self.var.cellArea],
                "lake3", False)

        # report(decompress(runoff_LR), "C:\work\output3/run.map")

        return outLdd, lakeResOutflowDis

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
