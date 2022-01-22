import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xlrd
import xlwt
import os



class load_CMG_simulation():
    def __init__(self, file_folder, file_name):
        self.folder = file_folder
        self.file = file_name
        

        container = STARScontainer()
        #container.parse_runfile(str(self.folder)+'/'+str(self.file))
        container.parse_runfile(os.path.join(self.folder, self.file))
        self.time = np.insert(container.t, 0, 0) #[hr]
        df_sim = pd.DataFrame({})
        df_sim['Time'] = self.time
        self.TEMP = container.TEMP_output
        self.SPEC_VALs = container.SPEC_VALS
        self.TEMP1 = container.TEMP1
        self.TEMP2 = container.TEMP2
        self.TEMP3 = container.TEMP3
        self.TEMP4 = container.TEMP4
        self.TEMP5 = container.TEMP5

        self.dfsim = df_sim

    def print_info(self):
        print('The following simulation file is read: '+str(self.folder)+'/'+str(self.file))





# Kinetic cell matching code from Tim
class RTO:
    def __init__(self, NAME, HR, TIME, TEMP, O2, X_O2, GRAD_X_O2):
        self.NAME = NAME
        self.HR = HR
        self.TIME = TIME
        self.TEMP = TEMP
        self.O2 = O2
        self.X_O2 = X_O2
        self.GRAD_X_O2 = GRAD_X_O2        
        
    def getNAME(self):
        return self.NAME
    def getHR(self):
        return self.HR
    def getTIME(self):
        return self.TIME
    def getTEMP(self):
        return self.TEMP
    def getO2(self):
        return self.O2
    def getX_O2(self):
        return self.X_O2
    def getGRAD_X_O2(self):
        return self.GRAD_X_O2    
    
    def setNAME(self, value):
        self.NAME = value
    def setHR(self, value):
        self.HR = value
    def setTIME(self, value):
        self.TIME = value
    def setTEMP(self, value):
        self.TEMP = value
    def setO2(self, value):
        self.O2 = value
    def setX_O2(self, value):
        self.X_O2 = value
    def setGRAD_X_O2(self, value):
        self.GRAD_X_O2 = value      
        
"""
Change the numSPHIST and self.SPEC_VALS according to dat files
"""
class STARScontainer():
    def __init__(self):
        self.numSPHIST = 5
        
    def parse_runfile(self, input_file_base, CMG_type = 'JABS'):

        # Open STARS output files
        FID = [line.decode('utf-8', errors='ignore') for line in open(input_file_base + '.irf', 'rb')]
        FID_bin = open(input_file_base + '.mrf', 'rb')

        # Initialize variables/objects for parsing data
        i = 0
        self.UNIT = {} # Create units dictionary
        self.COMP = {} # Compositions dictionary
        self.TIME = {} # Time dictionary
        self.GRID = {} # Grid dictionary
        self.GRID['TEMP'] = []
        self.SP = {}
        self.SP['VAL'] = []
        self.SPHIST = {}

        self.TIME['VECT'] = []
        self.TIME['CHR'] = []
        self.TIME['CHR_UNIT'] = []

        self.REC = {}
        REC_list = ['WELL-REC', 'LAYER-REC', 'GROUP-REC', 'SECTOR-REC', 'RSTSPEC01-REC', 'RSTSPEC02-REC', 'RSTSPEC03-REC',
            'RSTSPEC04-REC', 'RSTSPEC05-REC', 'RSTSPEC06-REC', 'RSTSPEC07-REC', 'RSTSPEC08-REC', 'RSTSPEC09-REC',
            'RSTSPEC10-REC', 'RSTSPEC11-REC', 'RSTSPEC12-REC', 'RSTSPEC13-REC', 'RSTSPEC14-REC', 'RSTSPEC15-REC',
            'RSTSPEC16-REC', 'RSTSPEC17-REC', 'RSTSPEC18-REC', 'RSTSPEC19-REC', 'RSTSPEC20-REC', 'RSTSPEC21-REC',
            'RSTSPEC22-REC']

        skip_list = ['WELL-ARRAY', 'LAYER-ARRAY', 'GROUP-ARRAY']

        # Iterate over lines of file
        while i < len(FID):
            line_split = FID[i].split()

            # Parse case-by-case quantities
            if line_split[0] == 'INTERNAL-UNIT-TABLE':
                i+=1 
                self.UNIT['INT'] = [str(FID[i+j].split()[1]) for j in range(21)]
                self.UNIT['DESCP'] = [str(FID[i+j].split()[3]) for j in range(21)]
                i+=21 

            elif line_split[0] == 'OUTPUT-UNIT-TABLE':
                i+=1
                self.UNIT['OUT'] = [str(FID[i+j].split()[1]) for j in range(21)]
                i+=21

            elif line_split[0] == 'NCOMP':
                self.COMP['NUM'] = int(line_split[1])
                i+=1

            elif line_split[0] == 'COMPNAME':
                i+=1
                self.COMP['NAME'] = [str(FID[i+j].split()[1]).replace("'","") for j in range(self.COMP['NUM'])]
                i+=self.COMP['NUM']

            elif line_split[0] == 'COMP-PHASE-TEMPLATE':
                i+=1
                self.COMP['TEMPL'] = [str(FID[i+j].split()[2]) for j in range(self.COMP['NUM'])]
                i+=self.COMP['NUM']

            elif line_split[0] == 'TIME' and 'SPEC-HISTORY' in [FID[i-1].split()[0], FID[i-3].split()[0]]:
                self.TIME['VECT'].append(line_split[1:])
                i+=1

            elif line_split[0] == 'TIMCHR':
                self.TIME['CHR'].append(line_split[2])
                self.TIME['CHR_UNIT'].append(line_split[3].replace("'",""))
                i+=1


            # Parse grid values from binary file
            elif line_split[0] == 'GRID':
                item_num = int(FID[i].split()[2])
                for j in range(item_num):
                    num_bytes = np.fromfile(FID_bin, np.int64, count=1).byteswap()
                    _ = np.fromfile(FID_bin, np.int8, count=np.asscalar(num_bytes)).byteswap() 
                _, i = self.parse_nobin(FID, i)

            elif line_split[0] == 'GRID-VALUE':
                props_list, i = self.parse_nobin(FID, i)
                for prop in props_list[3:]:
                    num_bytes = np.fromfile(FID_bin, np.int64, count=1).byteswap()
                    if prop == 'TEMP':
                        grid_temp = np.fromfile(FID_bin, np.float64, count=int(num_bytes/8)).byteswap()
                        self.GRID['TEMP'].append(grid_temp)
                    else:
                        try:
                            _ = np.fromfile(FID_bin, np.int8, count=np.asscalar(num_bytes)).byteswap()
                        except:
                            pass


            # Parse species names and values from binary file
            elif line_split[0] == 'SPHIST-NAMES':
                self.SPHIST['NUM'] = []
                self.SPHIST['NAME'] = []
                i+=1
                for j in range(self.numSPHIST):
                    line_split = FID[i+j].split()
                    self.SPHIST['NUM'].append(line_split[0])
                    self.SPHIST['NAME'].append(' '.join([line.replace("'","") for line in line_split[3:]]))
                i+=self.numSPHIST

            elif line_split[0] == 'SPEC-HISTORY':
                props_list, i = self.parse_nobin(FID, i)
                for prop in props_list[3:]:
                    num_bytes = np.fromfile(FID_bin, np.int64, count=1).byteswap()
                    # print(prop)
                    # print(np.asscalar(num_bytes))
                    if prop == 'SPVALS':
                        spvals_temp = np.fromfile(FID_bin, np.float64, count=self.numSPHIST).byteswap()
                        self.SP['VAL'].append(spvals_temp)
                    else:
                        _ = np.fromfile(FID_bin, np.byte, count=np.asscalar(num_bytes)).byteswap()


            # Variables stored in binary but skipped in parsing
            elif line_split[0] in skip_list:
                item_num = int(line_split[2])
                for j in range(item_num):
                    num_bytes = np.fromfile(FID_bin, np.int64, count=1).byteswap()
                    try:
                        _ = np.fromfile(FID_bin, np.byte, count=np.asscalar(num_bytes)).byteswap()
                    except:
                        pass
                i+=1


            # Other variables parsed from text file
            elif line_split[0] in REC_list:
                list_temp = FID[i].split()
                while list_temp[-1] != '/':
                    i+=1
                    list_temp += FID[i].split()
                self.REC[line_split[0]] = list_temp[1:-1]
                i+=1

            else:
                i+=1

        FID_bin.close()

        # OUTPUT VARIABLES
        self.t = (np.asarray(self.TIME['VECT'])[:,1]).astype(np.float)
        TEMP_VALS = np.asarray(self.GRID['TEMP'])
        self.TEMP_output = TEMP_VALS
        #self.lin_HR = np.asarray(TEMP_VALS)[:,8]
        self.SPEC_VALS = np.asarray(self.SP['VAL'])
        

        if CMG_type == 'JABS':
            self.TEMP1 = self.SPEC_VALS[:,0]
            self.TEMP2 = self.SPEC_VALS[:,1]
            self.TEMP3 = self.SPEC_VALS[:,2]
            self.TEMP4 = self.SPEC_VALS[:,3]
            self.TEMP5 = self.SPEC_VALS[:,4]
            
            
        else:
            print('This type of CMG result reader has not been implemented yet .....')

            
    @staticmethod
    def parse_nobin(file_in, i):
        '''
        file_in - file to be parsed
        i - current line number
        '''
        
        list_temp = file_in[i].split()
        while list_temp[-1] != '/':
            i+=1
            list_temp += file_in[i].split()
        i+=1

        return list_temp[1:-1], i
