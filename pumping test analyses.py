#####################################################################
#
# PumpingTest.py
#
# aquifer test analysis by various methods
#
# by Walt McNab
#
#####################################################################

from functools import partial
from numpy import *
from scipy.integrate import quad
from scipy.integrate import odeint
from scipy.special import *
import matplotlib.pyplot as plt
import sys
from PyQt5 import QtCore, QtWidgets, uic

# support for user interfacce
qtCreatorFile = 'pumping_test_interface.ui'
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


class DataSet:
    
    def __init__(self):
        # pumping test data (time and drawdown arrays)
        self.t = []
        self.s = []
        lineInput = []        
        inputFile = open('transducer.txt','r')
        for line in inputFile: lineInput.append(line.split())
        inputFile.close()
        for i in range(1, len(lineInput)):
            self.t.append(float(lineInput[i][0]))
            self.s.append(float(lineInput[i][1]))
        self.t = array(self.t)
        self.s = array(self.s)
        print('Read test data set.')


class Aquifer:
    
    def __init__(self):
         # aquifer characteristics
        lineInput = []        
        inputFile = open('aquifer.txt','r')
        for line in inputFile: lineInput.append(line.split())
        inputFile.close()
        self.K = float(lineInput[0][1])     # aquifer properties
        self.Ss = float(lineInput[1][1])
        self.Sy = float(lineInput[2][1])
        self.b = float(lineInput[3][1])     # used as saturated thickness for unconfined aquifer
        self.bc = float(lineInput[4][1])        
        self.Kc = float(lineInput[5][1])    # 'c' refers to clay/aquitard
        self.Ssc = float(lineInput[6][1])
        self.S = self.Ss * self.b           # derive storage coefficient from specific storage
        print('Read aquifer characteristics.')

    def WriteValues(self):
        # update parameter file with current values
        output_file = open('aquifer.txt','w')
        output_file.writelines(['K', '\t', str(self.K),'\n'])
        output_file.writelines(['Ss', '\t', str(self.Ss), '\n'])
        output_file.writelines(['Sy', '\t', str(self.Sy), '\n'])
        output_file.writelines(['b', '\t', str(self.b), '\n'])
        output_file.writelines(['bc', '\t', str(self.bc), '\n'])        
        output_file.writelines(['Kc', '\t', str(self.Kc), '\n'])        
        output_file.writelines(['Ssc', '\t', str(self.Ssc), '\n'])
        output_file.close()         
        
        
class Well:
    
    def __init__(self, t0, tEnd):
        # well properties
        lineInput = []        
        inputFile = open('well.txt','r')
        for line in inputFile: lineInput.append(line.split())
        inputFile.close()
        self.r = float(lineInput[0][1])     # well radius; assume radial distance for monitoring drawdown
        self.Q = float(lineInput[1][1])     # pumping rate from well (negative value = extraction)
        self.tArray = logspace(log10(t0), log10(tEnd), num=60, endpoint=True)     # evaluation times

    def WriteValues(self):
        # update parameter file with current values
        output_file = open('well.txt','w')
        output_file.writelines(['r', '\t', str(self.r),'\n'])
        output_file.writelines(['Q', '\t', str(self.Q),'\n'])
        output_file.close()     

        
class Hantush:            # Hantush and Jacob (1955) solution

    def __init__(self, aquifer, well):
        self.B = sqrt(aquifer.bc*aquifer.K*aquifer.b/aquifer.Kc)
        self.aquifer = aquifer
        self.well = well

    def Integrand(self, y):
        # integral term for the Hantush well function
        x = exp(-y - self.well.r**2/(4.*self.B**2*y))/y
        return x

    def W(self, u):
        # Hantush well function
        x = quad(self.Integrand, u, +inf)[0]
        return x
        
    def Drawdown(self):
        s = zeros(len(self.well.tArray), float)
        for i, t in enumerate(self.well.tArray):        
            u = self.well.r**2*self.aquifer.Ss/(4*self.aquifer.K*t)
            s[i] = -self.well.Q/(4*pi*self.aquifer.K*self.aquifer.b) * self.W(u)
        return s

        
class Theis:    # Theis (1935) solution

    def __init__(self, aquifer, well):
        self.aquifer = aquifer
        self.well = well
        
    def W(self, u):
        # Theis well function
        return expn(1, u)        

    def Drawdown(self, mode):
        s = zeros(len(self.well.tArray), float)
        if mode == 0:       # confined aquifer
            for i, t in enumerate(self.well.tArray):    
                u = self.well.r**2 * self.aquifer.Ss/(4*self.aquifer.K*t)
                s[i] = -self.well.Q/(4*pi*self.aquifer.K*self.aquifer.b) * self.W(u)
        else:               # unconfined aquifer (assuming ~ constant saturated thickness)
            for i, t in enumerate(self.well.tArray):    
                u = self.well.r**2 * self.aquifer.Sy/(4*self.aquifer.K*self.aquifer.b*t)
                s[i] = -self.well.Q/(4*pi*self.aquifer.K*self.aquifer.b) * self.W(u)                
        return s
       

class MOL:  # numerical (method-of-lines) solution for an unconfined aquifer
    
    def __init__(self, aquifer, well):
        self.aquifer = aquifer
        self.well = well
        self.N = 70                                                 # default number of radial grid cells       
        self.rFace = self.Gridder()                                 # array of grid cell interface radii
        self.r = 0.5*self.rFace[1:] + 0.5*self.rFace[:-1]           # radius of node point associated with each cell
        self.r = insert(self.r, 0, self.well.r)                     # cell representing well
        self.A = pi*(self.rFace[1:]**2 - self.rFace[:-1]**2)        # base areas associated with individual grid cells
        self.A = insert(self.A, 0, pi*self.rFace[0]**2)
        self.Sy = zeros(self.N, float) + aquifer.Sy                 # assign storage coefficient of 1.0 to wellbore cell
        self.Sy = insert(self.Sy, 0, 1.0)
        self.S = zeros(self.N, float) + aquifer.S
        self.S = insert(self.S, 0, 1.0)       
    
    def Gridder(self):
        # generate radial grid
        rb = self.aquifer.b * 100.                   # set fixed boundary condition = 10X the available drawdown        
        index = arange(0, self.N+1, 1)
        f = 10.**(log10((rb/self.well.r))/self.N)   # sequential scaling factor
        r = self.well.r * f**index
        return r

    def Dupuit(self, h, t):
        # ordinary differential equations (volumetric balance for water) for grid cells; variable saturated thickness
        J = 2. * pi * self.aquifer.K * self.rFace[:-1] * (0.5*h[1:] + 0.5*h[:-1]) * (h[1:] - h[:-1]) / (self.r[1:] - self.r[:-1])
        J = insert(J, 0, -self.well.Q)  
        J = append(J, 2.*pi*self.aquifer.K*self.rFace[-1]*(0.5*h[-1]+0.5*self.aquifer.b)
            *(self.aquifer.b-h[-1])/(self.rFace[-1]-self.r[-1]))            # append flux from across exterior boundary
        dhdt = (J[1:] - J[:-1]) / (self.A * self.Sy)
        return dhdt       
    
    def Theis(self, h, t):
        # ordinary differential equations (volumetric balance for water) for grid cells; fixed saturated thickness
        J = 2. * pi * self.aquifer.K * self.rFace[:-1] * self.aquifer.b * (h[1:] - h[:-1]) / (self.r[1:] - self.r[:-1])
        J = insert(J, 0, -self.well.Q)                                      # express pumping as extraction from well
        J = append(J, 2.*pi*self.aquifer.K*self.rFace[-1]*self.aquifer.b
            *(self.aquifer.b-h[-1])/(self.rFace[-1]-self.r[-1]))            # append flux from across exterior boundary
        dhdt = (J[1:] - J[:-1]) / (self.A * self.S)
        return dhdt 
    
    def Drawdown(self, mode):
        # solve the transient unconfined aquifer test problem using the numerical method-of-lines
        h = zeros(self.N+1,float) + self.aquifer.b
        if mode == 0: h_t = odeint(self.Dupuit, h, self.well.tArray)
        else: h_t = odeint(self.Theis, h, self.well.tArray)
        h_t = transpose(h_t)
        s = self.aquifer.b - h_t[0]         # drawdown vector for cell representing well bore
        return s      
        

class GUI(QtWidgets.QMainWindow, Ui_MainWindow):
    
    def __init__(self, aquifer, well, tests, data):
        
        # initiate GUI
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        
        # load values read from files
        self.KhInput.setText(str(aquifer.K))
        self.SsInput.setText(str(aquifer.Ss)) 
        self.SyInput.setText(str(aquifer.Sy)) 
        self.bInput.setText(str(aquifer.b)) 
        self.bcInput.setText(str(aquifer.bc)) 
        self.KcInput.setText(str(aquifer.Kc))  
        self.SscInput.setText(str(aquifer.Ssc)) 
        self.rInput.setText(str(well.r))
        self.QInput.setText(str(well.Q))        

        # button functionality
        self.pushUpdate.clicked.connect(partial(self.Update, aquifer, well))
        self.pushEval.clicked.connect(partial(self.Evaluate, aquifer, well, tests, data))       
        self.pushSave.clicked.connect(partial(self.SaveFiles, aquifer, well))
        self.pushExit.clicked.connect(self.close)
       
    def Update(self, aquifer, well):
        # update aquifer and well objects with values on form
        aquifer.K = float(self.KhInput.text())
        aquifer.Ss = float(self.SsInput.text())
        aquifer.Sy = float(self.SyInput.text())
        aquifer.b = float(self.bInput.text())
        aquifer.bc = float(self.bcInput.text())
        aquifer.Kc = float(self.KcInput.text()) 
        aquifer.Ssc = float(self.SscInput.text())
        aquifer.S = aquifer.Ss * aquifer.b      
        well.r = float(self.rInput.text())
        well.Q = float(self.QInput.text())
    
    def SaveFiles(self, aquifer, well):
        # write current model to aquifer and well files
        aquifer.WriteValues()
        well.WriteValues()
    
    def Evaluate(self, aquifer, well, tests, data):
    
        # unpack test objects
        theis = tests[0]
        hantush = tests[1]
        numericWaterTable = tests[2]    
        
        # plot transducer data  
        plt.scatter(data.t, data.s, s=10, facecolors='none', edgecolors='black', label = 'Data') 

        # run checked models and add to plot
        if self.checkTheisConf.checkState():
            sTheisC = theis.Drawdown(0)
            plt.plot(well.tArray, sTheisC, color = 'red', label = 'Confined (Theis)')
        if self.checkMOLTheis.checkState():
            sMOLt = numericWaterTable.Drawdown(1)
            plt.plot(well.tArray, sMOLt, color = 'magenta', label = 'Confined (wellbore storage)')             
        if self.checkHantush.checkState():
            sHantushS = hantush.Drawdown()
            plt.plot(well.tArray, sHantushS, color = 'green', label = 'Leaky (Hantush & Jacob)')              
        if self.checkTheisUnconf.checkState():
            sTheisU = theis.Drawdown(1)
            plt.plot(well.tArray, sTheisU, color = 'blue', label = 'Unconfined (Theis, with Sy)')            
        if self.checkMOLDupuit.checkState():
            sMOLd = numericWaterTable.Drawdown(0)
            plt.plot(well.tArray, sMOLd, color = 'cyan', label = 'Unconfined (Dupuit; numerical)')
        
        plt.xscale('log')
        plt.yscale('log')        
        plt.xlabel('Time')
        plt.ylabel('Drawdown')
        plt.legend(loc=4)
        plt.show()    
        

### main script ###

def PumpTest():

    # read parameters
    data = DataSet()
    well = Well(data.t.min(), data.t.max())
    aquifer = Aquifer()

    # set up test method objects
    theis = Theis(aquifer, well)
    hantush = Hantush(aquifer, well)
    numericWaterTable = MOL(aquifer, well)
    tests = [theis, hantush, numericWaterTable]
    
    # set up GUI
    app = QtCore.QCoreApplication.instance()
    if app is None: app = QtWidgets.QApplication(sys.argv)
    window = GUI(aquifer, well, tests, data)
    window.show()
    sys.exit(app.exec_())
    
    
# run script    
PumpTest()
    


    
