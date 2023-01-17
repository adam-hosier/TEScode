#GUI for use with TES during data collection.
#Author: Grant Mondeel (gmondee) -- gmondee@g.clemson.edu
#Research advisor: Endre Takacs -- etakacs@clemson.edu

from click import progressbar
import numpy as np
import pylab as plt
import tkinter as tk
#from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import FileDialog, askdirectory         #from guiScripts
import mass
from mass.off import ChannelGroup, getOffFileListFromOneFile, Channel, labelPeak, labelPeaks
import os
import sys
import progress
import io
from matplotlib.lines import Line2D
import matplotlib
#from mass.calibration import _highly_charged_ion_lines         #this is added conditionally later on


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        #initialize window
        self.geometry("1500x1500")
        #add main title
        main_title=tk.Label(self, text="TES GUI\t\tVersion 0.1\tBy: Grant Mondeel    gmondee@g.clemson.edu")
        main_title.grid(row=1, column=0, columnspan=2)

        self.quitButton=tk.Button(self, text="QUIT", command=self.quitApp)     #FullCalStatesClick() displays the energy plots with Mass
        self.quitButton.grid(row=1, column=3)

        self.fullCalTF=False    #this is set to True after the user has added a calibration line. See TButtonClick() function for usage.

        #add the terminal output box at the bottom. this displays most of what is sent to stdout and stderr, but not the progress bars and some other Mass outputs (can't figure out why).
        self.errOutLabel = tk.Label(self,text='Output:')
        self.errOutLabel.grid(row=25, column=0, sticky='w')
        self.errOut = ConsoleText(master=self)
        self.errOut.grid(row=26,column=0, columnspan=5)
        self.errOut.start()
        self.load()
        ###prepare to choose file for import.  fileIndicator is red if no file is loaded, green if file is loaded
        #creates label to display folder and a button to send it to the calibration.
           
        #self.folder_path.set("Folder path here...")
        fileImportButton = tk.Button(self, text="Click to choose data folder", command=self.openFile)  #openFile() prompts user to select folder. 
        fileImportButton.grid(row=2, column=0)
        self.fileImportReadout = tk.Label(self, text=self.folder_path.get())
        self.fileImportReadout.grid(row=2, column=1)

        fileSubmitButton = tk.Button(self, text="SUBMIT", command=self.sendFile)    #sendFile() saves the selected folder into the variable self.d
        fileSubmitButton.grid(row=2, column=3)
        self.fileIndicator = tk.Label(self, text="", bg="red", width=10)    #red if nothing has been saved, green if something has been saved. not necessarily a functional input.
        self.fileIndicator.grid(row=2, column=4, sticky="e")


        ###Gets date, run, calibration state(s), calibration channel, and bin size. If calibration channel is left as 'None', Mass will select the first good channel.
        
        self.Date=guiField(self, labelText="Enter date as YYYYMMDD",myRow=3, defVal=self.dateDef)
        self.Run=guiField(self, labelText="Enter run as ####",myRow=4, defVal=self.runDef)
        self.CalState=guiField(self, labelText="Enter Calibration State(s) as A B C ...",myRow=5, defVal=self.stateDef)
        self.CalChannel=guiField(self, labelText="Enter Calibration Channel or use None to get First Good Channel",myRow=6, defVal=None)   
        self.BinSize=guiField(self, labelText="Enter Bin Size (Default: 0.7)",myRow=7, defVal="0.7")

        fileOpenButton = tk.Button(self, text="Click to view calibration channel", command=self.openRunFile)    #self.openRunFile starts the filtVal calibration. The next steps are called in here.
        fileOpenButton.grid(row=8, column=0)
        submitAllButton = tk.Button(self, text="SUBMIT ALL", command=self.submitAllClicked) #submitAllClicked functionally clicks all of the above buttons
        submitAllButton.grid(row=8, column=2)

        self.HCionBox() #function to initiate the tickbox to opt in/out of highly charged ions for the line selection step





    #FUNCTIONS
    def quitApp(self):  #quits out of the gui without crashing ipython
        sys.stdout = self.errOut.original_stdout
        sys.stderr = self.errOut.original_stderr
        self.destroy()

    def openFile(self): #opens the user's data folder (not any specific file). This should be where the date folders are ex. 20220310
        self.filename = tk.StringVar()
        self.filename.set(askdirectory()) # show an "Open" dialog box and return the path to the selected file
        self.folder_path.set(self.filename.get())
        print(self.folder_path.get())
        self.fileImportReadout.grid_forget()
        self.fileImportReadout = tk.Label(self, text=self.folder_path.get())
        self.fileImportReadout.grid(row=2, column=1)

    def sendFile(self): #sets file path variable "d". used for Mass.
        try:
            self.d= (self.folder_path.get())          #.replace("/","\\")      #windows-specific hack (not needed anymore but leaving it here)
        except:
            print("unable to set file")
            self.fileIndicator.configure(bg="red")
            #self.file.setIndicator(newColor="red")
        else:
            self.fileIndicator.configure(bg="green")
            #self.file.setIndicator(newColor="green")
            #print("file set successfully")

    def openRunFile(self):  #takes the user inputs, sets variables for Mass and opens the run files.
        #self.d set in sendFile
        self.date   =   str(self.Date.Val)                          #date
        self.rn     =   str(self.Run.Val)                           #run number, can be found in dastard command 
        self.state1 =   str(self.CalState.Val).split()              #states to view during calibration. can be just one 'A' or many 'A B C D ...'
        self.fl     =   getOffFileListFromOneFile(os.path.join(self.d, f"{self.date}", f"{self.rn}", f"{self.date}_run{self.rn}_chan1.off"), maxChans=192)  #Mass' way of finding all .off files in the folder.
        self.data   =   ChannelGroup(self.fl)                       #self.data is used when we do something to all of the channels. see self.ds below for individual channels.
        self.data.setDefaultBinsize(float(self.BinSize.Val))        #this bin size is used to display the calibrated plots.
        if self.CalChannel.Val == "None":    #selection of individual channel. If left as "None", Mass will select the first good channel
            self.ds =   self.data.firstGoodChannel()
        else:
            self.ds =   self.data[int(self.CalChannel.Val)]       
        self.save()                
        self.viewCalPlot()

    def save(self): #saves most recent entries to be loaded for the next time the program is started
        states=" ".join(self.state1)
        np.savetxt(fname='TESGUI_save.txt',X=np.asarray([str(self.d), self.date, self.rn, states]), fmt='%s')#np.column_stack([self.d, self.date, self.rn, self.state1]))  #saves directory, date, run#, calibration state(s)

    def load(self): #loads 'TESGUI_save.txt' into the user interface
        try:    #if the file exists
            d, date, run, calstate = np.loadtxt(fname='TESGUI_save.txt',dtype=str, delimiter='\n')

        except: #if the file doesn't exist, keep using (hard-coded) defaults
            print("Failed to load saved presets. Presets are saved by viewing calibration channel.")
            date, run, calstate = ['20220310','0002','A']
            d=("Folder path here...")

        self.folder_path = tk.StringVar()   
        self.folder_path.set(d)
        self.dateDef=date
        self.runDef=run
        self.stateDef=calstate




    def HCionBox(self):     #checkbox to toggle whether or not highly charged ions are available for selection during manual calibration
        self.HCvar=tk.IntVar()
        HCionCheckBox=tk.Checkbutton(self, text='View HC Ions?',variable=self.HCvar, onvalue=1, offvalue=0)
        HCionCheckBox.grid(row=8, column=3)

    def viewCalPlot(self):  #displays the filter value plot of the calibration channel on the inputted states and feeds mouse clicks back for calibration.
        try:
            def mouse_event(event): #gets x-position of mouse clicks on the filtVal plot
                self.calXValue = event.xdata
                self.showCalValue()
            self.myFig = plt.figure()
            self.myAxis=plt.gca()
            self.calXValue=0    #variable for most recently-clicked x-position
            self.ds.plotHist( np.arange(0, 60000, 20), "filtValue", coAddStates=False, states=self.state1, axis=self.myAxis )
            self.myFig.canvas.mpl_connect('button_press_event', mouse_event)
            self.initCalDC()    #initializes calibration plan and does drift correction of cal channel  
            self.showCalValue() #shows the x-position (filter value) that will be submitted for the calibration line
            self.showTOptions() #transition options dropdown
            self.myFig.show()
        except:
            print("Failed to open Calibration Plot. Try a different channel?")  #some channels are bad

    def submitAllClicked(self): #all initial "SUBMIT" buttons at once, for convenience
        self.Date.buttonClick()
        self.Run.buttonClick()
        self.CalState.buttonClick()
        self.CalChannel.buttonClick()
        self.BinSize.buttonClick()
        self.sendFile()

    def showCalValue(self): #displays the most recently clicked calibration value after the channel has been plotted
        #Display currently selected calibration x-value
        try:
            self.calXDisplay = tk.Label(self, text=np.round(self.calXValue,1))
            self.calXDisplay.grid_forget()
            self.calXDisplay = tk.Label(self, text=np.round(self.calXValue,1))
            self.calXDisplay.grid(row=9, column=0)
        except TypeError:   #throws type error if you click on the edge of the plot, where plt returns None instead of an x value.
            pass

    def showTOptions(self): #a dropdown with the transitions in the mass spectra library
        if self.HCvar.get() == 1:       #optional import of highly charged ions to the dropdown
            from mass.calibration import _highly_charged_ion_lines
        self.TOptionsDict=list(mass.spectra.keys())     #mass.spectra.keys() is a list of transitions stored in Mass
        self.TOption=tk.StringVar()
        self.TOption.set("MgKAlpha")    #default value for selected element
        self.TOptionsDropdown=tk.OptionMenu(self, self.TOption, *self.TOptionsDict)
        self.TOptionsDropdown.grid(row = 9, column = 1, sticky="w")
        self.TButton=tk.Button(self, text="SEND TO CALIBRATION", command=self.TButtonClick)
        self.TButton.grid(row=9, column=2)
        self.linesList = []     #linesList array is used later to show the line fits for each calibration transition

    def TButtonClick(self): #Sends the selected transition and filter value to calibration plan. Don't change channel after assigning one; click "view calibration channel" again to reset this
        self.ds.calibrationPlanAddPoint(int(self.calXValue), str(self.TOption.get()), states=self.state1)
        print("added", self.TOption.get(), "at", self.calXValue)
        self.linesList.append(self.TOption.get())
        #self.ds.plotHist( np.arange(0, 60000, 20), "filtValueDC", coAddStates=False, states=None )   #states=None by default uses all states
        if self.fullCalTF == False:     #creates the "START CALIBRATION" button after there's at least one calibration line
            self.fullCalTF==True
            self.startFullCalibration() 
    
    def initCalDC(self): #starts calibration of one channel with drift correction and makes the calibration plan
        self.ds.learnDriftCorrection(overwriteRecipe=True, states=self.state1)
        self.ds.calibrationPlanInit("filtValueDC")
        print("Finished drift correction.")

    def startFullCalibration(self): #makes button to initiate calibration of all channels. fullCalClick() actually does the calibrating steps and starts further steps.
        self.fullCalButton=tk.Button(self, text="START CALIBRATION", command=self.fullCalClick)
        self.fullCalButton.grid(row=10, column=0)

    def fullCalClick(self): #do the calibration for all channels. takes a while.
        plt.close(self.myFig)   #closes clickable calibration plot
        
        ###Single channel calibration
        # self.ds.learnPhaseCorrection(uncorrectedName="filtValueDC", linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
        # print("Starting calibration.")
        # self.ds.calibrateFollowingPlan("filtValueDCPC") 
        # self.ds.diagnoseCalibration()

        ###Full channel calibration. prints don't happen until the end, not sure why. this is a basic Mass calibration that could be improved if needed, starting with data.learnTimeDriftCorrection()
        print("Starting Drift Correction.")
        self.data.learnDriftCorrection(overwriteRecipe=True, states=self.state1)
        print("Aligning Channels.")
        self.data.alignToReferenceChannel(self.ds, "filtValueDC", np.arange(0,40000,30))
        print("Starting Phase Correction.")
        self.data.learnPhaseCorrection(uncorrectedName="filtValueDC", linePositionsFunc = lambda ds: ds.recipes["energyRough"].f._ph)
        print("Starting Calibration.")
        self.data.calibrateFollowingPlan("filtValueDCPC")
        print("Finished Calibration.")
        self.lineFitsChecker()  #button to generate data.linefit calls for each calibration transition

        ###GUI elements for displaying the calibrated energy plots.
        self.CoAdd=tk.IntVar()  #True or False; toggles if states are added in an energy plot
        self.FullCalCoAdd=tk.Checkbutton(self, text='Add states together?',variable=self.CoAdd, onvalue=True, offvalue=False)
        self.FullCalCoAdd.grid(row=12, column=1)
        self.FullCalStates=guiField(self, labelText="Enter States to Display",myRow=11, defVal=self.CalState.Val)
        self.FullCalStatesButton=tk.Button(self, text="PLOT THESE STATES", command=self.FullCalStatesClick)     #FullCalStatesClick() displays the energy plots with Mass
        self.FullCalStatesButton.grid(row=12, column=2)

        ###Button to start real-time plotting routine.
        self.updateIndex=0  #shows how many times plot has been updated
        self.StartRTP=tk.Button(self, text="START REAL-TIME PLOTTING",  command=self.initRTP)   
        self.StartRTP.grid(row=13, column=0)
        #Box to enter RTP update delay, in seconds
        self.RTPdelay=guiField(self, labelText="Enter RTP Update Delay (s):", myRow=14, defVal=15)
        #Button to stop real-time plotting routine.
        self.StopRTP=tk.Button(self, text="STOP RPT", command=self.stopRTP)     #FullCalStatesClick() displays the energy plots with Mass
        self.StopRTP.grid(row=13, column=1)

        ###Button to start 2D RTP routine.
        self.Start2dRTP=tk.Button(self, text="START 2D REAL-TIME PLOTTING",  command=self.init2dRTP)   
        self.Start2dRTP.grid(row=16, column=0)
        #Button to stop real-time plotting routine.
        self.Stop2dRTP=tk.Button(self, text="STOP 2D RPT", command=self.stop2dRTP)     #FullCalStatesClick() displays the energy plots with Mass
        self.Stop2dRTP.grid(row=16, column=1)
        #Input for DT timings
        self.DTcook=tk.Entry(self, text='Initial Cook (ms)')
        self.DTcook.grid(row=17,column=0)
        self.DTvolt=tk.Entry(self, text='DT voltage points (kV)')
        self.DTvolt.grid(row=17,column=1)
        self.DTtime=tk.Entry(self, text='DT time points (ms)')
        self.DTtime.grid(row=18,column=2)



    def initRTP(self): #creates axes, clears variables, and starts real-time plotting routine
        plt.ion()
        self.updateIndex=0                  #tracks how many updates have happened
        self.alphas=[]                      #transparencies of lines
        self.rtpline=[]                     #list of all plotted lines
        self.plottedStates=[]               #used for the RTP legend
        self.energyPlot=plt.figure()        #everything is plotted onto this figure
        self.energyAxis = plt.gca()
        plt.grid()
        plt.title('Real-time energy')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Counts per'+str(self.BinSize.Val)+'eV bin')
        self.rate = 0.1                     #how much to lower the transparency of lines each update
        self.RTPdelay.buttonClick()         #this does make the submit button redundant, but it makes the default value usable w/o submitting
        self.updateFreq = int(self.RTPdelay.Val)*1000             #in ms after the multiplication
        self.continueRTP = True             #toggle to keep updating plots
        self.UpdatePlots()

    def stopRTP(self):
        self.continueRTP = False

    def UpdatePlots(self):  #real-time plotting routine. refreshes data, adjust alphas, and replots graphs
        print(f"iteration {self.updateIndex}")
        self.data.refreshFromFiles()                #Mass function to refresh .off files as they are updated
        States = str.split(self.FullCalStates.Val)  #Gets last-submitted states from the lower states text box
        self.rtpdata=[]                             #temporary storage of data points
        self.rtpline.append([])                     #stores every line that is plotted according to the updateIndex

        for s in range(len(States)):    #looping over each state passed in
            self.rtpdata.append(self.data.hist(np.arange(0, 14000, float(self.BinSize.Val)), "energy", states=States[s]))   #[x,y] points of current energy spectrum, one state at a time
            self.energyAxis.plot(self.rtpdata[s][0],self.rtpdata[s][1],alpha=1, color=self.getColorfromState(States[s]))    #plots the [x,y] points and assigns color based on the state
            if States[s] not in self.plottedStates:     #new states are added to the legend; old ones are already there                      
                self.plottedStates.append(States[s])
            self.rtpline[self.updateIndex].append(self.energyAxis.lines[-1])    #stores most recent line

        self.alphas.append(1)   #lines in the same refresh cycle share an alpha (transparency) value. a new one is made for the newest lines.
        customLegend=[]         #temporary list to store the legend
        for s in self.plottedStates:    #loops over all states called during the current real-time plotting routine
            customLegend.append(Line2D([0], [0], color=self.getColorfromState(s)))      #each state is added to the legend with the state's respective color
        self.energyAxis.legend(customLegend, list(self.plottedStates))                  #makes the legend

        ###change transparency of current elements, plot adjusted lines
        for lineI in range(len(self.alphas)):   #loops over the number of refresh cycles, which is also the length of the alphas list
            if self.alphas[lineI] > 0.1:        #alpha values cannot be below 0
                self.alphas[lineI] = self.alphas[lineI] - self.rate
            for setI in range(len(self.rtpline[lineI])):    #loops over the states within one refresh cycle
                self.rtpline[lineI][setI].set_alpha(self.alphas[lineI])     #sets adjusted alpha values
        plt.draw()
        self.updateIndex=self.updateIndex+1
        if self.continueRTP==True:  #if off button hasn't been pressed
            self.after(self.updateFreq, self.UpdatePlots)   #calls the UpdatePlots function again after time has passed

    def init2dRTP(self): #creates axes, clears variables, and starts real-time plotting routine
        plt.ion()
        self.updateIndex=0                  #tracks how many updates have happened
        self.alphas=[]                      #transparencies of lines
        self.rtpline=[]                     #list of all plotted lines
        self.plottedStates=[]               #used for the RTP legend
        self.energyPlot=plt.figure()        #everything is plotted onto this figure
        self.energyAxis = plt.gca()
        plt.grid()
        plt.title('Real-time energy')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Counts per'+str(self.BinSize.Val)+'eV bin')
        self.rate = 0.1                     #how much to lower the transparency of lines each update
        self.RTPdelay.buttonClick()         #this does make the submit button redundant, but it makes the default value usable w/o submitting
        self.updateFreq = int(self.RTPdelay.Val)*1000             #in ms after the multiplication
        self.continueRTP = True             #toggle to keep updating plots
        self.UpdatePlots()

    def stop2dRTP(self):
        self.continueRTP = False

    def Update2dPlots(self):  #real-time plotting routine. refreshes data, adjust alphas, and replots graphs
        print(f"iteration {self.updateIndex}")
        self.data.refreshFromFiles()                #Mass function to refresh .off files as they are updated
        States = str.split(self.FullCalStates.Val)  #Gets last-submitted states from the lower states text box
        self.rtpdata=[]                             #temporary storage of data points
        self.rtpline.append([])                     #stores every line that is plotted according to the updateIndex

        for s in range(len(States)):    #looping over each state passed in
            self.rtpdata.append(self.data.hist(np.arange(0, 14000, float(self.BinSize.Val)), "energy", states=States[s]))   #[x,y] points of current energy spectrum, one state at a time
            self.energyAxis.plot(self.rtpdata[s][0],self.rtpdata[s][1],alpha=1, color=self.getColorfromState(States[s]))    #plots the [x,y] points and assigns color based on the state
            if States[s] not in self.plottedStates:     #new states are added to the legend; old ones are already there                      
                self.plottedStates.append(States[s])
            self.rtpline[self.updateIndex].append(self.energyAxis.lines[-1])    #stores most recent line

        self.alphas.append(1)   #lines in the same refresh cycle share an alpha (transparency) value. a new one is made for the newest lines.
        customLegend=[]         #temporary list to store the legend
        for s in self.plottedStates:    #loops over all states called during the current real-time plotting routine
            customLegend.append(Line2D([0], [0], color=self.getColorfromState(s)))      #each state is added to the legend with the state's respective color
        self.energyAxis.legend(customLegend, list(self.plottedStates))                  #makes the legend

        ###change transparency of current elements, plot adjusted lines
        for lineI in range(len(self.alphas)):   #loops over the number of refresh cycles, which is also the length of the alphas list
            if self.alphas[lineI] > 0.1:        #alpha values cannot be below 0
                self.alphas[lineI] = self.alphas[lineI] - self.rate
            for setI in range(len(self.rtpline[lineI])):    #loops over the states within one refresh cycle
                self.rtpline[lineI][setI].set_alpha(self.alphas[lineI])     #sets adjusted alpha values
        plt.draw()
        self.updateIndex=self.updateIndex+1
        if self.continueRTP==True:  #if off button hasn't been pressed
            self.after(self.updateFreq, self.UpdatePlots)   #calls the UpdatePlots function again after time has passed

    def getColorfromState(self, s): #pass in a state label string like 'A' or 'AC' (for states past Z) and get a color index using plt.colormaps() 
        c = plt.cm.get_cmap('gist_rainbow')     #can be other colormaps, rainbow is most distinct
        maxColors=8                             #how many unique colors there are. lower values make it easier to distinguish between neighboring states.
        cinter=np.linspace(0,1,maxColors)       #colormaps use values between 0 and 1. this interpolation lets us assign integers to colors easily.
        cIndex=0
        cIndex=cIndex+ord(s[-1])-ord('A')       #uses unicode values of state labels (e.g. 'A') to get an integer 
        if len(s)!=1:   #for states like AA, AB, etc., 27 is added to the value of the ones place
            cIndex=cIndex+26+ord(s[0])-ord('A')+1 #26 + value of first letter (making A=1 so AA != Z)

        while cIndex >= maxColors:       #colors repeat after maxColors states. loops until there is a valid interpolated index from cinter
            cIndex=cIndex-maxColors
        #print(c(cinter[cIndex]))
        return c(cinter[cIndex])        #returns values that can be assigned as a plt color, looks like (X,Y,Z,...) or something similar       

    def lineFitsChecker(self):  #shows button to "view line fits"
        lineFitsButton=tk.Button(self, text="VIEW LINE FITS", command=self.lineFitsClick)
        lineFitsButton.grid(row=10, column=1)

    def lineFitsClick(self):    #generates Mass linefit plots for each line added in the calibration step
        for LineName in self.linesList: #loops over each line
            self.data.linefit(LineName, states=str.split(self.FullCalStates.Val), dlo=20, dhi=20)   #calls Mass linefit. dlo/dhi are how far from the centroid to fit
        plt.show(block=False)

    def FullCalStatesClick(self, eAxis=None):   #generates Mass calibrated energy plots for inputted states. 
        ###self.CoAdd is True or False and determines whether plots are added into one curve or displayed individually. eAxis ("energy axis") can be used to plot multiple histograms on one figure
        States=str.split(self.FullCalStates.Val)
        self.data.plotHist( np.arange(0,14000,1), "energy", coAddStates=self.CoAdd.get(), states=States, axis=eAxis)    
        plt.show(block=False)

        


class guiField(tk.Frame):   #makes a row with elements: label, entry, button(for submitting) and indicator on row myRow.
    def __init__(self, root, labelText, showLabel=True,showEntry=True, showButton=True, showIndicator=True, myRow=0, defVal=None):
        super().__init__()
        self.root=root      #tells tkinter to use the same window as the other elements
        self.Row=int(myRow) #used to keep all elements on the same row
        self.Column=0       #index for the first unused column. this should be incremented after each element is shown
        self.defVal=defVal  #default value for text boxes, shows up when it is created. ex: defVal for calibration channel is likely 1.
        if showLabel==True:
            self.labelText=labelText
            self.makeLabel()

        if showEntry==True:
            self.makeEntry()

        if showButton==True:
            self.makeButton()

        if showIndicator==True:
            self.makeIndicator()

    def makeLabel(self):        #creates a label with tkinter that displays the self.labelText passed into the object.
        self.label=tk.Label(self.root, text=self.labelText)
        self.label.grid(row=self.Row, column=self.Column)
        self.Column = self.Column + 1

    def makeEntry(self):        #creates a text box with tkinter that the user can type in. the value is read with the buttonClick function or with self.entry.get()
        self.entry=tk.Entry(self.root)  #this is what the user types in and where the self.Val input comes from. THIS IS A STRING BY DEFAULT!!!
        self.entry.insert(0, str(self.defVal))
        self.entry.grid(row=self.Row, column=self.Column)
        self.Column = self.Column + 1

    def makeButton(self):       #creates a button with tkinter that is linked to the buttonClick function
        self.button=tk.Button(self.root, text="SUBMIT", command=self.buttonClick)
        self.button.grid(row=self.Row, column=self.Column)
        self.Column = self.Column + 1

    def buttonClick(self):      #When the guiField's button is clicked, it loads the entry value into self.Val
        try:
            self.Val=self.entry.get()   #reads the entry's current contents as a string. you will often need to typecast to an int/float. e.g. int(self.Date.Val)
        except:
            print("unable to submit")
            self.setIndicator("red")
        else:
            self.setIndicator("green")

    def makeIndicator(self):    #Makes an indicator to show submission status. Looks like a red/green rectangle on the right.
        self.indicator = tk.Label(self.root, text="", bg="red", width=10)
        self.indicator.grid(row=self.Row, column=self.Column, sticky="e")

    def setIndicator(self, newColor):  #sets the indicator of one of the inputs to some color (red or green)
        self.indicator.configure(bg=newColor)


###These classes are copied from https://stackoverflow.com/questions/18517084/how-to-redirect-stdout-to-a-tkinter-text-widget
###they help display the stdout/stderr messages to the GUI
class ConsoleText(tk.Text):
    '''A Tkinter Text widget that provides a scrolling display of console
    stderr and stdout.'''

    class IORedirector(object):
        '''A general class for redirecting I/O to this Text widget.'''
        def __init__(self,text_area):
            self.text_area = text_area

    class StdoutRedirector(IORedirector):
        '''A class for redirecting stdout to this Text widget.'''
        def write(self,str):
            self.text_area.write('\r'+str,False)

        def flush(*args, **kwargs):
            pass

    class StderrRedirector(IORedirector):
        '''A class for redirecting stderr to this Text widget.'''
        def write(self,str):
            self.text_area.write(str,True)

        def flush(*args, **kwargs):
            pass


    def __init__(self, master=None, cnf={}, **kw):
        '''See the __init__ for Tkinter.Text for most of this stuff.'''
        tk.Text.__init__(self, master, cnf, **kw)

        self.started = False

        self.tag_configure('STDOUT',background='white',foreground='black')
        self.tag_configure('STDERR',background='white',foreground='red')


    def start(self):
        
        if self.started:
            return

        self.started = True

        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr

        #self.pOrig_stdout=progress.stderr

        stdout_redirector = ConsoleText.StdoutRedirector(self)
        stderr_redirector = ConsoleText.StderrRedirector(self)

        sys.stdout = stdout_redirector
        sys.stderr = stderr_redirector
        # progress.Infinite.check_tty = False   
        #progress.Infinite.file = sys.stdout    #this redirects the progress bars to the "Output" box, but it doesn't delete old bars and doesnt go to the next line. disabled for now.


    def stop(self):

        if not self.started:
            return

        self.started = False

        sys.stdout = self.original_stdout
        sys.stderr = self.original_stderr

    def write(self,val,is_stderr=False):
        #self.write_lock.acquire()
        self.insert('end',val,'STDERR' if is_stderr else 'STDOUT')
        self.see('end')
        #self.write_lock.release()


if __name__=="__main__":
    plt.ion()
    matplotlib.use('TkAgg')     #changes backend so matplotlib and tkinter work at the same time. will display a frozen black window for plt plots otherwise.
    app=App()                   #starts the GUI
    app.mainloop()              #mainloop keeps the GUI running until the window is closed
    












    



