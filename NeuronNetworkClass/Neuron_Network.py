from neuron import Neuron
from Adjacency_Matrix import Adjacency_Matrix
import random
from scipy.stats import invgauss
import subprocess
import sys



class Neuron_Network:
    '''
    Important Notes on Adjacency Matrix Map: 
    
        Currently the convention is that interneurons are loaded first, then the excitatory neurons, and then a decoder neuron (if applicable). 
        This order must be retained so that connection probabilities are assigned correctly. 
    '''
    
    import os
    
    
    def __init__(self, duration, timeStep,resolution,networkID):
        '''
        # Initializes the network. It starts of as empty and takes in the variables:
        #   Duration
        #   Samplerate aka timestep
        #   The resolution of data for matlab to graph
            Network ID will be used for coordinating cross network communication
        '''
        self.listOfNeurons = []
        self.adjacencyMatrix = Adjacency_Matrix()
        self.duration = duration
        self.timeStep = timeStep
        self.arrayOfInitialSpikeTimes = []
        self.deadCopyOfArrayOfInitialSpikeTimes = []
        self.mapImported = False
        
        self.resolution = resolution
        self.fakeExciteAmplitude = 0
        self.fakeInhibAmplitude = 0
        self.firstWrite = True
        self.networkID = networkID  # store as a string. literally like 'CA3'
        self.lastNeuronToFire = 0
       
        
        # Data needed for run a single dt method (NOT THE SAME ONES IN THE RUN METHOD)
        self.numberOfStepsNeeded = 0
        self.actualStepsTaken = 0
        self.totalTime = 0
        self.stepsPer5Percent = 0
        self.percentComplete = 0           
        self.mailStation = []
        self.runFinished = False
        
        
        # Data Library List
        self.listOfVoltages = []
        
    
    def getNetworkID(self):
        '''
        Returns network ID
        '''
        
        return self.networkID
    
    
    def addNeuron(self, spikeAmplitude, spikeAmplitude2, exciteTimeConstant,inhibTimeConstant,designation,delay1,delay2,wantBGCurrent):
        '''
        # Adds a single neuron to the network
        # It also automatically adds a section for itself in the adjacency matrix
        '''
        iD = len(self.listOfNeurons)
        neuron = Neuron(spikeAmplitude, spikeAmplitude2, exciteTimeConstant,inhibTimeConstant,designation,iD,delay1,delay2,wantBGCurrent,self.networkID)
        neuron.setTime(self.timeStep,self.duration,self.resolution)
        self.listOfNeurons.append(neuron)
        self.adjacencyMatrix.addNeuron()
    
    def getNumberOfNeurons(self):
        return len(self.listOfNeurons)
    
    
    def initialSpikeGeneration(self, numberOfSpikesPerNeuron, meanTime, synchronny,typeOfDraw,a,b,exciteExciteAmplitude,exciteInhibAmplitude):
        '''
        Generates the intital excitation spikes aka fake neurons
        Also will write the spike times to a file. This is accomplished here
        instead of the write function because the master list of spike times
        gets destroyed in the run process (Note: as of the latest update, we
        are now keeping a dead copy of the initial spike time list, because
        it proves better to do so for the synchronny trials, however, because
        the current write function works perfectly, I have not changed it to
        use the dead copy).
        This method is now also responsible for the amplitudes of the fake
        neurons
        
        Note: Fake neurons are essentially the ones created in the initial spike generation
        they aren't really modelable neurons. The are quite literally a list of when to send currents
        to the 'real' neurons of the entire system
        '''
        self.fakeExciteAmplitude = exciteExciteAmplitude
        self.fakeInhibAmplitude = exciteInhibAmplitude
        
        # Cleans out pre-existing fake neuron data
        for file in self.os.listdir(self.networkID+'/FakeNeuronSpikes/'):
            self.os.remove("FakeNeuronSpikes/" + file)   
            
        # for type of draw:
        #     1 = Gaussian
        #     2 = Exponential 
        #     3 = Uniform
        #     4 = Inverse Gauss
        if meanTime > self.duration:
            print("\nError: Your mean time is greater than the duration, \nand you can ignore the next error statement")
            return
        self.arrayOfInitialSpikeTimes = []
        self.deadCopyOfArrayOfInitialSpikeTimes = []
        
        # Essentially the "for loop's" job is to iterate once for every single
        # neuron in the list of all Neurons
        # The "while loop" will loop the amount of times as spikes desired per neuron
        # for every "for" iteration, one temporary list will store the spike times
        # and the id number for the neuron to recieve those spikes. At then end 
        # of the single for iteration, the temporary list will be appended to the
        # master list of spike times.
        for x in range(len(self.listOfNeurons)-1):
            # Fist Write is used to help keep a csv format
            self.firstWrite = True
            excitationSpikesdataFile = open("FakeNeuronSpikes/excitationSpike" + str(x+1) + ".txt", "w")
            numberOfSpikesGenerated = 0
            individualNeuronSpikeTimeList = []
            while(numberOfSpikesGenerated != numberOfSpikesPerNeuron):
                
                # Drawing excitation spikes from the chosen distribution
                if   typeOfDraw == 1:
                    spikeTime = random.gauss(meanTime,synchronny)
                elif typeOfDraw == 2:
                    spikeTime = random.expovariate(meanTime)
                elif typeOfDraw == 3:
                    spikeTime = random.uniform(a, b)
                elif typeOfDraw == 4:
                    spikeTime = invgauss.rvs(meanTime)
                else:
                    raise ValueError("Initial Spike Gen Error: your requested draw type does not exist")
                
                # Makes sure the spike time is not out of bounds
                if (spikeTime > 0 and spikeTime < self.duration):
                    
                    individualNeuronSpikeTimeList.append(spikeTime)
                    
                    # This if else structure is meant to help keep the csv format
                    # the first time a number is written to a file, on the number
                    # char is written. Spaces and commas are only added from the 
                    # second write and on.
                    if self.firstWrite == True:
                        excitationSpikesdataFile.write(str(spikeTime))
                        self.firstWrite = False
                    else:
                        excitationSpikesdataFile.write(", " + str(spikeTime))
                    individualNeuronSpikeTimeList.append(x)
                    self.arrayOfInitialSpikeTimes.append(individualNeuronSpikeTimeList)
                    numberOfSpikesGenerated += 1
                    individualNeuronSpikeTimeList = []
                    
            excitationSpikesdataFile.close()
            
        # Creating the dead copy, so we can still use the information even after the loss of the original spike time list
        # I am taking the copy before sorting because that way, its easier to process for the required data because it's
        # ordered for the neurons and not the spike times.
        self.deadCopyOfArrayOfInitialSpikeTimes = list(self.arrayOfInitialSpikeTimes)        
        # Sorts the array, but because the addresses have been appended to the 
        # spike times, the user specified distribution is still preserved. 
        self.arrayOfInitialSpikeTimes.sort()
        
     
        
        
    
    def getNeuronConnectionsOf(self,neuronid):
        '''
        Returns the actual neuron objects that are connected to the given id
        '''
        
        listOfConnectingNeuronsIDs = self.adjacencyMatrix.getConnectionIDs(neuronid)
        listOfConnectingNeurons = []
        for neuronID in listOfConnectingNeuronsIDs:
            listOfConnectingNeurons.append(self.listOfNeurons[neuronID-1])
        
        return listOfConnectingNeurons
    
    def getSpikeTimeList(self):
        '''
        returns the spike times of the fake neurons
        '''
        
        return self.deadCopyOfArrayOfInitialSpikeTimes
    
    
    def getNeuron(self,neuronID):
        '''
        returns the actual neuron object provided its id #
        '''
        
        # Uses the id's from 1 and on
        if neuronID == "last":
            return self.listOfNeurons[len(self.listOfNeurons)-1]
        else:
            return self.listOfNeurons[neuronID-1]
    
    
    
    
    def rapidAddNeuron(self, spikeAmplitude, spikeAmplitude2, exciteTimeConstant,inhibTimeConstant,designation,numberOfNeurons,delay1,delay2,wantBGCurrent):
        '''
        # Adds multiple neurons of the provided parameters to the network
        # and makes room for them in the adjacency matrix
        '''
        
        while (numberOfNeurons != 0):
            self.addNeuron(spikeAmplitude,spikeAmplitude2, exciteTimeConstant,inhibTimeConstant,designation,delay1,delay2,wantBGCurrent)
            numberOfNeurons -= 1
    
    
    
    
    def getAdjacencyMatrix(self):
        '''
        # Returns the adjancency matrix
        '''
        
        return self.adjacencyMatrix
    
    
    
    
    def fileWipe(self,name = False):
        '''
        This method will remove all voltage, current, and spike time data 
        for the real neurons
        '''
        if name == False:
            name = self.networkID
            
        for file in self.os.listdir(name+'/Time/'):
            self.os.remove(name+'/Time/' + file) 
            
        # Cleaning up the excitatory neuron data    
        for file in self.os.listdir(name+'/ExcitatoryNeuron/Voltage/'):
            self.os.remove(name+"/ExcitatoryNeuron/Voltage/" + file)   
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/InhibitoryCurrent/'):
            self.os.remove(name+'/ExcitatoryNeuron/InhibitoryCurrent/' + file)
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/ExcitatoryCurrent/'):
            self.os.remove(name+'/ExcitatoryNeuron/ExcitatoryCurrent/' + file)
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/SpikeTime/'):
            self.os.remove(name+'/ExcitatoryNeuron/SpikeTime/' + file)
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/BackgroundCurrent/'):
            self.os.remove(name+'/ExcitatoryNeuron/BackgroundCurrent/' + file) 
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/NMDA/'):
            self.os.remove(name+'/ExcitatoryNeuron/NMDA/' + file)  
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/NMDA_Diagnostics/NMDA_Cond/'):
            self.os.remove(name+'/ExcitatoryNeuron/NMDA_Diagnostics/NMDA_Cond/' + file)
            
        for file in self.os.listdir(name+'/ExcitatoryNeuron/NMDA_Diagnostics/ddeactdt/'):
            self.os.remove(name+'/ExcitatoryNeuron/NMDA_Diagnostics/ddeactdt/' + file)        
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/NMDA_Diagnostics/mgBlock/'):
            self.os.remove(name+'/ExcitatoryNeuron/NMDA_Diagnostics/mgBlock/' + file)        
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/NMDA_Diagnostics/dactdt/'):
            self.os.remove(name+'/ExcitatoryNeuron/NMDA_Diagnostics/dactdt/' + file)        
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/NMDA_Diagnostics/self.deact/'):
            self.os.remove(name+'/ExcitatoryNeuron/NMDA_Diagnostics/self.deact/' + file)   
            
        for file in self.os.listdir(name+'/ExcitatoryNeuron/NMDA_Diagnostics/self.act/'):
            self.os.remove(name+'/ExcitatoryNeuron/NMDA_Diagnostics/self.act/' + file)        
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/gate_diagnostics/h/'):
            self.os.remove(name+'/ExcitatoryNeuron/gate_diagnostics/h/' + file)  
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/gate_diagnostics/m/'):
            self.os.remove(name+'/ExcitatoryNeuron/gate_diagnostics/m/' + file)    
            
        for file in self.os.listdir(name+'/ExcitatoryNeuron/gate_diagnostics/n/'):
            self.os.remove(name+'/ExcitatoryNeuron/gate_diagnostics/n/' + file)     
    
        for file in self.os.listdir(name+'/ExcitatoryNeuron/gate_diagnostics/dv/'):
            self.os.remove(name+'/ExcitatoryNeuron/gate_diagnostics/dv/' + file)
            
        for file in self.os.listdir(name+'/ExcitatoryNeuron/gate_diagnostics/dv1/'):
            self.os.remove(name+'/ExcitatoryNeuron/gate_diagnostics/dv1/' + file)
            
        for file in self.os.listdir(name+'/ExcitatoryNeuron/gate_diagnostics/dv2/'):
            self.os.remove(name+'/ExcitatoryNeuron/gate_diagnostics/dv2/' + file) 
        
        for file in self.os.listdir(name+'/ExcitatoryNeuron/gate_diagnostics/dv3/'):
            self.os.remove(name+'/ExcitatoryNeuron/gate_diagnostics/dv3/' + file)             
    
        for file in self.os.listdir(name+'/ExcitatoryNeuron/gate_diagnostics/dv4/'):
            self.os.remove(name+'/ExcitatoryNeuron/gate_diagnostics/dv4/' + file)              
        
        
        # Cleaning up the inhibitory neuron data    
        for file in self.os.listdir(name+'/InhibitoryNeuron/Voltage/'):
            self.os.remove(name+"/InhibitoryNeuron/Voltage/" + file)   
        
        for file in self.os.listdir(name+'/InhibitoryNeuron/InhibitoryCurrent/'):
            self.os.remove(name+'/InhibitoryNeuron/InhibitoryCurrent/' + file)
        
        for file in self.os.listdir(name+'/InhibitoryNeuron/ExcitatoryCurrent/'):
            self.os.remove(name+'/InhibitoryNeuron/ExcitatoryCurrent/' + file)
        
        for file in self.os.listdir(name+'/InhibitoryNeuron/SpikeTime/'):
            self.os.remove(name+'/InhibitoryNeuron/SpikeTime/' + file)  
        
        for file in self.os.listdir(name+'/InhibitoryNeuron/BackgroundCurrent/'):
            self.os.remove(name+'/InhibitoryNeuron/BackgroundCurrent/' + file)
            
        for file in self.os.listdir(name+'/InhibitoryNeuron/NMDA/'):
            self.os.remove(name+'/InhibitoryNeuron/NMDA/' + file) 
            
        for file in self.os.listdir(name+'/InhibitoryNeuron/NMDA_Diagnostics/NMDA_Cond/'):
            self.os.remove(name+'/InhibitoryNeuron/NMDA_Diagnostics/NMDA_Cond/' + file)
            
        for file in self.os.listdir(name+'/InhibitoryNeuron/NMDA_Diagnostics/ddeactdt/'):
            self.os.remove(name+'/InhibitoryNeuron/NMDA_Diagnostics/ddeactdt/' + file)        
        
        for file in self.os.listdir(name+'/InhibitoryNeuron/NMDA_Diagnostics/mgBlock/'):
            self.os.remove(name+'/InhibitoryNeuron/NMDA_Diagnostics/mgBlock/' + file)        
        
        for file in self.os.listdir(name+'/InhibitoryNeuron/NMDA_Diagnostics/dactdt/'):
            self.os.remove(name+'/InhibitoryNeuron/NMDA_Diagnostics/dactdt/' + file)        
        
        for file in self.os.listdir(name+'/InhibitoryNeuron/NMDA_Diagnostics/self.deact/'):
            self.os.remove(name+'/InhibitoryNeuron/NMDA_Diagnostics/self.deact/' + file)   
            
        for file in self.os.listdir(name+'/InhibitoryNeuron/NMDA_Diagnostics/self.act/'):
            self.os.remove(name+'/InhibitoryNeuron/NMDA_Diagnostics/self.act/' + file)
            
        for file in self.os.listdir(name+'/InhibitoryNeuron/gate_diagnostics/h/'):
            self.os.remove(name+'/InhibitoryNeuron/gate_diagnostics/h/' + file)  
        
        for file in self.os.listdir(name+'/InhibitoryNeuron/gate_diagnostics/m/'):
            self.os.remove(name+'/InhibitoryNeuron/gate_diagnostics/m/' + file)    
            
        for file in self.os.listdir(name+'/InhibitoryNeuron/gate_diagnostics/n/'):
            self.os.remove(name+'/InhibitoryNeuron/gate_diagnostics/n/' + file)     
    
        for file in self.os.listdir(name+'/InhibitoryNeuron/gate_diagnostics/dv/'):
            self.os.remove(name+'/InhibitoryNeuron/gate_diagnostics/dv/' + file)     
            
        for file in self.os.listdir(name+'/InhibitoryNeuron/gate_diagnostics/dv1/'):
            self.os.remove(name+'/InhibitoryNeuron/gate_diagnostics/dv1/' + file)
            
        for file in self.os.listdir(name+'/InhibitoryNeuron/gate_diagnostics/dv2/'):
            self.os.remove(name+'/InhibitoryNeuron/gate_diagnostics/dv2/' + file) 
        
        for file in self.os.listdir(name+'/InhibitoryNeuron/gate_diagnostics/dv3/'):
            self.os.remove(name+'/InhibitoryNeuron/gate_diagnostics/dv3/' + file)             
    
        for file in self.os.listdir(name+'/InhibitoryNeuron/gate_diagnostics/dv4/'):
            self.os.remove(name+'/InhibitoryNeuron/gate_diagnostics/dv4/' + file)         
            
        # Cleaning up the decoder neuron data    
        for file in self.os.listdir(name+'/DecoderNeuron/Voltage/'):
            self.os.remove("DecoderNeuron/Voltage/" + file)   
        
        for file in self.os.listdir(name+'/DecoderNeuron/InhibitoryCurrent/'):
            self.os.remove(name+'/DecoderNeuron/InhibitoryCurrent/' + file)
        
        for file in self.os.listdir(name+'/DecoderNeuron/ExcitatoryCurrent/'):
            self.os.remove(name+'/DecoderNeuron/ExcitatoryCurrent/' + file)
        
        for file in self.os.listdir(name+'/DecoderNeuron/SpikeTime/'):
            self.os.remove(name+'/DecoderNeuron/SpikeTime/' + file)        
    
    
    
    
    def write(self):
        '''
        # Used to make every real neuron write their data to text files
        # Note, this call is not recursive. "neuron.write() is a
        # call to the "write()" method of the neuron class
        '''
        
        for neuron in self.listOfNeurons:
            neuron.write()
            
           
    
    
    def getFrequency(self,neuronID):
        '''
        # Returns the frequency of the requested neuron
        # For simplicity, the user can write "last" instead of providing
        # a number ID to get the frequency of the very last neuron aka
        # the decoder neuron
        '''
        if neuronID == "last":
            return self.listOfNeurons[len(self.listOfNeurons)-1].getFrequency()
        else:
            return self.listOfNeurons[neuronID-1].getFrequency()
    
    
    
    def getNumberOfSpikes(self,neuronID):
        if neuronID == "last":
            return self.listOfNeurons[len(self.listOfNeurons)-1].getNumberOfSpikes()
        else:
            return self.listOfNeurons[neuronID-1].getNumberOfSpikes()
        
        
    
    
    def reset(self):
        '''
        # Resets every neuron so that the neuron network cycle can go once again
        # This is not recursive. "neuron.reset()" is a call to the "neuron" classes
        # "reset()" function. 
        '''
        for neuron in self.listOfNeurons:
            neuron.reset()          
    
        self.numberOfStepsNeeded = 0
        self.actualStepsTaken = 0
        self.totalTime = 0
        self.stepsPer5Percent = 0
        self.percentComplete = 0           
        self.mailStation = []
        self.runFinished = False
    
    def importAdjacencyMatrix(self,adjacencyMatrix):
        
        '''
        # Allows the user to bring in an adjacency matrix
        '''
        
        self.adjacencyMatrix = adjacencyMatrix
        self.mapImported = True
        
    def calculateExternalConnections(self,dictOfInfo,listOfNetworks,neuronID):
        listOfConnections = []
        for adjMap in dictOfInfo[self.networkID][3]: # this is the mapbook
            listOfConnectingNeuronIDs = dictOfInfo[self.networkID][3][adjMap].getConnectionIDs(neuronID)
            for cnIDs in listOfConnectingNeuronIDs:
                neuronObj = listOfNetworks[dictOfInfo[adjMap][4]].getNeuron(cnIDs)
                listOfConnections.append(neuronObj)
        
        return listOfConnections
    
    def run(self, wantStatusReport):   
        '''
        # Run is basically what it is. It will step every single neuron through
        # all the required cycles and facillitate neuron to neuron communication
        # using the adjacnency matrix as a map.
        # Want status report is provided as boolean. Basically when true, it provides
        # a loading percentage. The option to turn it off is for the synchronny
        # because trials will "run" multiple times and then you would get 
        # "100% Complete" for each of the trials. 
        '''
        
        # Erase all previous data
        self.fileWipe()
        
        ## Makes sure you maintain a 1 to 1 excite and inhib ratio
        #if (len(self.arrayOfInitialSpikeTimes) + 1 != len(self.listOfNeurons)):
            #print("\nYour excitation to inhibition neuron ratio is not 1:1.")
            #return
        
        # Makes sure your time step is not larger than resolution
        if self.timeStep > self.resolution:
            print("\nYour time step cannot be greater than the resolution")
            return
        
        # This is basically a failsafe for when the user imports
        # their own adjacency matrix. It makes sure that the adjacency
        # matrix size is correct for the network.
        if len(self.adjacencyMatrix) != len(self.listOfNeurons):
            print("\nAdjacency Matrix and Neural Network mismatch")
            return
        
          
                
        
        # Calculates the number of times cycle has to be called
        # also used to calculate which inhibitory neuron to spike
        numberOfStepsNeeded = float(self.duration)/float(self.timeStep)
        actualStepsTaken = 0
        totalTime = 0
        stepsPer5Percent = numberOfStepsNeeded//19
        percentComplete = 0           
        mailStation = []  # Mail Station Structure: [[amp,from neuron designation,spiketime,to neuron],[amp,from neuron designation,spiketime, to neuron]...]
        
        # takes every neuron through a time step
        while (actualStepsTaken != numberOfStepsNeeded+1):
            
            ## Responsible for fake neuron current spikes
            # Basically, while the array spike times array is not empty and the 
            # spike time is less than or equal to the current time send a spike
            while (len(self.arrayOfInitialSpikeTimes) != 0 and self.arrayOfInitialSpikeTimes[0][0] <= totalTime):
                
                # If sending current to an excitatory neuron
                if self.listOfNeurons[self.arrayOfInitialSpikeTimes[0][1]].getDesignation() == 0:
                    
                    self.listOfNeurons[self.arrayOfInitialSpikeTimes[0][1]].sendConductance(self.fakeExciteAmplitude,0)
                
                # If sending current to an inhibitory neuron
                else:
                    self.listOfNeurons[self.arrayOfInitialSpikeTimes[0][1]].sendConductance(self.fakeInhibAmplitude,0)
                
                
                # Sends a spike to the decoder neuron. If the user wants, the decoder's spike amplitude
                self.listOfNeurons[len(self.listOfNeurons)-1].sendConductance(self.fakeExciteAmplitude,0)
                self.arrayOfInitialSpikeTimes.pop(0)
            
            # Redistribution of packages
            for mail in mailStation:
                mail[3].recieveSpike(mail)
            
            # Cycles every neuron 
            mailStation = []
            for neuron in self.listOfNeurons:
                neuron.cycle(self.getNeuronConnectionsOf(neuron.getID()), mailStation)
                      
            #if len(mailStation) != 0:
                #print(self.networkID+'/Hello')
            
            
            
            if actualStepsTaken%stepsPer5Percent == 0 and wantStatusReport == True:
                print(percentComplete, "% done")
                percentComplete += 5
                
            totalTime += self.timeStep
            actualStepsTaken += 1
            
        
        # Writes the data of each neuron to files in a csv format for matlab
        # to read
        self.write()
        
        if wantStatusReport == True:  
            print("100% Complete")         
        
    def runADt(self,dictOfInfo,listOfNetworks):   
        '''
        This method will basically allow the incrementation of the network by a single dt. The network can be run by 
        calling this method again and again. This will be used for incrementing multiple networks using another class as a manager. 
        
        mapbook and list of networks allow the nn to calculate external connections. Mail station is responsible for delivering spikes,
        and while each network is responsible for handling only its own neurons, mail station technically can deliver anywhere so we only have to add 
        the external neurons to the list of connections for cross network support
        '''
        
        # There are several values that need to be caclculated once and saved. This if structure allows this method to be called repeatedly
        # and have the values calculated only once and preserved 
        if self.actualStepsTaken == 0:
            
            # Erase all previous data
            self.fileWipe()
            if self.networkID == 'CA3':
                #self.adjacencyMatrix.presetThree(50,.25,.15,.25,.05)     # II,IE,EI,EE - Abid's
                if self.mapImported == False:
                    self.adjacencyMatrix.presetThree(50,.25,.20,.20,.05)
                
            elif self.networkID == 'CA1':
                self.adjacencyMatrix.presetThree(50,0,.35,0,0)
            else:
                raise ValueError('Unrecognized Network For Map Determination')
            
            ## Makes sure you maintain a 1 to 1 excite and inhib ratio
            #if (len(self.arrayOfInitialSpikeTimes) + 1 != len(self.listOfNeurons)):
                #print("\nYour excitation to inhibition neuron ratio is not 1:1.")
                #return
            
            
            
            # This is basically a failsafe for when the user imports
            # their own adjacency matrix. It makes sure that the adjacency
            # matrix size is correct for the network.
            if len(self.adjacencyMatrix) != len(self.listOfNeurons):
                print("\nAdjacency Matrix and Neural Network mismatch")
                return

            # Calculates the number of times cycle has to be called
            # also used to calculate which inhibitory neuron to spike
            self.numberOfStepsNeeded = float(self.duration)/float(self.timeStep)
            self.actualStepsTaken = 0
            self.totalTime = 0
            self.stepsPer5Percent = self.numberOfStepsNeeded//19
            self.percentComplete = 0           
            self.mailStation = []  # Mail Station Structure: [[amp,from neuron designation,spiketime,to neuron],[amp,from neuron designation,spiketime, to neuron]...]
            
            ## Responsible for fake neuron current spikes
            # Basically, while the array spike times array is not empty and the 
            # spike time is less than or equal to the current time send a spike
            while (len(self.arrayOfInitialSpikeTimes) != 0 and self.arrayOfInitialSpikeTimes[0][0] <= totalTime):
                
                # If sending current to an excitatory neuron
                if self.listOfNeurons[self.arrayOfInitialSpikeTimes[0][1]].getDesignation() == 0:
                    
                    self.listOfNeurons[self.arrayOfInitialSpikeTimes[0][1]].sendConductance(self.fakeExciteAmplitude,0)
                
                # If sending current to an inhibitory neuron
                else:
                    self.listOfNeurons[self.arrayOfInitialSpikeTimes[0][1]].sendConductance(self.fakeInhibAmplitude,0)
                
                
                # Sends a spike to the decoder neuron. If the user wants, the decoder's spike amplitude
                self.listOfNeurons[len(self.listOfNeurons)-1].sendConductance(self.fakeExciteAmplitude,0)
                self.arrayOfInitialSpikeTimes.pop(0)
            
            # Redistribution of packages
            for mail in self.mailStation:
                mail[3].recieveSpike(mail)
                
            for neuron in self.listOfNeurons:
                externalConnections = self.calculateExternalConnections(dictOfInfo,listOfNetworks,neuron.getID())
                neuron.setConnections(self.getNeuronConnectionsOf(neuron.getID())+externalConnections)            
            
            # Cycles every inhibitory neuron and the decoder neuron
            self.mailStation = []
            self.lastNeuronToFire = 0
            for neuron in self.listOfNeurons:
                #externalConnections = self.calculateExternalConnections(dictOfInfo,listOfNetworks,neuron.getID()) 
                neuron.cycle(self.mailStation,dictOfInfo[self.networkID][5])
                self.lastNeuronToFire += 1
            
            if self.actualStepsTaken%self.stepsPer5Percent == 0:
                print(self.percentComplete, "% done")
                self.percentComplete += 5
                
            self.totalTime += self.timeStep
            self.actualStepsTaken += 1            
        
        elif self.actualStepsTaken == self.numberOfStepsNeeded+1:
            # Writes the data of each neuron to files in a csv format for matlab
            # to read
            self.write()
            self.runFinished = True
            
        
        else:
            ## Responsible for fake neuron current spikes
            # Basically, while the array spike times array is not empty and the 
            # spike time is less than or equal to the current time send a spike
            while (len(self.arrayOfInitialSpikeTimes) != 0 and self.arrayOfInitialSpikeTimes[0][0] <= totalTime):
                
                # If sending current to an excitatory neuron
                if self.listOfNeurons[self.arrayOfInitialSpikeTimes[0][1]].getDesignation() == 0:
                    
                    self.listOfNeurons[self.arrayOfInitialSpikeTimes[0][1]].sendConductance(self.fakeExciteAmplitude,0)
                
                # If sending current to an inhibitory neuron
                else:
                    self.listOfNeurons[self.arrayOfInitialSpikeTimes[0][1]].sendConductance(self.fakeInhibAmplitude,0)
                
                
                # Sends a spike to the decoder neuron. If the user wants, the decoder's spike amplitude
                self.listOfNeurons[len(self.listOfNeurons)-1].sendConductance(self.fakeExciteAmplitude,0)
                self.arrayOfInitialSpikeTimes.pop(0)
            
            # Redistribution of packages
            for mail in self.mailStation:
                mail[3].recieveSpike(mail)
            
            # Cycles every neuron
            self.mailStation = []
            self.lastNeuronToFire = 0
            for neuron in self.listOfNeurons:
                #externalConnections = self.calculateExternalConnections(dictOfInfo,listOfNetworks,neuron.getID())
                #if len(externalConnections) != 0:
                    #print(self.networkID+'/hello')
                neuron.cycle(self.mailStation,dictOfInfo[self.networkID][5])
                self.lastNeuronToFire += 1
                      
           
            
            if self.actualStepsTaken%self.stepsPer5Percent == 0:
                print(self.percentComplete, "% done")
                self.percentComplete += 5
                
            self.totalTime += self.timeStep
            self.actualStepsTaken += 1
                
            
        
        
          
    def isRunFinished(self):
        return self.runFinished
    def getLastNeuronToFireID(self):
        return self.lastNeuronToFire
    def createDataLibrary(self):
        
        for neuron in self.listOfNeurons:
            self.listOfVoltages.append(neuron.getVoltageList())
        
        return self.listOfVoltages
    def getTimeList(self):
        return self.listOfNeurons[0].getTimeList()
    def setNMDAparam(self,b_act_step, b_deact_step,p_act_step, p_deact_step, tau_nmda_act_1, tau_nmda_deact_1, tau_nmda_act_2, tau_nmda_deact_2):
        for neuron in self.listOfNeurons:
            neuron.setNMDAparam(b_act_step, b_deact_step,p_act_step, p_deact_step, tau_nmda_act_1, tau_nmda_deact_1, tau_nmda_act_2, tau_nmda_deact_2)
        
    
if __name__ == '__main__':
    print("\n Wrong File!!")
    #import os
    #preset = 3

    ## Preset for 2 interneurons
    #if preset == 1: 

        ### Remember 0 for excitatory and 1 for inhibitory
        
        ## For the inhibition Neurons
        #spikeAmplitude      = 3           # Is 3x stronger than EPSP to P cell
        #exciteTimeConstant  = .00081      # Takes Approx 13 ms to decay halfway
        #inhibTimeConstant   = 5           ## As of now, there is no inhib to inhib, so this doesn't do anything
        #refractDuration     = 2           ## 
        #designation         = 1           # Good
        #numberOfNeurons     = 2          # Good
        #delay               = 1.5         # Adjusted
        
        
        ## For the decoder Neuron
        #spikeAmplitude2     = 1           ##
        #exciteTimeConstant2 = .000063     # Takes Approx 47 ms to decay halfway
        #inhibTimeConstant2  = .00034      # Takes Approx 20 ms to decay halfway
        #refractDuration2    = 2           ##
        #designation2        = 0           # Good
        #delay2              = 0           # Good
        
        
        ## Neuron Network Initial Parameters
        ### Excite Amplitude  = the current that the initial spikes send to an excitatory neuron
        ### Inhib Amplitude   = the current that the initial spikes send to an inhibitory neuron
        #duration         = 100            # Good
        #timeStep         = .1             # Good
        #exciteAmplitude  = 1              # Adjusted, is 1/3 of inhib amplitude
        #inhibAmplitude   = 1              ## Won't do anything at this time since no inhib to inhib neuron connections
        #resolution       = .2             # Good
        
        
        ## Excitation Spike Parameters
        ### Mean time and synchronny are used in the gauss distribution
        ### RangeMin and RangeMax are used during the uniform distribution draw
        ### Only mean time is used for the exponential distribution 
        ### Number of Spikes per Neuron is the ammount of times the excitatory neurons
        ###     will fire each. The ratios are kept even so one inhib cannot recieve 
        ###     3 spikes while another inhib recieves one spikes. Both would recive 
        ###     the same amount
        ### See line ??? to adjust ratio limiter
        #numberOfSpikesPerNeuron = 1
        #meanTime                = 30
        #synchronny              = 2
        #RangeMin                = 1
        #RangeMax                = 4
        #typeOfDraw              = 1
        ###    1      2       3           4
        ###  Gauss  Expon  Uniform  Inverse Gauss
        
        
        ## Synchronny Hz, Mean, Jitter Grapher UI
        #'''
        #Special Mapping = when on, the user must set the special mapping
        #order to select one of the other special maps
        #'''
        #Want_To_Run          = False
        #SyncMax              = 20  
        #Trials_Per_Sync      = 5      # Needed for Jitter Calculations
        #Special_Mapping      = True
        #Special_Mapping_Order = 0
        
  
    ## Preset for 100 neurons
    #elif preset == 2: #---------------------------------------------------------------------
          ##---------------------------------------------------------------------
          ##---------------------------------------------------------------------
        
        ### Remember 0 for excitatory and 1 for inhibitory
                
        ## For the inhibition Neurons
        #spikeAmplitude      = 3           # Is 3x stronge than EPSP to P cell
        #exciteTimeConstant  = .00081      # Takes Approx 13 ms to decay halfway
        #inhibTimeConstant   = 5           ## As of now, there is no inhib to inhib, so this doesn't do anything
        #refractDuration     = 1.75        ## 
        #designation         = 1           # Good
        #numberOfNeurons     = 100         # Good
        #delay               = 2         # Adjusted
        
        
        ## For the decoder Neuron
        #spikeAmplitude2     = 1           ##
        #exciteTimeConstant2 = .000063     # Takes Approx 47 ms to decay halfway
        #inhibTimeConstant2  = .00034      # Takes Approx 20 ms to decay halfway
        #refractDuration2    = 2           ##
        #designation2        = 0           # Good
        #delay2              = 0           # Good
        
        
        ## Neuron Network Initial Parameters
        ### Excite Amplitude  = the current that the initial spikes send to an excitatory neuron
        ### Inhib Amplitude   = the current that the initial spikes send to an inhibitory neuron
        #duration         = 100            # Good
        #timeStep         = .1             # Good
        #exciteAmplitude  = 1              # Adjusted, is 1/3 of inhib amplitude
        #inhibAmplitude   = 1              ## Won't do anything at this time since no inhib to inhib neuron connections
        #resolution       = .2             # Good
        
        
        ## Excitation Spike Parameters
        ### Mean time and synchronny are used in the gauss distribution
        ### RangeMin and RangeMax are used during the uniform distribution draw
        ### Only mean time is used for the exponential distribution 
        ### Number of Spikes per Neuron is the ammount of times the excitatory neurons
        ###     will fire each. The ratios are kept even so one inhib cannot recieve 
        ###     3 spikes while another inhib recieves one spikes. Both would recive 
        ###     the same amount
        ### See line 138 to remove ratio limiter
        #numberOfSpikesPerNeuron = 1
        #meanTime                = 30
        #synchronny              = 0
        #RangeMin                = 1
        #RangeMax                = 4
        #typeOfDraw              = 1
        ###    1      2       3           4
        ###  Gauss  Expon  Uniform  Inverse Gauss
        
        
        ## Synchronny Hz, Mean, Jitter Grapher UI
        #'''
        #Special Mapping = when on, the user must set the special mapping
        #order to select one of the other special maps
        #'''        
        #Want_To_Run           = False
        #SyncMax               = 10
        #Trials_Per_Sync       = 10     # Needed for Jitter Calculations
        #Special_Mapping       = False
        #Special_Mapping_Order = 0
    
    
    
    
    
    #else:
        #print("\n\n-----------------------\n      Wrong File\n-----------------------\n\n")
        #sys.exit(0)
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    ##----------------------------------------------------------------------------
    ##----------------------------------------------------------------------------
    ##----------------------------------------------------------------------------
    
    ### This section is a bunch of code to for the above interface to work. There is
    ### no need for the user to look at this area.
    
    ### Sets up a neuron network. Users can interface with it using the above section
    ## Neuron Network Initialization
    #neuralNetwork = Neuron_Network(duration,timeStep,exciteAmplitude,resolution, inhibAmplitude)
    
    ## Adding the inhibitory neurons
    #neuralNetwork.rapidAddNeuron(spikeAmplitude , exciteTimeConstant ,inhibTimeConstant ,refractDuration , designation , numberOfNeurons,delay)
    
    ## Adding the decoder neuron
    #neuralNetwork.addNeuron(     spikeAmplitude2, exciteTimeConstant2,inhibTimeConstant2,refractDuration2, designation2, delay2)
    
    ## Generating the initial excitation spike data
    #if Want_To_Run == False:
        #neuralNetwork.initialSpikeGeneration(numberOfSpikesPerNeuron, meanTime, synchronny,typeOfDraw,RangeMin,RangeMax)    
    
    
    ## Sets up the neuron network connection map in the adjacency Matrix. 
    ## Preset One means attach all neurons to the decoder.
    #x = neuralNetwork.getAdjacencyMatrix()
    #x.presetOne()
    #if Want_To_Run == False:
        ## True for print loading bar, false for holding jitter data
        #neuralNetwork.run(True)
    
        ## Prints the decoder Frequency
        #print("\nDecoder Frequency:", neuralNetwork.getFrequency("last"))





    ## Synchronny Hz Grapher Functions
    
    ## Creates a list of synchronnies. For examples, if the user
    ## enters 20, then it will make a list with 1, 2, 3, 4, ..., 20 
    #def setSynchronnyRange(Max):
        #synchronnyList = []
        #for x in range(0, Max+1):
            #synchronnyList.append(x)
        
        #return synchronnyList
    
    
    
    ## Responsible for writing the frequency data to a file
    ## Because synchronny always runs from 0 to some number n
    ## there is no need to store synchronny data because we can
    ## calculate which frequency data member belongs to which 
    ## synchronny
    #def FJM_Write(frequencyList,jitterList,meanList):
        
        ## Wipe all files in the Synchronny Folder
        #os.remove("SynchronnyVSHz/HzDataFile.txt")
        #os.remove("SynchronnyVSJitter/JitterDataFile.txt")
        #os.remove("SynchronnyVSMean/MeanDataFile.txt")
        
        ## Write the list of frequencies to the data file:
        #frequencyDataFile = open("SynchronnyVSHz/HzDataFile.txt", "w")
        #frequencyDataFile.write((self.networkID+'/[%s]' % ', '.join(map(str, frequencyList)))[1:-1])
        #frequencyDataFile.close()
        
        ## Write the list of jitter to the data file:
        #jitterDataFile = open("SynchronnyVSJitter/JitterDataFile.txt", "w")
        #jitterDataFile.write((self.networkID+'/[%s]' % ', '.join(map(str, jitterList)))[1:-1])
        #jitterDataFile.close()        
        
        ## Write the list of mean initial spike times to the data file
        #meanDataFile = open("SynchronnyVSMean/MeanDataFile.txt", "w")
        #meanDataFile.write((self.networkID+'/[%s]' % ', '.join(map(str, meanList)))[1:-1])
        #meanDataFile.close()        
        
        
        
    ## It will run the neural network for each synchronny in the synchronny list
    #def runSyncGrapher(Max):
        #'''
        #Master Frequency list = will store the mean frequencies for the respective 
                                #synchronnies
        #Synchronny list       = list of desired synchronnies to test
        #iterationPer5Percent  = calculation for loading bar
        #iteration             = tracker used for loading bar
        #percentComplete       = percent complete for loading bar
        #trialsDone            = used to run multiple trials per synchronny to 
                                #gather jitter data vs synchronny
        #small Frequency List  = stores the the frequencies for every trial per
                                #synchronny. The mean will be attached the the
                                #master list
        #'''
        
        #neuralNetwork.reset(False) # False meaning delete jitter data.
        #masterFrequencyList      = []
        #jitterDataList           = []
        #meanInitialSpikeTimeList = []
        #synchronnyList           = setSynchronnyRange(Max)
        #iterationsPer5Percent    = len(synchronnyList)//19
        #iteration                = 0
        #percentComplete          = 0   
        
        
        #if Special_Mapping == True:
            
            #if Special_Mapping_Order == 0:
                #pass
        #else:
            #x.presetOne()
        ## increment through every synchronny
        #for sync in synchronnyList:
            
            ## Stored here in order to be reset for every synchronny
            #smallFrequencyList = []
            #trialsDone         = 0
            #attempts           = 0
            #ignore             = False
            #warnAfterAttempts  = 10
            #neuralNetwork.reset(False) # Delete Jitter data because we are 
                                       ## collecting a new set for a new sync
            
            
            ## run multiple trials per synchronny to collect jitter data
            #while trialsDone != Trials_Per_Sync:
                #neuralNetwork.initialSpikeGeneration(numberOfSpikesPerNeuron, meanTime, sync,1,RangeMin,RangeMax)
                #neuralNetwork.run(False, True) # False meaning no loading bar
                                               ## True meaning store data necessary for jitter calculation
                
                
                ## As for frequency, we still want the data even if the decoder
                ## didn't spike, so it's left out of the trial discount mechanism
                #smallFrequencyList.append(neuralNetwork.getFrequency("last"))
                
                ## Ensures that we use trials in which the decoderSpikes
                ## The Jitter Data collector will not add any info for
                ## a trial that did not spike. As a result, all we have to
                ## do here is just not count that trial.
                #decoderNumberOfSpikes = neuralNetwork.getNumberOfSpikes("last")
                #if decoderNumberOfSpikes != 0:
                    
                    #trialsDone += 1
 
                #attempts += 1
                
                ## Reset the network to run again for the same synchronny
                #neuralNetwork.reset(True) # True meaning keep jitter data                   
                
                ## This structure is basically a check agains an infinite loop
                ## the thing is the sync data trials will keep running until all
                ## the trials wanted are done. But if the decoder never fires,
                ## then the user will get stuck. So after 10 failed trials, 
                ## the option to quit, snooze, or ignore becomes available
                #if trialsDone - attempts != 0 and ignore == False: 
                    #if abs(trialsDone - attempts) % warnAfterAttempts == 0:
                        #print("\nYou are on Synchronny:", sync, "\nThis is attempt:", attempts, "\nYour number of successful trials:", trialsDone, "\n\nYou might either be stuck in an infinite loop, or have a very low chance of the decoder firing. \nPlease choose an option:\n1) Quit and save what you have\n2) Try again " + str(warnAfterAttempts) + " more times \n3) Ignore this warning for the rest of this trial\n4) Warn you after x attempts\n(Enter 1, 2, 3, or 4)")
                        #x = int(input())
                        
                        #if x == 1:
                            #return
                        #elif x == 2:
                            #attempts = trialsDone
                        #elif x == 3:
                            #ignore = True
                        #elif x == 4:
                            #y = input("Warn you after how many inputs? Note: this will reset back to 10 for synchronny test\n")
                            #warnAfterAttempts = int(y)
                        #else:
                            #print("Thats not a valid option. Snoozing for " + str(warnAfterAttempts) + " more tries")
                            #attempts = trialsDone
            
            ## Appends the mean of the small frequency list
            #masterFrequencyList.append(statistics.mean(smallFrequencyList))
            #jitterDataList.append(neuralNetwork.getJitterData())
            #meanInitialSpikeTimeList.append(neuralNetwork.getJitterData("mean"))
            
                       
           
            ## Calculations for the loading bar
            #if iteration%iterationsPer5Percent == 0:
                #print(percentComplete, "% done")
                #percentComplete += 5            
            
            ## Incrementation for the loading bar
            #iteration += 1
            
            
        
        ## Writing sync data to file
        #FJM_Write(masterFrequencyList, jitterDataList, meanInitialSpikeTimeList)
        
        ## For the loading bar
        #print("\n\n-----------------------------------------------------------")
        #print("100% Complete")
        
        ## For Complementary Data on the final synchronny:
        #print("\nYour final data set is on synchronny:", synchronnyList[-1])
        #print("Decoder Frequency:", neuralNetwork.getFrequency("last"))
            
    ## This script is just for calling a pop-up window when sync grapher is on
    ## Sync grapher takes a while, and we need to make sure the user doesn't try
    ## to look for data before the process finishes. 
    #if Want_To_Run == True:
        #applescript = """
        #display dialog "Wait for the Go Message" ¬
        #with title "Warning" ¬
        #with icon caution ¬
        #buttons {"OK"}
        #"""
        
        #subprocess.call("osascript -e '{}'".format(applescript), shell=True)        
        #runSyncGrapher(SyncMax)
        
        #applescript = """
        #display dialog "Go" ¬
        #with title "Warning" ¬
        #with icon caution ¬
        #buttons {"OK"}
        #"""        
        #subprocess.call("osascript -e '{}'".format(applescript), shell=True)