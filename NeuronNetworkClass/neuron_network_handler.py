import os 
from neuron import Neuron
from Adjacency_Matrix import Adjacency_Matrix
from Neuron_Network import Neuron_Network
import statistics
import subprocess
import multiprocessing
import shutil
import time
from scipy.fftpack import fft
from scipy.signal import savgol_filter
import numpy
#from collections import OrderedDict
class neuron_network_handler():     
    
    def __init__(self,CA1_ON):
        '''
        list of networks contains all the network type objects.
        networks info is a dictionary that uses the network names as keys and each entry has info on each network
        '''
        self.listOfNetworks = []
        self.networksInfo = dict() # The structure:   " key : [number of inhibitory neurons, totalNumberOfNeurons, [[ii,ie,ee,ei, network name],[ii,ie,ee,ei, network name],...], [[1 -> 2],[ 1 -> 3], [1 -> 4],...], location in list ] "
        self.CA1_ON = CA1_ON
        self.gateDiagnostics = False
        
        # Data processing when we aren't using matlab 
        self.totalFieldPotential = []
        self.fourrierPower = []
        self.fourrierFrequency = []
        self.timeStep = 0
        self.listOfTime = []
        self.peakTimes = []
        self.ca3ListOfCyclesPerHz = []
        self.totalHzList = [] 
        self.smoothFieldPotential = []
        
    
    # Must create all networks before assigning probabilities
    def createNetwork(self, networkID,                                                                                                                                # The network id
                 spikeAmplitude1,spikeAmplitude2,exciteTimeConstant,inhibTimeConstant,designation,numberOfNeurons,delay1,delay2,wantBGCurrent1,       # Inhib Neurons
                 spikeAmplitude3,spikeAmplitude4,exciteTimeConstant2,inhibTimeConstant2,designation2,numberOfNeurons2,delay3,delay4,wantBGCurrent2,  # Excit Neurons                                                                                            # Delays Times
                 wantADecoder, spikeAmplitude5,spikeAmplitude6, exciteTimeConstant3,inhibTimeConstant3,designation3,delay5,delay6,                   # Decoder
                 duration,timeStep,resolution,                                                                                                                        # Netowork Time Parameters
                 wantFakeNeurons,exciteAmplitude,inhibAmplitude,numberOfSpikesPerNeuron,meanTime,synchronny,RangeMin,RangeMax,typeOfDraw,                             # Fake Neurons
                 adjacencyMatrix = False                                                                                                                              # Adjacency Matrix is False by default  
                 ):                                                                                                                                                                         
        
        '''
        # These are all the initial variables needed to create the Neural Network
        # These vars don't include the info needed for the synchronny trials. Those
        # will be added and adjusted via additional methods
        
        See next to the variables for a description of what each is
        '''
        
        
        #self.spikeAmplitude1         = spikeAmplitude1         # Real neuron spike amplitude sent to excitatory neurons
        #self.spikeAmplitude2         = spikeAmplitude2         # Real neuron spike amplitude sent to inhibitory neurons
        #self.exciteTimeConstant      = exciteTimeConstant      # Real neuron excitatory current decay time     
        #self.inhibTimeConstant       = inhibTimeConstant       # Real neuron inhibitory current decay time   
        #self.refractDuration         = refractDuration         # Real neuron refractory duration  
        #self.designation             = designation             # Real neuron designation
        #self.numberOfNeurons         = numberOfNeurons         # Real neuron, number of inhibitory interneurons
        #self.delay1                  = delay1                  # Real neuron AP delay: From I to I
        #self.delay2                  = delay2                  # Real neuron AP delay: From I to E
        
        #self.spikeAmplitude3         = spikeAmplitude3         # Real neuron spike amplitude sent to excitatory neurons
        #self.spikeAmplitude4         = spikeAmplitude4         # Real neuron spike amplitude sent to inhibitory neurons
        #self.exciteTimeConstant2     = exciteTimeConstant2     # Real neuron excitatory current decay time     
        #self.inhibTimeConstant2      = inhibTimeConstant2      # Real neuron inhibitory current decay time   
        #self.refractDuration2        = refractDuration2        # Real neuron refractory duration  
        #self.designation2            = designation2            # Real neuron designation
        #self.numberOfNeurons2        = numberOfNeurons2        # Real neuron, number of inhibitory interneurons
        #self.delay3                  = delay3                  # Real neuron AP delay: From E to I
        #self.delay4                  = delay4                  # Real neuron AP delay: From E to E        
        
        #self.wantADecoder = wantADecoder
        #if self.wantADecoder == True:
            #raise ValueError('At this time, the decoder neuron function is down for maintentance')
            #self.spikeAmplitude5         = spikeAmplitude5         # Real neuron spike amplitude sent to excitatory neurons by decoder
            #self.spikeAmplitude6         = spikeAmplitude6         # Real neuron spike amplitude sent to excitatory neurons by decoder
            #self.exciteTimeConstant3     = exciteTimeConstant3     # Real neuron decoder excitatory current decay time
            #self.inhibTimeConstant3      = inhibTimeConstant3      # Real neuron decoder inhibitory current decay time
            #self.refractDuration3        = refractDuration3        # Real neuron decoder refractory period   
            #self.designation3            = designation3            # Real neuron decoder designation 
            #self.delayDecoder1           = delay5                  # Real neuron decoder delay: to I
            #self.delayDecoder2           = delay6                  # Real neuron decoder delay: to E
        
        #self.duration                = duration                # Neural Network duration
        #self.timeStep                = timeStep                # Neural Network timestep (aka sample rate)
        #self.resolution              = resolution              # Neural Network data collection resolution
        
        #self.wantFakeNeurons = wantFakeNeurons
        #if self.wantFakeNeurons == True:
            #raise ValueError('At this time, the fake neuron function is down for maintentance')
            #self.fakeExciteAmplitude     = exciteAmplitude         # Fake neuron spike amplitude sent to excitatory neurons
            #self.fakeInhibAmplitude      = inhibAmplitude          # Fake neuron spike amplitude sent to inhibitory neurons
            #self.numberOfSpikesPerNeuron = numberOfSpikesPerNeuron # Number of fake neurons per inhibitory interneurons
            #self.meanTime                = meanTime                # Mean time at which fake neurons spike  (applies for certain distributions only)
            #self.synchronny              = synchronny              # Standard deviation of when fake neurons spike  (applies for certain distributions only)
            #self.RangeMin                = RangeMin                # Minimum time of when fake neurons spike  (applies for certain distributions only)
            #self.RangeMax                = RangeMax                # Maximum time of when fake neurons spike  (applies for certain distributions only)
            #self.typeOfDraw              = typeOfDraw              # Type of distribution for generating fake neuron spike times
        
        #self.specialMap              = False                   # Boolean value that tell whether a custom map was uploaded or not
        
        
        # System of checks against bad input:
        #if wantFakeNeurons == False and numberOfNeurons2 == 0:
            #raise ValueError("You have no mode of excitement")
        
        if wantFakeNeurons == True and numberOfNeurons2 > 0:
            raise ValueError("You can't have both fake and real excitatory neurons")        
        
        # Makes sure your time step is not larger than resolution
        if timeStep > resolution:
            raise ValueError("\nYour time step cannot be greater than the resolution")        
        
        # For data library calculation
        self.timeStep = timeStep
        
        # Neuron Network Initialization
        neuralNetwork = Neuron_Network(duration,timeStep,resolution,networkID)
        
        # Adding the inhibitory neuron
        if numberOfNeurons > 0:
            neuralNetwork.rapidAddNeuron(spikeAmplitude1, spikeAmplitude2, exciteTimeConstant ,inhibTimeConstant ,designation , numberOfNeurons,delay1,delay2,wantBGCurrent1)
        
        # Adding the excitatory neuron      
        if numberOfNeurons2 > 0:
            neuralNetwork.rapidAddNeuron(spikeAmplitude3, spikeAmplitude4, exciteTimeConstant2 ,inhibTimeConstant2 , designation2 , numberOfNeurons2,delay3,delay4,wantBGCurrent2)
        
        # Adding the decoder neuron
        if wantADecoder == True:
            neuralNetwork.addNeuron(spikeAmplitude5, spikeAmplitude6, exciteTimeConstant3,inhibTimeConstant3,designation3, delayDecoder1, delayDecoder2,False)   
        
        # Adding fake neuron spikes
        if wantFakeNeurons == True:
            self.neuralNetwork.initialSpikeGeneration(numberOfSpikesPerNeuron, meanTime, synchronny,typeOfDraw,RangeMin,RangeMax,self.fakeExciteAmplitude, self.fakeInhibAmplitude)        
        
        # Adding the Adjacency Matrix
        if adjacencyMatrix != False:
            neuralNetwork.importAdjacencyMatrix(adjacencyMatrix)
    
        
        # Making filling out the info dictionary
        self.networksInfo[networkID] = [numberOfNeurons,neuralNetwork.getNumberOfNeurons(),[],dict(),len(self.listOfNetworks),dict()]
                                       #      0                       1                     2   3               4             5
        
        # 0 = Number of inhibitory neurons
        # 1 = Total Number Of Neurons
        # 2 = list of connection probabilities (see update info dict for more info on this structure)
        # 3 = Map book (more on structure below)
        # 4 = Index of the network in the list of networks
        # 5 = Spike parameters ie what amplitude and time constant to use for out network neurons
                        
            #             MAPBOOK STRUCTURE
            # Mapbook = to network name : adjacency map
            #           to network name : adjacency map
            #                        ...
            
            # Note: adjacency maps refer to forward going maps, ie the neurons connection probability 
            # towards sending signals, not recieving them
            
            # The mapbook essentially has a map for every network that the from network is connected to
            
            
            #          SPIKE PARAMETER STRUCTURE
            # spikeParams = to network name : [iTOiAmp,iTOeAmp,eTOiAmp,eTOeAmp,i_inhibCurrentTau,i_exciteCurrentTau,e_inhibCurrentTau,e_exciteCurrentTau]
            #               to network name : [iTOiAmp,iTOeAmp,eTOiAmp,eTOeAmp,iTOiTau,iTOeTau,eTOiTau,eTOeTau]

        self.listOfNetworks.append(neuralNetwork)
    
    def mapGenerator(ii,ie,ee,ei,fromSize,toSize,howManyFromINs,howManyToINs):
        '''
        Gets all the info needed for preset 4 of adjacency matrix, and applies preset 4
        '''
        
        adjacencyMatrix =  Adjacency_Matrix(True)
        adjacencyMatrix.generate(fromSize,toSize)
        adjacencyMatrix.presetFour(ii,ie,ei,ee,howManyFromINs,howManyToINs)
        return adjacencyMatrix
    
    
    def updateInfoDictionary(self,networkID,listOfConnectionProbabilities,spikeParameters): # [[ii,ie,ee,ei, network name],[ii,ie,ee,ei, network name],...]
                                                                                            # is the structure for list of connection pr's 
                                                                                            # Spike Parameters looks like:  [[network name,iTOiAmp,iTOeAmp,eTOiAmp,eTOeAmp,iTOiTau,iTOeTau,eTOiTau,eTOeTau],[network name,iTOiAmp,iTOeAmp,eTOiAmp,eTOeAmp,iTOiTau,iTOeTau,eTOiTau,eTOeTau],.....]
                                        
        '''
        This method is used to assign important cross network communication variables. This method will
        also create the maps needed for cross network communication
        '''
                                              
        if len(listOfConnectionProbabilities) != len(self.listOfNetworks) - 1:
            raise ValueError('Number of Networks and Connection Probabilities Mismatch')
        
        if len(spikeParameters) != len(self.listOfNetworks) - 1:
            raise ValueError('Number of Networks and Number of Spike Parameters Mismatch')        
        
        self.networksInfo[networkID][2] = listOfConnectionProbabilities
        
        for parameter in spikeParameters:
            self.networksInfo[networkID][5][parameter[0]] = parameter[1:]
         
        # the map is a forward map, ie the neurons connection probability towards sending signals, not recieving them
        for cP in listOfConnectionProbabilities:
            self.networksInfo[networkID][3][cP[4]] = neuron_network_handler.mapGenerator(cP[0],cP[1],cP[2],cP[3],self.networksInfo[networkID][1],self.networksInfo[cP[4]][1],self.networksInfo[networkID][0],self.networksInfo[cP[4]][0])
            # First part is the location of map book in the dictionary of info, and the second part is the map generator function running
            
    def printMap(self):
        '''
        Will print all the maps in map book. 
        '''
        for key in self.networksInfo:
            
            for networkMap in self.networksInfo[key][3]:
                print(key + ' to ' + networkMap)
                print(self.networksInfo[key][3][networkMap])
    

    def normalRun(self): # provide access to all adjacency maps and nn? That way the network can calculate what each neuron connects to?
        '''
        # Normal run is basically just running once and not worrying about the syncTrials 
        # As of the latest update, the sync trials will no longer work. Normal run will be 
          the only working option
          
          How it works: so basically neuron networks is responsible for managing its own network.
          Neuron network has a method called run a dt which basially runs the network for one time step
          The data is not deleted so essentially the simulation can be run by calling run a dt as many
          times as needed. Next for cross network communication, all the information needed is in the 
          info dictionary. For each run a dt, the dictionary is passed to the neuron network. the network
          program will find it's out network connections and then attach those to the list of thing to
          send spikes to. From there the communication is done and it all runs like normal. this works
          because the send spike que uses the actual neuron objects for sending data. 
        '''
        
        if len(self.listOfNetworks) == 0:
            raise ValueError("You have 0 Networks")
        
        # The actual running part
        
        # Calculations for finding presynaptic CA3 Neurons.
        # Basically for every CA1 neuron, we will find which
        # CA3 neurons send current to it using the reverse 
        # connections method. Then we will print it out
        
        # Clear the directories first
        for file in os.listdir('CA1'+'/ExcitatoryNeuron/ReverseConnections/'):
            os.remove('CA1'+"/ExcitatoryNeuron/ReverseConnections/" + file)   
        
        for file in os.listdir('CA1'+'/InhibitoryNeuron/ReverseConnections/'):
            os.remove('CA1'+"/InhibitoryNeuron/ReverseConnections/" + file)           
        
        if self.CA1_ON == True:
        # For the inhib neurons
            for x in range(1, self.networksInfo['CA1'][0]+1): # For every neuron id in range 1 to number of inhib neurons in CA1
                a = self.networksInfo['CA3'][3]['CA1'].getReverseConnection(x) # Get reverse connections of CA3 -> CA1 Map for id x
                connectionDataFile = open("CA1"+"/InhibitoryNeuron/ReverseConnections/neuron" + str(x) + ".txt", "w")
                connectionDataFile.write(('[%s]' % ', '.join(map(str, a)))[1:-1])
                connectionDataFile.close()
            
            # For the excite neurons
            for x in range(self.networksInfo['CA1'][0]+1, self.networksInfo['CA1'][1]+1): # For every neuron id in range first pn to number of all neurons in CA1
                a = self.networksInfo['CA3'][3]['CA1'].getReverseConnection(x) # Get reverse connections of CA3 -> CA1 Map for id x
                connectionDataFile = open("CA1"+"/ExcitatoryNeuron/ReverseConnections/neuron" + str(x) + ".txt", "w")
                connectionDataFile.write(('[%s]' % ', '.join(map(str, a)))[1:-1]) 
                connectionDataFile.close()
        else:
            self.clearNetwork('CA1')
        self.printMap()
        
        self.gateDiagnostics = self.isGateDiagnosticsOn()
        
        while self.listOfNetworks[0].isRunFinished() == False:
            if self.gateDiagnostics == True:
                try:
                    for neuronNet in self.listOfNetworks:
                        neuronNet.runADt(self.networksInfo,self.listOfNetworks)
                except ArithmeticError:
                    for neuronNet in self.listOfNetworks:
                        neuronNet.write()
                    lastNrn = self.listOfNetworks[0].getLastNeuronToFireID()
                    if lastNrn <= 50:   
                        print('The last neuron to properly fire was: IN ' + str(lastNrn))
                    else:
                        print('The last neuron to properly fire was: PN ' + str(lastNrn-50))
                    raise RuntimeError('You Had an arithmetic error. Data has been force written. Go check it out!')
            else:
                for neuronNet in self.listOfNetworks:
                    neuronNet.runADt(self.networksInfo,self.listOfNetworks)                
        
        # Makes sure that every neuron network is finished. All should be done.
        for neuronNet in self.listOfNetworks:
            if neuronNet.isRunFinished() == False:
                raise ValueError("Multi Network Time Mismatch. " + str(neuronNet.getNetworkID()) + " did not finish")
        
        ## Final touches upon finishing
        # Experimental Data Library:
        self.createDataLibrary()
        #self.writeDataLibrary()
        
        # Sound Alert that we are done
        os.system('afplay /System/Library/Sounds/Glass.aiff')        
        
        ## Prints the decoder Frequency
        #if self.wantADecoder == True:
            #print("\nDecoder Frequency:", self.neuralNetwork.getFrequency("last"))        
        
        # Prints a finishing message
        print("\n!!! Run Finished !!!\n")
    
   
      
    def setSynchronnyRange(self,Min, Max):
        '''
        # Creates a list of synchronnies. For examples, if the user
        # enters 20, then it will make a list with 1, 2, 3, 4, ..., 20   
        '''
        synchronnyList = []
        for x in range(Min, Max+1):
            synchronnyList.append(x)
        return synchronnyList    
    
    
    
   
    def FJM_Write(self,frequencyList,jitterList,meanList):
        '''
        # Responsible for writing the frequency data to a file
        # Because synchronny always runs from 0 to some number n
        # there is no need to store synchronny data because we can
        # calculate which frequency data member belongs to which
        # synchronny
        '''
        
        # Wipe all files in the Synchronny Folder
        os.remove("SynchronnyVSHz/HzDataFile.txt")
        os.remove("SynchronnyVSJitter/JitterDataFile.txt")
        os.remove("SynchronnyVSMean/MeanDataFile.txt")
        
        # Write the list of frequencies to the data file:
        frequencyDataFile = open("SynchronnyVSHz/HzDataFile.txt", "w")
        frequencyDataFile.write(('[%s]' % ', '.join(map(str, frequencyList)))[1:-1])
        frequencyDataFile.close()
        
        # Write the list of jitter to the data file:
        jitterDataFile = open("SynchronnyVSJitter/JitterDataFile.txt", "w")
        jitterDataFile.write(('[%s]' % ', '.join(map(str, jitterList)))[1:-1])
        jitterDataFile.close()
        
        # Write the list of mean initial spike times to the data file
        meanDataFile = open("SynchronnyVSMean/MeanDataFile.txt", "w")
        meanDataFile.write(('[%s]' % ', '.join(map(str, meanList)))[1:-1])
        meanDataFile.close() 
       
    
    
    
    def ensure_dir(self,file_path):
        '''
        # Essentially, a file folder path will be provided to this function.
        # If the folder exists, it will delete it (and all contents) then recreate it
        # If it doesn't exist, then it will create the folder.
        '''
        if os.path.exists(file_path):
            shutil.rmtree(file_path)
            os.makedirs(file_path)
        else:
            os.makedirs(file_path)

        
    
    
    
    def eSpikesAndINSpikesHarvester(self,synchronny,storeHowManyInterNeurons):
        '''
        # Special method for harvesting the IN spikes and E spikes data for the 
        # synchronny trials. This method will also be responsible for writing the
        # data to the file rather than FJM write because writing these ones is 
        # kinda complex.
        '''
        
        if storeHowManyInterNeurons > self.neuralNetwork.getNumberOfNeurons()-1:
            raise ValueError("You requested more rasters than there are interneurons")
        
        
        # Already has all the excite spikes for a given synchronny
        listOfExciteSpikes = self.neuralNetwork.getSpikeTimeList()
        
        # Unfortunately we can't just collect all the inter N spikes with one method.
        # We have to iterate through every neuron to collect each individiual interN's
        # spikes. We will collect up to the amount the user specifies, and at the end
        # we will add the decoder neuron
        listOfListOfinterNSpikes = []
        # Adding the interneurons to the list
        for neuron in range(0,storeHowManyInterNeurons): # we don't need the very last neuron (its the decoder)
            listOfListOfinterNSpikes.append(self.neuralNetwork.getNeuron(neuron).getSpikeTimes())
        
        # Adding the decoder to the list
        listOfListOfinterNSpikes.append(self.neuralNetwork.getNeuron(self.neuralNetwork.getNumberOfNeurons()-1).getSpikeTimes())
            
        
            
        # Writing to the Folders
        eSpikesFolderPath = "/Users/Abid/Documents/MATLAB/NeuronNetworkClass/SynchronnyESpikes/S" + str(synchronny)
        iSpikesFolderPath = "/Users/Abid/Documents/MATLAB/NeuronNetworkClass/SynchronnyINSpikes/S" + str(synchronny)
        
        self.ensure_dir(eSpikesFolderPath)
        self.ensure_dir(iSpikesFolderPath)
        
        
        
        
        # Write the IN spikes to the data file
        iterationNumber = 1 # Used to help name the files
        for spikeTimeListOfInterneuron in listOfListOfinterNSpikes:
            iSpikesFilePath = "SynchronnyINSpikes/S" + str(synchronny) + "/IN" + str(iterationNumber) + ".txt"
            inDataFiles = open( iSpikesFilePath , "w")
            inDataFiles.write(('[%s]' % ', '.join(map(str, spikeTimeListOfInterneuron)))[1:-1])
            inDataFiles.close()
            iterationNumber += 1
            
        # Write the E spikes to the data file
        for x in range(0,storeHowManyInterNeurons):
            eSpikesFilePath = "SynchronnyESpikes/S" + str(synchronny) + "/E" + str(x+1) + ".txt"
            EDataFiles = open(eSpikesFilePath, "w")
            firstWrite = True
            for y in range(0,self.numberOfSpikesPerNeuron):
                time = listOfExciteSpikes[x*self.numberOfSpikesPerNeuron+y][0]
                neuronToGoTo = listOfExciteSpikes[x+y][1]
                if firstWrite == True:
                    EDataFiles.write(str(time))
                    firstWrite = False
                else:
                    EDataFiles.write(", " + str(time))                
            EDataFiles.close()
   
    def didNeuronSpike(self, neuronID):
        if neuronID == 'last':
            decoder = self.neuralNetwork.getNeuron('last')
            if decoder.getNumberOfSpikes() != 0:
                return True
            else:
                return False
        else:
            someNeuron = self.neuralNetwork.getNeuron(neuronID)
            if someNeuron.getNumberOfSpikes() != 0:
                return True
            else:
                return False            
    
    
    # This will find the average of the positional values of a list
    # of lists. Ie if there are 5 lists, it will find the average of 
    # position 1, position 2, etc. It handles uneven sized lists by
    # simply not including the posistions not present in the average
    # calculation 
    
    ## This function is no longer used
    
    #def matrixAvg(self,unevenMatrix):
        #listOfAvgs = []
        
        ## Finding the length of the longest list
        #maxLength = 0
        #for lists in unevenMatrix:
            #if len(lists) > maxLength:
                #maxLength = len(lists)
        
        
        #for position in range(0,maxLength):
            #positionSum   = 0
            #count         = 0
            #for lists in unevenMatrix:
                ## Since the matrix is uneven we are just going to iterate
                ## by the length of the longest list so we need to make sure
                ## we don't overiterate 
                #if position < len(lists):
                    #positionSum   += float(lists[position])
                    #count += 1
            
            #listOfAvgs.append(positionSum/count)
        
        #return listOfAvgs
                
         
          
    
    
    def syncRun(self,syncMin, syncMax, trialsPerSync, storeHowManyInterNeurons, specialMap = False):
        
        applescript = """
        display dialog "Wait for the Go Message" ¬
        with title "Warning" ¬
        with icon caution ¬
        buttons {"OK"}
        """
        subprocess.call("osascript -e '{}'".format(applescript), shell=True)
        
        # If structure so that we can do stuff with special maps
        # Note specialMap and self.specialMap are two different vars. In the case
        # of specialMap, True means the user wants to use a custom map. In the case
        # of self.specialMap, True means a map was uploaded, and false means one
        # wasn't. Essentially its a check that the user actually uploaded a map
        if specialMap == False:
            nn_internal_matrix = self.neuralNetwork.getAdjacencyMatrix()
            nn_internal_matrix.presetOne()
        else:
            if self.specialMap == False:
                print("\n-----------------------------------\nYou Haven't Uploaded a Special Map!\n-----------------------------------\n\n")
                return
        
        
        '''
        Master Frequency list = will store the mean frequencies for the respective 
                                synchronnies
        Synchronny list       = list of desired synchronnies to test
        iterationPer5Percent  = calculation for loading bar
        iteration             = tracker used for loading bar
        percentComplete       = percent complete for loading bar
        trialsDone            = used to run multiple trials per synchronny to 
                                gather jitter data vs synchronny
        small Frequency List  = stores the the frequencies for every trial per
                                synchronny. The mean will be attached the the
                                master list
        maxIgnore             = basically a variable that allows a hidden 
                                option to avoid all warning of a low decoder
                                firing rate or possible infinite loop given
                                that at least one trial was succesful. this
                                is accomplished if the user types in '96' 
                                when a prompt is given about what to do. so
                                if a trial fails ten times, and the decoder
                                never spiked even once, only then will the user
                                get a warning
        Jitter Data List      = Var that means standard deviation of the firts
                                spike time of the decoder. For every trial per
                                sync (ones in which the decoder spiked) we
                                store the intial spike times of the decoder
                                and after the trials (ie exit the while loop)
                                we calculate the mean and standard dev. The
                                Standard dev gets appended to jitter
        Mean Initial....List  = (see above description) the means get appended
                                here
        '''  
 
        SyncMax                     = syncMax
        SyncMin                     = syncMin
        Trials_Per_Sync             = trialsPerSync        
        synchronnyList              = self.setSynchronnyRange(SyncMin,SyncMax)
        iterationsPer5Percent       = 1 #len(synchronnyList)//19
        iteration                   = 0
        percentComplete             = 0
        
        masterFrequencyList         = []
        jitterDataList              = []
        meanInitialSpikeTimeList    = []
        listOfListOfinterNSpikes    = []
        listOfExciteSpikes          = []

        maxIgnore                   = False
        
        
        # increment through every synchronny
        for sync in synchronnyList:
            
            # Stored here in order to be reset for every synchronny
            smallFrequencyList            = []
            trialsDone                    = 0
            attempts                      = 0
            ignore                        = False
            warnAfterAttempts             = 20
            decoderFirstSpikeTimes        = []
               
            
            
            # run multiple trials per synchronny to collect jitter data
            while trialsDone != Trials_Per_Sync:
                self.neuralNetwork.reset()
                self.neuralNetwork.initialSpikeGeneration(self.numberOfSpikesPerNeuron, self.meanTime, sync, 1, 0, 0)
                self.neuralNetwork.run(False) # False meaning no loading bar
                
                # As for frequency, we still want the data even if the decoder
                # didn't spike, so it's left out of the trial discount mechanism
                smallFrequencyList.append(self.neuralNetwork.getFrequency("last"))
                
                # Ensures that we use trials in which the decoderSpikes
                if self.didNeuronSpike('last') == True:                    
                    trialsDone += 1
                    decoderFirstSpikeTimes.append(float(self.neuralNetwork.getNeuron('last').getFirstSpikeTime()))
                    
                    
                attempts += 1
                
                
                # This structure is basically a check agains an infinite loop
                # the thing is the sync data trials will keep running until all
                # the trials wanted are done. But if the decoder never fires,
                # then the user will get stuck. So after 10 failed trials, 
                # the option to quit, snooze, or ignore becomes available
                if trialsDone - attempts != 0 and ignore == False:
                    
                    if maxIgnore == False or trialsDone == 0: 
                        
                        if abs(trialsDone - attempts) % warnAfterAttempts == 0:
                            os.system('afplay /System/Library/Sounds/Glass.aiff')
                            print("\nYou are on Synchronny:", sync, "\nThis is attempt:", attempts, "\nYour number of successful trials:", trialsDone, "\n\nYou might either be stuck in an infinite loop, or have a very low chance of the decoder firing. \nPlease choose an option:\n1) Quit and save what you have\n2) Try again " + str(warnAfterAttempts) + " more times \n3) Ignore this warning for the rest of this trial\n4) Warn you after x attempts\n(Enter 1, 2, 3, or 4)")
                            x = int(input())
                            
                            if x == 1:
                                return
                            elif x == 2:
                                attempts = trialsDone
                            elif x == 3:
                                ignore = True
                            elif x == 4:
                                y = input("Warn you after how many inputs? Note: this will reset back to 10 for synchronny test\n")
                                warnAfterAttempts = int(y)
                            
                            # Hidden option from user
                            elif x == 96:
                                maxIgnore = True
                            else:
                                print("Thats not a valid option. Snoozing for " + str(warnAfterAttempts) + " more tries")
                                attempts = trialsDone
            
            
            # We collect the excite spikes and inter N spikes here since we only need one succesful trial per synchronny
            self.eSpikesAndINSpikesHarvester(sync,storeHowManyInterNeurons)
            
            
            
            
            
            
            
            
            
            
            
                
            # Calculations for the loading bar
            if iteration%iterationsPer5Percent == 0:
                print(percentComplete, "% done")
                percentComplete += 5            
            iteration += 1   # Incrementation for the loading bar            
            
            # Appends the mean of the small frequency list
            masterFrequencyList.append(statistics.mean(smallFrequencyList))
            jitterDataList.append(statistics.stdev(decoderFirstSpikeTimes))
            meanInitialSpikeTimeList.append(statistics.mean(decoderFirstSpikeTimes))
            
                       
            
            #for x in range(0, storeHowManyInterNeurons):
                
                #for spikeTime in smallListOfListOfinterNSpikes[x]:
                    
                    #if spikeTime != -1:
                        #listOfListOfinterNSpikes[x].append(spikeTime)
                        #listOfListOfinterNSpikes[x].append(sync)
            
            
            #tempList = []
            #for spikeTime in smallListOfExciteSpikes:
                #tempList.append(spikeTime[0])
                #tempList.append(spikeTime[1])
                
            #listOfExciteSpikes.append(tempList)
            
        
        # Writing sync data to file
        self.FJM_Write(masterFrequencyList, jitterDataList, meanInitialSpikeTimeList)
        
        
        ## Final touches upon finishing
        # Sound Alert that we are done
        os.system('afplay /System/Library/Sounds/Glass.aiff')        
        applescript = """
        display dialog "Go" ¬
        with title "Warning" ¬
        with icon caution ¬
        buttons {"OK"}
        """        
        subprocess.call("osascript -e '{}'".format(applescript), shell=True)        
        # For the loading bar
        print("\n\n-----------------------------------------------------------")
        print("100% Complete")
        # For Complementary Data on the final synchronny:
        print("\nYour final data set is on synchronny:", synchronnyList[-1])
        print("Decoder Frequency:", self.neuralNetwork.getFrequency("last"))
        
        
    def is_nmda_on(self):
        return self.listOfNetworks[0].getNeuron(0).isNMDA_on()
    def isGateDiagnosticsOn(self):
        return self.listOfNetworks[0].getNeuron(0).isGateDiagnosticsOn()
    def clearNetwork(self,name):
        self.listOfNetworks[0].fileWipe(name)
    def importMap(self,maP): # Only to be used for the single neuron thing and for CA3
        self.listOfNetworks[0].importAdjacencyMatrix(maP)
    
   
    def createDataLibrary(self):
        self.listOfTime = self.listOfNetworks[0].getTimeList()
        for x in range(0,len(self.listOfTime)):
            self.listOfTime[x] = float(self.listOfTime[x])
        
        # Total Field Potential
        listOfCA3Voltages = self.listOfNetworks[0].createDataLibrary()
        for dt in range(0,len(listOfCA3Voltages[0])):
            neuronVoltageOfASingleDtSum = 0
            for neuron in listOfCA3Voltages:
                neuronVoltageOfASingleDtSum += float(neuron[dt])
            self.totalFieldPotential.append(neuronVoltageOfASingleDtSum/len(listOfCA3Voltages))
            
        
        # Fourrier
        Fs = 1000/(float(self.listOfTime[1])-float(self.listOfTime[0]))
        #T = 1/Fs
        L = len(self.totalFieldPotential)
        #t = list(range(0,L))*T
        
        
        adjustedVoltages = []
        mean = sum(self.totalFieldPotential)/len(self.totalFieldPotential)
        for voltage in self.totalFieldPotential:
            adjustedVoltages.append(voltage-mean)
        
        rawFft = fft(adjustedVoltages)
        p2Fft = []
        p1Fft = []
        for element in rawFft:
            p2Fft.append(abs(element/L))
        
        p1Fft = p2Fft[0:(int(L/2)+1)]
        
        for x in range(1,len(p1Fft)):
            p1Fft[x] = 2*p1Fft[x]
        
        self.fourrierPower = p1Fft
        self.fourrierFrequency = list(range(0,int(L/2)+1))
        for x in range(0,len(self.fourrierFrequency)):
            self.fourrierFrequency[x] = (self.fourrierFrequency[x]* Fs)/L
         
            
        # List of Peak times
        heightFilter = -60
        timeFilter = 10
        self.smoothFieldPotential = (savgol_filter(self.totalFieldPotential, 61, 2)).tolist()
        
        for x in range(1,len(self.smoothFieldPotential)-1):
            # If statement checks that right and left values are less (ie its a peak), then checks that peak 
            # is tall enough, finally that its not too close to last peak
            if len(self.peakTimes) == 0:
                if self.smoothFieldPotential[x-1] <= self.smoothFieldPotential[x] and self.smoothFieldPotential[x+1] < self.smoothFieldPotential[x] and self.smoothFieldPotential[x] > heightFilter:
                    self.peakTimes.append(self.listOfTime[x])
            else:
                if self.smoothFieldPotential[x-1] <= self.smoothFieldPotential[x] and self.smoothFieldPotential[x+1] < self.smoothFieldPotential[x] and self.smoothFieldPotential[x] > heightFilter and (self.listOfTime[x] - self.peakTimes[-1]) > timeFilter:
                    self.peakTimes.append(self.listOfTime[x])                
                
        # Bin finder:
        binSize = 10
        maxHz = 200
        minHz = 0
        
        
        while minHz + binSize <= maxHz:
            
            numCycles = 0
            
            for x in range(0,len(self.peakTimes)-1):
                period = self.peakTimes[x+1]-self.peakTimes[x]
                frequency = 1000/period
                
                if frequency >= minHz and frequency < (minHz+binSize):
                    numCycles += 1
                    self.totalHzList.append(frequency)
            
            self.ca3ListOfCyclesPerHz.append(numCycles)
            minHz += 1
        
            
    def writeDataLibrary(self):
        
        
        # Erase the folder
        for file in os.listdir('CA3/DataLibrary'):
            os.remove('CA3/DataLibrary/' + file)
        
        # Write the list of total potential to the data file:
        fieldPotentialDataFile = open("CA3/DataLibrary/TotalFieldDataFile.txt", "w")
        fieldPotentialDataFile.write(('[%s]' % ', '.join(map(str, self.totalFieldPotential)))[1:-1])
        fieldPotentialDataFile.close()
        
        # Write list of smoothfield potential
        smoothPotentialDataFile = open("CA3/DataLibrary/SmoothFieldDataFile.txt", "w")
        smoothPotentialDataFile.write(('[%s]' % ', '.join(map(str, self.smoothFieldPotential)))[1:-1])
        smoothPotentialDataFile.close()
        
        # Write fourrier power vector
        fourrierPowerDataFile = open("CA3/DataLibrary/fourrierPowerDataFile.txt", "w")
        fourrierPowerDataFile.write(('[%s]' % ', '.join(map(str, self.fourrierPower)))[1:-1])
        fourrierPowerDataFile.close()   
        
        # Write fourrier frequency vector
        fourrierFrequencyDataFile = open("CA3/DataLibrary/fourrierFrequencyDataFile.txt", "w")
        fourrierFrequencyDataFile.write(('[%s]' % ', '.join(map(str, self.fourrierFrequency)))[1:-1])
        fourrierFrequencyDataFile.close() 
        
        # Write list of cycles per Hz
        histogramDataFile = open("CA3/DataLibrary/histogramDataFile.txt", "w")
        histogramDataFile.write(('[%s]' % ', '.join(map(str, self.ca3ListOfCyclesPerHz)))[1:-1])
        histogramDataFile.close() 
        
        # Write list of peak times
        peakDataFile = open("CA3/DataLibrary/peakDataFile.txt", "w")
        peakDataFile.write(('[%s]' % ', '.join(map(str, self.peakTimes)))[1:-1])
        peakDataFile.close()  
        
    def setNMDAparam(self,b_act_step, b_deact_step,p_act_step, p_deact_step, tau_nmda_act_1, tau_nmda_deact_1, tau_nmda_act_2, tau_nmda_deact_2):
        for network in self.listOfNetworks:
            network.setNMDAparam(b_act_step, b_deact_step,p_act_step, p_deact_step, tau_nmda_act_1, tau_nmda_deact_1, tau_nmda_act_2, tau_nmda_deact_2)
    def getTrialSD(self):
        return numpy.std(self.totalHzList)
    def getTrialFreq(self):
        return self.fourrierFrequency[self.fourrierPower.index(max(self.fourrierPower))]
if __name__ == '__main__':
    start_time = time.time()
    
    # Switches
    # Version 1: CA3 by itself
    # Version 2: CA3 and CA1
    # Version 3: A single neuron
    
    version = 3
    
    if version == 1:   # Model to use for CA3
        
        neuronNetwork = neuron_network_handler(False)
        
        # Preset for CA3
        
        ## Remember 0 for excitatory and 1 for inhibitory
        # Network ID
        networkName = 'CA3'
        
        # For the inhibition Neurons
        ## Spike Amplitude 1 = Amplitude sent to excitatory neurons
        ## Spike Amplitude 2 = Amplitude sent to inhibitory neurons
        ## Delay 1 = Delay to inhibitory neuron
        ## Delay 2 = Delay to excitation neuron        
        spikeAmplitude1     = .15         # .214 - table   More like excited conductance increase. Unit is: mS
        spikeAmplitude2     = .25         # .208 - table
        exciteTimeConstant  = .4          # .625 - prajal paper
        inhibTimeConstant   = 1           # .833 - table 
        designation         = 1           # Good
        numberOfNeurons     = 50          # Good
        delay1              = 1           # .6   - table
        delay2              = 1.4         # 1.1  - table
        wantBGCurrent       = True
        
        
        # For the excitatory Neurons
        ## Spike Amplitude 3 = Amplitude sent to excitatory neurons
        ## Spike Amplitude 4 = Amplitude sent to inhibitory neurons
        ## Delay 1 = Delay to inhibitory neuron
        ## Delay 2 = Delay to excitation neuron
        spikeAmplitude3     = .05           # .1   - prajal paper
        spikeAmplitude4     = .1            # .1   - prajal paper
        exciteTimeConstant2 = .5            # .588 - prajal paper    As time Constant goes up, the current falls faster
        inhibTimeConstant2  = .15           # .303 - table           As time Constant goes up, the current falls faster
        designation2        = 1             # Good
        numberOfNeurons2    = 200           # Good
        delay3              = 2.1           # 1.8  - prajal paper
        delay4              = 1.5           # .5   - prajal paper
        wantBGCurrent2      = True
        
        # For the decoder Neuron
        wantADecoder        = False
        spikeAmplitude5     = 1           ##
        spikeAmplitude6     = 1
        exciteTimeConstant3 = .909        # Takes Approx 47 ms to decay halfway
        inhibTimeConstant3  = .303        # Takes Approx 20 ms to decay halfway
        designation3        = 2           # Good
        delay5              = 0           # Good
        delay6              = 0           # Good
        
        
        # Neuron Network Initial Parameters
        duration         = 3000      # Good
        timeStep         = .02       # Good
        resolution       = .2        # Good
        
        
        # Fake Neuron Parameters
        ## Mean time and synchronny are used in the gauss distribution
        ## RangeMin and RangeMax are used during the uniform distribution draw
        ## Only mean time is used for the exponential distribution 
        ## Number of Spikes per Neuron is the ammount of times the excitatory neurons
        ##     will fire each. The ratios are kept even so one inhib cannot recieve 
        ##     3 spikes while another inhib recieves one spikes. Both would recive 
        ##     the same amount
        ## Excite Amplitude  = the current that the initial spikes send to an excitatory neuron
        ## Inhib Amplitude   = the current that the initial spikes send to an inhibitory neuron
        wantFakeNeurons = False
       
        exciteAmplitude  = .2            # Adjusted, is 1/3 of inhib amplitude
        inhibAmplitude   = .5             
        
        numberOfSpikesPerNeuron = 5
        meanTime                = 30
        synchronny              = 20
        RangeMin                = 1
        RangeMax                = 4
        typeOfDraw              = 1
        ##    1      2       3           4
        ##  Gauss  Expon  Uniform  Inverse Gauss
        
        
        # Synchronny Hz, Mean, Jitter Grapher UI
        '''
        Special Mapping = when on, the user must set the special mapping
        order to select one of the other special maps
        '''
        Want_To_Run          = False
        SyncMin              = 0
        SyncMax              = 20  
        Trials_Per_Sync      = 5      # Needed for Jitter Calculations
        HowManyNRasters      = 5
        
        #Custom Adjacency Matrix
        send_special_map = False
        #myMatrix = Adjacency_Matrix()
        #myMatrix.rapidAddNeuron(250)
        #myMatrix.presetThree(50,.25,.15,.25,.05)
        
        # CA3 -> CA1 Network Connection
        inhib_to_Inhib   = 0
        inhib_to_excite  = 0
        excite_to_excite = .50
        excite_to_inhib  = .50
        connect_to       = 'CA1'
        
        inhib_to_inhib_amp   = .65
        inhib_to_excite_amp  = .65
        excite_to_excite_amp = .01
        excite_to_inhib_amp  = .015         # Adjust this parameter in order to decrease inhib excitation?
        
        iNeuron_exciteTau = .625
        iNeuron_inhibTau  = .833
        eNeuron_exciteTau = .909
        eNeuron_inhibTau  = .303
        
        
        
        
        cp1 = [inhib_to_Inhib,inhib_to_excite,excite_to_excite,excite_to_inhib,connect_to]
        cp11 = ['CA1',inhib_to_inhib_amp,inhib_to_excite_amp,excite_to_inhib_amp,excite_to_excite_amp,iNeuron_inhibTau,iNeuron_exciteTau,eNeuron_inhibTau,eNeuron_exciteTau]
        
        neuronNetwork.createNetwork(         networkName,
                                             spikeAmplitude1,spikeAmplitude2,exciteTimeConstant,inhibTimeConstant,1,numberOfNeurons,delay1,delay2,wantBGCurrent,
                                               spikeAmplitude3,spikeAmplitude4,exciteTimeConstant2,inhibTimeConstant2,0,numberOfNeurons2,delay3,delay4,wantBGCurrent2,
                                               wantADecoder, spikeAmplitude5,spikeAmplitude6, exciteTimeConstant3,inhibTimeConstant3,2,delay5,delay6,
                                               duration,timeStep,resolution,
                                               wantFakeNeurons,exciteAmplitude,inhibAmplitude,numberOfSpikesPerNeuron,meanTime,synchronny,RangeMin,RangeMax,typeOfDraw,
                                               send_special_map)    
        
    
    
    
        
    elif version == 2: # Model for CA1 and CA3
        neuronNetwork = neuron_network_handler(True)
    
        # Preset for CA3
        
        ## Remember 0 for excitatory and 1 for inhibitory
        # Network ID
        networkName = 'CA3'
        
        # For the inhibition Neurons
        ## Spike Amplitude 1 = Amplitude sent to excitatory neurons
        ## Spike Amplitude 2 = Amplitude sent to inhibitory neurons
        ## Delay 1 = Delay to inhibitory neuron
        ## Delay 2 = Delay to excitation neuron        
        spikeAmplitude1     = .1         # More like excited conductance increase. Unit is: mS
        spikeAmplitude2     = .1
        exciteTimeConstant  = .625   #.303        #.625 - Abids       # Takes Approx 13 ms to decay halfway
        inhibTimeConstant   = .833        
        designation         = 1           # Good
        numberOfNeurons     = 50          # Good
        delay1              = .6     #1.1    #.6 - Abids   # Adjusted
        delay2              = .5     #.5 - Abids           # Adjusted
        wantBGCurrent       = True
        
        # For the excitatory Neurons
        ## Spike Amplitude 3 = Amplitude sent to excitatory neurons
        ## Spike Amplitude 4 = Amplitude sent to inhibitory neurons
        ## Delay 1 = Delay to inhibitory neuron
        ## Delay 2 = Delay to excitation neuron
        spikeAmplitude3     = .1           # Is 3x stronger than EPSP to P cell
        spikeAmplitude4     = .1
        exciteTimeConstant2 = .909 #.588     #.909 - Abids      # Takes Approx 13 ms to decay halfway
        inhibTimeConstant2  = .303 #.625     #.303       ## As of now, there is no inhib to inhib, so this doesn't do anything
        designation2        = 1          # Good
        numberOfNeurons2    = 200        # Good
        delay3              = 1.1 #1.8    #1.1 - Abids
        delay4              = 1.8 #.5     #1.8 - Abids   # Adjusted 
        wantBGCurrent2      = True
        
        # For the decoder Neuron
        wantADecoder        = False
        spikeAmplitude5     = 1           ##
        spikeAmplitude6     = 1
        exciteTimeConstant3 = .909        # Takes Approx 47 ms to decay halfway
        inhibTimeConstant3  = .303        # Takes Approx 20 ms to decay halfway
        designation3        = 2           # Good
        delay5              = 0           # Good
        delay6              = 0           # Good
        
        
        # Neuron Network Initial Parameters
        duration         = 600         # Good
        timeStep         = .04            # Good
        resolution       = .2            # Good
        
        
        # Fake Neuron Parameters
        ## Mean time and synchronny are used in the gauss distribution
        ## RangeMin and RangeMax are used during the uniform distribution draw
        ## Only mean time is used for the exponential distribution 
        ## Number of Spikes per Neuron is the ammount of times the excitatory neurons
        ##     will fire each. The ratios are kept even so one inhib cannot recieve 
        ##     3 spikes while another inhib recieves one spikes. Both would recive 
        ##     the same amount
        ## Excite Amplitude  = the current that the initial spikes send to an excitatory neuron
        ## Inhib Amplitude   = the current that the initial spikes send to an inhibitory neuron
        wantFakeNeurons = False
       
        exciteAmplitude  = .2            # Adjusted, is 1/3 of inhib amplitude
        inhibAmplitude   = .5             
        
        numberOfSpikesPerNeuron = 5
        meanTime                = 30
        synchronny              = 20
        RangeMin                = 1
        RangeMax                = 4
        typeOfDraw              = 1
        ##    1      2       3           4
        ##  Gauss  Expon  Uniform  Inverse Gauss
        
        
        # Synchronny Hz, Mean, Jitter Grapher UI
        '''
        Special Mapping = when on, the user must set the special mapping
        order to select one of the other special maps
        '''
        Want_To_Run          = False
        SyncMin              = 0
        SyncMax              = 20  
        Trials_Per_Sync      = 5      # Needed for Jitter Calculations
        HowManyNRasters      = 5
        
        #Custom Adjacency Matrix
        send_special_map = False
        #myMatrix = Adjacency_Matrix()
        #myMatrix.rapidAddNeuron(250)
        #myMatrix.presetThree(50,.25,.15,.25,.05)
        
        # CA3 -> CA1 Network Connection
        inhib_to_Inhib   = 0
        inhib_to_excite  = 0
        excite_to_excite = .50
        excite_to_inhib  = .50
        connect_to       = 'CA1'
        
        inhib_to_inhib_amp   = .65
        inhib_to_excite_amp  = .65
        excite_to_excite_amp = .01
        excite_to_inhib_amp  = .015         # Adjust this parameter in order to decrease inhib excitation?
        
        iNeuron_exciteTau = .625
        iNeuron_inhibTau  = .833
        eNeuron_exciteTau = .909
        eNeuron_inhibTau  = .303
        
        
        
        
        cp1 = [inhib_to_Inhib,inhib_to_excite,excite_to_excite,excite_to_inhib,connect_to]
        cp11 = ['CA1',inhib_to_inhib_amp,inhib_to_excite_amp,excite_to_inhib_amp,excite_to_excite_amp,iNeuron_inhibTau,iNeuron_exciteTau,eNeuron_inhibTau,eNeuron_exciteTau]
        
        neuronNetwork.createNetwork(         networkName,
                                             spikeAmplitude1,spikeAmplitude2,exciteTimeConstant,inhibTimeConstant,1,numberOfNeurons,delay1,delay2,wantBGCurrent,
                                               spikeAmplitude3,spikeAmplitude4,exciteTimeConstant2,inhibTimeConstant2,0,numberOfNeurons2,delay3,delay4,wantBGCurrent2,
                                               wantADecoder, spikeAmplitude5,spikeAmplitude6, exciteTimeConstant3,inhibTimeConstant3,2,delay5,delay6,
                                               duration,timeStep,resolution,
                                               wantFakeNeurons,exciteAmplitude,inhibAmplitude,numberOfSpikesPerNeuron,meanTime,synchronny,RangeMin,RangeMax,typeOfDraw,
                                               send_special_map)    
        
        
        
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    
        
        # Preset for CA1
        
        ## Remember 0 for excitatory and 1 for inhibitory
        # Network ID
        networkName2 = 'CA1'
        
        # For the inhibition Neurons
        ## Spike Amplitude 1 = Amplitude sent to excitatory neurons
        ## Spike Amplitude 2 = Amplitude sent to inhibitory neurons
        ## Delay 1 = Delay to inhibitory neuron
        ## Delay 2 = Delay to excitation neuron        
        spikeAmplitude1     = .85         # Is 3x stronger than EPSP to P cell
        spikeAmplitude2     = .65         ## Irrelevant
        exciteTimeConstant  = .285        # ?????????
        inhibTimeConstant   = 0           ## Irrelevant
        designation         = 1           # Good
        numberOfNeurons     = 50          # Good
        delay1              = 0           ## Irrelevant
        delay2              = .9          # Adjusted
        wantBGCurrent       = True
        
        # For the excitatory Neurons
        ## Spike Amplitude 3 = Amplitude sent to excitatory neurons
        ## Spike Amplitude 4 = Amplitude sent to inhibitory neurons
        ## Delay 1 = Delay to inhibitory neuron
        ## Delay 2 = Delay to excitation neuron
        spikeAmplitude3     = .1           # Is 3x stronger than EPSP to P cell
        spikeAmplitude4     = .1
        exciteTimeConstant2 = .588       # Takes Approx 13 ms to decay halfway
        inhibTimeConstant2  = .286       ## As of now, there is no inhib to inhib, so this doesn't do anything
        designation2        = 1          # Good
        numberOfNeurons2    = 200        # Good
        delay3              = 1.1
        delay4              = 1.8        # Adjusted 
        wantBGCurrent2      = True
        
        # For the decoder Neuron
        wantADecoder        = False
        spikeAmplitude5     = 1           ##
        spikeAmplitude6     = 1
        exciteTimeConstant3 = .909        # Takes Approx 47 ms to decay halfway
        inhibTimeConstant3  = .303        # Takes Approx 20 ms to decay halfway
        designation3        = 2           # Good
        delay5              = 0           # Good
        delay6              = 0           # Good
        
        
        # Neuron Network Initial Parameters
        #duration         = 200            # Good
        #timeStep         = .1             # Good
        #resolution       = .2             # Good
        
        
        # Fake Neuron Parameters
        ## Mean time and synchronny are used in the gauss distribution
        ## RangeMin and RangeMax are used during the uniform distribution draw
        ## Only mean time is used for the exponential distribution 
        ## Number of Spikes per Neuron is the ammount of times the excitatory neurons
        ##     will fire each. The ratios are kept even so one inhib cannot recieve 
        ##     3 spikes while another inhib recieves one spikes. Both would recive 
        ##     the same amount
        ## Excite Amplitude  = the current that the initial spikes send to an excitatory neuron
        ## Inhib Amplitude   = the current that the initial spikes send to an inhibitory neuron
        wantFakeNeurons = False
        
        exciteAmplitude  = .2            # Adjusted, is 1/3 of inhib amplitude
        inhibAmplitude   = .5             
        
        numberOfSpikesPerNeuron = 5
        meanTime                = 30
        synchronny              = 20
        RangeMin                = 1
        RangeMax                = 4
        typeOfDraw              = 1
        ##    1      2       3           4
        ##  Gauss  Expon  Uniform  Inverse Gauss
        
        
        # Synchronny Hz, Mean, Jitter Grapher UI
        '''
        Special Mapping = when on, the user must set the special mapping
        order to select one of the other special maps
        '''
        Want_To_Run          = False
        SyncMin              = 0
        SyncMax              = 20  
        Trials_Per_Sync      = 5      # Needed for Jitter Calculations
        HowManyNRasters      = 5
        
        #Custom Adjacency Matrix
        send_special_map = False
        #myMatrix = Adjacency_Matrix()
        #myMatrix.rapidAddNeuron(250)
        #myMatrix.presetThree(50,.25,.15,.25,.05)    
      
        # CA1 -> CA3 Network Connection
        inhib_to_Inhib   = 0
        inhib_to_excite  = 0
        excite_to_excite = 0
        excite_to_inhib  = 0
        connect_to       = 'CA3'    
        
        # Amps that CA1 sends to CA3
        inhib_to_inhib_amp   = .65
        inhib_to_excite_amp  = .65
        excite_to_excite_amp = .1
        excite_to_inhib_amp  = .1
        
        # Tau that CA1 will use for current sent from CA3
        iNeuron_exciteTau = .588             # decrease even further to reduce inhibitory current?
        iNeuron_inhibTau  = 0
        eNeuron_exciteTau = .588
        eNeuron_inhibTau  = .286 
            
            
        
        cp2 = [inhib_to_Inhib,inhib_to_excite,excite_to_excite,excite_to_inhib,connect_to]
        cp22 = ['CA3',inhib_to_inhib_amp,inhib_to_excite_amp,excite_to_inhib_amp,excite_to_excite_amp,iNeuron_inhibTau,iNeuron_exciteTau,eNeuron_inhibTau,eNeuron_exciteTau]
        
        if send_special_map == True:
            send_special_map = myMatrix
     
    
        neuronNetwork.createNetwork(         networkName2,
                                               spikeAmplitude1,spikeAmplitude2,exciteTimeConstant,inhibTimeConstant,1,numberOfNeurons,delay1,delay2,wantBGCurrent,
                                               spikeAmplitude3,spikeAmplitude4,exciteTimeConstant2,inhibTimeConstant2,0,numberOfNeurons2,delay3,delay4,wantBGCurrent2,
                                               wantADecoder, spikeAmplitude5,spikeAmplitude6, exciteTimeConstant3,inhibTimeConstant3,2,delay5,delay6,
                                               duration,timeStep,resolution,
                                               wantFakeNeurons,exciteAmplitude,inhibAmplitude,numberOfSpikesPerNeuron,meanTime,synchronny,RangeMin,RangeMax,typeOfDraw,
                                               send_special_map)
    
        neuronNetwork.updateInfoDictionary(networkName,[cp1],[cp11])
        
        neuronNetwork.updateInfoDictionary(networkName2,[cp2],[cp22])        
    
    elif version == 3: # Model for one neuron
        neuronNetwork = neuron_network_handler(False)
    
        # Preset for CA3
    
        ## Remember 0 for excitatory and 1 for inhibitory
        # Network ID
        networkName = 'CA3'
    
        # For the inhibition Neurons
        ## Spike Amplitude 1 = Amplitude sent to excitatory neurons
        ## Spike Amplitude 2 = Amplitude sent to inhibitory neurons
        ## Delay 1 = Delay to inhibitory neuron
        ## Delay 2 = Delay to excitation neuron        
        spikeAmplitude1     = .1         # More like excited conductance increase. Unit is: mS
        spikeAmplitude2     = .1
        exciteTimeConstant  = .625   #.303        #.625 - Abids       # Takes Approx 13 ms to decay halfway
        inhibTimeConstant   = .833        
        designation         = 1           # Good
        numberOfNeurons     = 0          # Good
        delay1              = .6     #1.1    #.6 - Abids   # Adjusted
        delay2              = .5     #.5 - Abids           # Adjusted
        wantBGCurrent       = True
    
        # For the excitatory Neurons
        ## Spike Amplitude 3 = Amplitude sent to excitatory neurons
        ## Spike Amplitude 4 = Amplitude sent to inhibitory neurons
        ## Delay 1 = Delay to inhibitory neuron
        ## Delay 2 = Delay to excitation neuron
        spikeAmplitude3     = .1           # Is 3x stronger than EPSP to P cell
        spikeAmplitude4     = .1
        exciteTimeConstant2 = .909 #.588     #.909 - Abids      # Takes Approx 13 ms to decay halfway
        inhibTimeConstant2  = .303 #.625     #.303       ## As of now, there is no inhib to inhib, so this doesn't do anything
        designation2        = 1          # Good
        numberOfNeurons2    = 1        # Good
        delay3              = 1.1 #1.8    #1.1 - Abids
        delay4              = 1.8 #.5     #1.8 - Abids   # Adjusted 
        wantBGCurrent2      = True
    
        # For the decoder Neuron
        wantADecoder        = False
        spikeAmplitude5     = 1           ##
        spikeAmplitude6     = 1
        exciteTimeConstant3 = .909        # Takes Approx 47 ms to decay halfway
        inhibTimeConstant3  = .303        # Takes Approx 20 ms to decay halfway
        designation3        = 2           # Good
        delay5              = 0           # Good
        delay6              = 0           # Good
    
    
        # Neuron Network Initial Parameters
        duration         = 300         # Good
        timeStep         = .01            # Good
        resolution       = .2            # Good
    
    
        # Fake Neuron Parameters
        ## Mean time and synchronny are used in the gauss distribution
        ## RangeMin and RangeMax are used during the uniform distribution draw
        ## Only mean time is used for the exponential distribution 
        ## Number of Spikes per Neuron is the ammount of times the excitatory neurons
        ##     will fire each. The ratios are kept even so one inhib cannot recieve 
        ##     3 spikes while another inhib recieves one spikes. Both would recive 
        ##     the same amount
        ## Excite Amplitude  = the current that the initial spikes send to an excitatory neuron
        ## Inhib Amplitude   = the current that the initial spikes send to an inhibitory neuron
        wantFakeNeurons = False
    
        exciteAmplitude  = .2            # Adjusted, is 1/3 of inhib amplitude
        inhibAmplitude   = .5             
    
        numberOfSpikesPerNeuron = 5
        meanTime                = 30
        synchronny              = 20
        RangeMin                = 1
        RangeMax                = 4
        typeOfDraw              = 1
        ##    1      2       3           4
        ##  Gauss  Expon  Uniform  Inverse Gauss
    
    
        # Synchronny Hz, Mean, Jitter Grapher UI
        '''
            Special Mapping = when on, the user must set the special mapping
            order to select one of the other special maps
            '''
        Want_To_Run          = False
        SyncMin              = 0
        SyncMax              = 20  
        Trials_Per_Sync      = 5      # Needed for Jitter Calculations
        HowManyNRasters      = 5
    
        #Custom Adjacency Matrix
        send_special_map = True
        myMatrix = Adjacency_Matrix()
        myMatrix.rapidAddNeuron(1)
        #myMatrix.presetThree(50,.25,.15,.25,.05)
    
        # CA3 -> CA1 Network Connection
        inhib_to_Inhib   = 0
        inhib_to_excite  = 0
        excite_to_excite = .50
        excite_to_inhib  = .50
        connect_to       = 'CA1'
    
        inhib_to_inhib_amp   = .65
        inhib_to_excite_amp  = .65
        excite_to_excite_amp = .01
        excite_to_inhib_amp  = .015         # Adjust this parameter in order to decrease inhib excitation?
    
        iNeuron_exciteTau = .625
        iNeuron_inhibTau  = .833
        eNeuron_exciteTau = .909
        eNeuron_inhibTau  = .303
    
    
    
    
        cp1 = [inhib_to_Inhib,inhib_to_excite,excite_to_excite,excite_to_inhib,connect_to]
        cp11 = ['CA1',inhib_to_inhib_amp,inhib_to_excite_amp,excite_to_inhib_amp,excite_to_excite_amp,iNeuron_inhibTau,iNeuron_exciteTau,eNeuron_inhibTau,eNeuron_exciteTau]
    
        neuronNetwork.createNetwork(         networkName,
                                                 spikeAmplitude1,spikeAmplitude2,exciteTimeConstant,inhibTimeConstant,1,numberOfNeurons,delay1,delay2,wantBGCurrent,
                                                 spikeAmplitude3,spikeAmplitude4,exciteTimeConstant2,inhibTimeConstant2,0,numberOfNeurons2,delay3,delay4,wantBGCurrent2,
                                                   wantADecoder, spikeAmplitude5,spikeAmplitude6, exciteTimeConstant3,inhibTimeConstant3,2,delay5,delay6,
                                                   duration,timeStep,resolution,
                                                   wantFakeNeurons,exciteAmplitude,inhibAmplitude,numberOfSpikesPerNeuron,meanTime,synchronny,RangeMin,RangeMax,typeOfDraw,
                                                   send_special_map)        
        
    else:
        raise ValueError('Unrecognized version number')
    
    # Prints out the Adjacency Matrices for the user to see
    neuronNetwork.printMap()
    
    # Want to run is for the synchronny trials. As of this version, want to run is unavailable.
    # Special Map allows the user to create their own map. See the adjacency matrix class to help
    # create your own matrix. This function is still available.
    
    if Want_To_Run == False:
        if send_special_map == True:
            neuronNetwork.importMap(myMatrix)
            neuronNetwork.normalRun()
        else:
            neuronNetwork.normalRun()
    else:
        if send_special_map == True:
            if myMatrix.requestPreset() == "Three" and wantADecoder == True:
                raise ValueError("You cant have a decoder with map preset 3")
            if myMatrix.requestPreset() == "Three" and wantFakeNeurons == True:
                raise ValueError("You cant have fake neurons with map preset 3")
            neuronNetwork.sendSpecialMap(myMatrix)
            neuronNetwork.syncRun(SyncMin,SyncMax, Trials_Per_Sync,HowManyNRasters,True)             
        else:
            neuronNetwork.syncRun(SyncMin,SyncMax, Trials_Per_Sync,HowManyNRasters,False)
    
    # Prints some basic information.
    print("This simulation took: ", (time.time() - start_time), 'seconds')
    if time.time() - start_time > 60:
        print("Which is also equal to: ", (time.time() - start_time)/60, 'minutes')
        if time.time() - start_time > 3600:
            print("Which is also equal to: ", (time.time() - start_time)/3600, 'hours')
    if version == 1:
        CA3_ON = True
        CA1_ON = False
    elif version == 2:
        CA3_ON = True
        CA1_ON = True
    else:
        CA3_ON = False
        CA1_ON = False        
    
        
    print('\nSome Basic Run Stats:\n-----------------------------\nDuration: '+ str(duration) +'\nCA3 On? : '+str(CA3_ON) +'\nCA1 On? : '+str(CA1_ON) +'\nNMDA On?: '+ str(neuronNetwork.is_nmda_on()))
    
    
    