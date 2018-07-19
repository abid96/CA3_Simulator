import random

class Adjacency_Matrix:
    
    
    def __init__(self,multiNetwork = False):
        '''
        # Basic initializer. It starts the object off with an empty list in a list,
        # and a counter for the number of neurons. The counter starts at 0
        
        Rows are the lists
        Columns are the elements of the lists: 
        
        From and To Neurons are only applicable if nxn is not the same
        '''
        self.doubleMatrix = [[]]
        self.numberOfNeurons = 0
        self.numberOfFromNeurons = 0
        self.numberOfToNeurons = 0
        self.presetStyle = "None Set"
        self.multiNetwork = multiNetwork
    
    
    def addNeuron(self):
        '''
        # Basicaly adds a neuron by adding a section to the interior list, and adds
        # a list to the double matrix of the proper size to make it a symmetrical
        # matrix
        '''
        
        self.numberOfNeurons += 1
        for row in self.doubleMatrix:
            row.append(0)
        if self.numberOfNeurons != 1:
            self.doubleMatrix.append([0]*self.numberOfNeurons)
    
   
    def rapidAddNeuron(self,howMany):
        '''
        # Add 'n' neurons at a time by calling 'addNeuron()' n times
        '''
        
        while howMany != 0:
            self.addNeuron()
            howMany -= 1
    
    
    def drawConnection(self,neuronFrom, to):
        '''
        # Makes a 'connection' by putting '1' for the adjacency matrix
        '''
        
        if neuronFrom == to and self.multiNetwork == False:
            raise ValueError('A Neuron Cannot Connect to Itself!')
        self.doubleMatrix[neuronFrom-1][to-1] = 1
    
   
    def deleteConnection(self,neuronFrom, to):
        '''
         # doesn't really work yet. 
        '''
        self.doubleMatrix[neuronFrom-1][to-1] = 0
        
    
    def drawMultiConnections(self,neuronFrom,listOfNeuronsTo):
        '''
        # provide the neuron from, and a list of neurons to in order
        # to draw multiple connections
        '''
        
        for to in listOfNeuronsTo:
            self.drawConnection(neuronFrom,to)
    def rapidConnection(self,neuronFrom,small,big):
        for x in range(small,big+1):
            if neuronFrom != x:
                self.drawConnection(neuronFrom,x)
    
    def deleteMultiConnections(self,neuronFrom,listOfNeuronsTo):
        '''
        # still experimental
        '''
        
        for to in listOfNeuronsTo:
            self.deleteConnection(neuronFrom,to)        
    
    
    def getConnectionIDs(self,neuronID):
        '''
        # will return a list of neurons that recieve input from a given neuron (neuronID)
        '''
        listOfConnectingNeurons = []
        if self.multiNetwork == True:
            for x in range(self.numberOfToNeurons):
                
                if self.doubleMatrix[neuronID-1][x] == 1:
                    listOfConnectingNeurons.append(x+1)
        else:
            for x in range(self.numberOfNeurons):
                
                if self.doubleMatrix[neuronID-1][x] == 1:
                    listOfConnectingNeurons.append(x+1)
        
        return listOfConnectingNeurons
    
    
    def __len__(self):
        '''
        # returns number of neurons
        '''
        return self.numberOfNeurons
    
    
    def plainText(self,matrixName):
        '''
        # it will essentially churn out code for how to make the exact given matrix
        # useful for when one wants to random generate a matrix, but keep that random
        # generated one
        '''
        instructionString = matrixName+'.rapidAddNeuron('+ str(self.numberOfNeurons)+')\n'
        for x in range(1,self.numberOfNeurons):
            instructionString += matrixName + '.drawMultiConnections(' + str(x) + ',' + str(self.getConnectionIDs(x))+')\n'
        return instructionString
            
    
    
    def __str__(self):
        '''
        # Prints out the matrix in a form that is (relatively) highly readable
        '''
        totalString = "   "
        lineString = "  -"
        rowNumber = 1
        printReference = 0
        if self.numberOfNeurons == 0:
            if self.numberOfToNeurons != 0:
                printReference = self.numberOfToNeurons
        else:
            printReference = self.numberOfNeurons
        
        for x in range(1,printReference+1):
            if x < 10:
                totalString += "0" + str(x) + " "
            else:
                totalString += str(x) + " "
            lineString += "---"
        
        totalString += "\n" + lineString + "\n"
        
        for row in self.doubleMatrix:
            if rowNumber < 10:
                totalString += "0" + str(rowNumber) + "| " 
            else: 
                totalString += str(rowNumber) + "| "
            
            for column in row:
                totalString += str(column) + "  "
            
            totalString += "\n"
            rowNumber += 1
        
        return totalString
    
    
    
    def requestPreset(self):
        
        '''
        If a preset was requested on the matrix, 
        it will return that preset name as "One" or "Two", etc.
        If no preset was requested, then it will return "None"
        '''
        
        return self.presetStyle
    
    def presetOne(self):
        '''
        # All inhibitory neurons are connected to the decoder
        '''
        self.presetStyle = "One"
        
        for fromNeuron in range(1,self.numberOfNeurons):
            self.drawConnection(fromNeuron,self.numberOfNeurons)
    
    
    
    def presetTwo(self):
        '''
        # All inhibitory neurons are connected to the decoder and
        # each inhibitory neuron has a 50% chance of connecting to
        # other inhibitory neurons
        '''
        self.presetStyle = "Two"
        
        for fromNeuron in range(1, self.numberOfNeurons):
            self.drawConnection(fromNeuron,self.numberOfNeurons)
            
            for fromNeuronsConnections in range(1, self.numberOfNeurons):
                if random.getrandbits(1) == 1 and fromNeuron != fromNeuronsConnections:
                    self.drawConnection(fromNeuron,fromNeuronsConnections)
    
    def presetThree(self,howManyInterneurons,II,IE,EI,EE):
        '''
        This map has a decoder. Its a combination of
        '''
        self.presetStyle = "Three"
        
        for fromNeuron in range(1, howManyInterneurons+1):                         # This selects the from interneurons
             
            for fromNeuronsConnections in range(1,howManyInterneurons+1):          # This selects the connections of the interneuron that are inhibitory
                if random.random() < II and fromNeuron != fromNeuronsConnections:  # If the chance is good and the connection isn't itself
                    self.drawConnection(fromNeuron,fromNeuronsConnections)
            
            for fromNeuronsConnections in range(howManyInterneurons+1,self.numberOfNeurons+1):  # This selects the connections of the interneuron that are excitatory
                if random.random() < IE:                                                        # No need to check if it's itself. Interneuron cant be an excitatory neuron
                    self.drawConnection(fromNeuron,fromNeuronsConnections)
        
        for fromNeuron in range(howManyInterneurons+1,self.numberOfNeurons+1):     # This selects the from excitatory neurons
            
            #self.drawConnection(fromNeuron,self.numberOfNeurons)
            
            for fromNeuronsConnections in range(1,howManyInterneurons+1):          # This selects the connections of the excitatory neuron that are inhibitory
                if random.random() < EI:
                    self.drawConnection(fromNeuron,fromNeuronsConnections)
            
            for fromNeuronsConnections in range(howManyInterneurons+1,self.numberOfNeurons+1):  # This selects the connections of the excitatory neuron that are excitatory
                if random.random() < EE and fromNeuron != fromNeuronsConnections:
                    self.drawConnection(fromNeuron,fromNeuronsConnections)
    def presetFour(self,II,IE,EI,EE,howManyFromINs,howManyToINs):
        '''
        This will be used for the multi network map
        '''
        self.presetStyle = 'Four'
        for fromNeuron in range(1, howManyFromINs+1):                         # This selects the from interneurons
            
            for fromNeuronsConnections in range(1,howManyToINs+1):          # This selects the connections of the interneuron that are inhibitory
                if random.random() < II :                                   # If the chance is good
                    self.drawConnection(fromNeuron,fromNeuronsConnections)            
            
            for fromNeuronsConnections in range(howManyToINs+1,self.numberOfToNeurons+1):  # This selects the connections of the interneuron that are excitatory
                if random.random() < IE:                                                   
                    self.drawConnection(fromNeuron,fromNeuronsConnections) 
        
        for fromNeuron in range(howManyFromINs+1,self.numberOfFromNeurons+1):     # This selects the from excitatory neurons
            
            for fromNeuronsConnections in range(1,howManyToINs+1):         # This selects the connections of the excitatory neuron that are inhibitory
                if random.random() < EI:
                    self.drawConnection(fromNeuron,fromNeuronsConnections)
            
            for fromNeuronsConnections in range(howManyToINs+1,self.numberOfToNeurons+1):  # This selects the connections of the excitatory neuron that are excitatory
                if random.random() < EE:
                    self.drawConnection(fromNeuron,fromNeuronsConnections)            
        
    def generate(self,row,column):
        '''
        Used to generate an NxN sized matrix. Useful if one already
        knows the needed size or uses an uneven matrix. 
        '''
        
        self.doubleMatrix = [None] * row
        for i in range(row):
            self.doubleMatrix[i] = [0] * column        
        
        
        if self.multiNetwork == False:
            if row == column:
                self.numberOfNeurons = row
            else: 
                raise ValueError("Map size mismatch. You can't have different sizes between rows and columns for non multi-network")
        else:
            self.numberOfFromNeurons = row
            self.numberOfToNeurons = column
    def getReverseConnection(self,idNumber):
        listOfReverseCons = []
        for x in range(0,len(self.doubleMatrix)):
            if self.doubleMatrix[x][idNumber-1] == 1:
                listOfReverseCons.append(x+1)
        
        return listOfReverseCons
           
if __name__ == '__main__':
    myMatrix = Adjacency_Matrix(True)
    
    myMatrix.generate(20,15)
    myMatrix.presetFour(.25,.25,.25,.25,10,6)
    x = myMatrix.getReverseConnection(4)
    print(myMatrix)
    print(x)
    
    