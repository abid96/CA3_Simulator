from collections import deque
class Neuron:
    
    import math
    import os
    import random
    
    
    def __init__(self,spikeAmplitude,spikeAmplitude2, exciteTimeConstant,inhibTimeConstant,designation,idNumber,delay1,delay2,wantBGCurrent,networkID):
        '''
        # Note on ID Number: The program provides the neuron ID numbers as 1,2,3,4... etc, but the stored ID numbers
        # start from 0. In other words, what the user knows as neuron 1 has "0" stored as its ID. Another example, 
        # running getID() on the neuron with ID = 0 will return 1.
        # Am I decoder is used because there is special data we want to store for the decoder only. That way we can store the special
        # data for that neuron only. It is false by default.
        
        # Excite Spike Amplitude = Amplitude to send to excitatory neurons
        # Inhib Spike Amplitude  = Amplitude to send to inhibitory neurons
        # Delay 1 = delay to I
        # Delay 2 = delay to E
        '''
        
        self.exciteSpikeAmplitude = spikeAmplitude  
        self.inhibSpikeAmplitude  = spikeAmplitude2 
        self.exciteTimeConstant   = exciteTimeConstant
        self.inhibTimeConstant    = inhibTimeConstant
        self.designation          = designation      # 0 for excitatory, 1 for inhibitory, 2 decoder
        self.idNumber             = idNumber
        self.delay1               = delay1
        self.delay2               = delay2
        if self.designation == 2:
            self.wantBGCurrent    = False
            if wantBGCurrent == True:
                print("Warning: You can't have background current for a decoder.\nThe simulation will run with no background current for the decoder")
        else:
            self.wantBGCurrent    = wantBGCurrent
        self.duration = 0
        self.sampleRate = 0
        self.exciteReversal = 0.0  #mV
        self.inhibReversal = -70.0 #mV
        self.excitatoryConductance = 0
        self.inhibitoryConductance = 0
        self.excitatoryCurrent = 0
        self.inhibitoryCurrent = 0
        self.backGroundCurrent = 0
        self.previousVoltage = 0
        
        paramSet = 1
        
        if paramSet == 1:
            self.voltage = -65
            self.NaG = 120    # mS/cm^2
            self.KG  = 36     # mS/cm^2
            self.LG  = 0.3    # mS/cm^2
            self.Vna = 50     # mV
            self.VK  = -77    # mV
            self.VL  = -65  # mV 
            self.m   = 0.0    # 0-1 scale
            self.h   = 0.0     # 0-1 scale
            self.n   = 0.0     # 0-1 scale
            self.cap = 1   #uF/cm^2
        elif paramSet == 2:
            self.voltage = 0
            self.NaG = 120    # mS/cm^2
            self.KG  = 36     # mS/cm^2
            self.LG  = 0.3    # mS/cm^2
            self.Vna = 115     # mV
            self.VK  = -12    # mV
            self.VL  = 10.6  # mV 
            self.m   = 0.1    # 0-1 scale
            self.h   = 0     # 0-1 scale
            self.n   = .9     # 0-1 scale
            self.cap = 1   #uF/cm^2            
        
        self.totalTime = 0.0
        self.dv = 0
        self.delayStepsForInhib  = 0           
        self.delayStepsForExcite = 0            
        self.resolutionSteps = 0
        self.resolutionCounter = 0
        self.numberOfSpikes = 0
        self.sendSpikeQueue = []
        self.recievedSpikes = deque()
        self.cleanTime = 20
        self.networkID = networkID
        self.listOfConnections = []
        
        
        # Program Switches
        self.nmdaON = True
        self.nmdaDiagnostics = False
        self.gateDiagnostics = False
        self.rungKutta = False
        
        # NMDA Experimental Junk
        self.b_act_step, self.b_deact_step = 0.001,0.001           
        self.p_act_step, self.p_deact_step = 0.01, 0.01 
        
        self.gVd = 0.0
        self._nmda_act, self._nmda_deact = 0.0,0.0
        self.nmdaCurrent = 0.0
        self.E_nmda = 0 #137.04
        
        if self.designation == 0:   # excite
            self.tau_nmda_act = 2.2340 #2.2340
            self.tau_nmda_deact = 365 #62.5808
        else:                       # inhib
            self.tau_nmda_act = 2.2340
            self.tau_nmda_deact = 500
       
       
        self.listofVoltages            = []
        self.listofInhibCurrent        = []
        self.listofExciteCurrent       = []
        self.listofTime                = []
        self.listofSpikeTimes          = []
        self.listOfbackGroundCurrent   = []
        self.listOfNMDACurrent         = []
        
        # Will only apply if NMDA diagnostics is on
        self.listOfNMDAInfo            = [] # will be used to transfer info between methods
        self.listOfNMDAConductance     = []
        self.listOfMgBlock             = []
        self.listOfdactdt              = []
        self.listOfddeactdt            = []
        self.listOfselfAct             = []
        self.listOfselfDeact           = []
        
        # Will only apply if gate diagnostics is on
        self.listOfm                   = []  # Note: List of m,n,h have been changed to hold the full conductances. H holds leak conductance, M hold Na conductance,
        self.listOfh                   = []  #       and n holds potassium conductance
        self.listOfn                   = []
        self.listOfdv                  = []
        self.dvPt1                     = []
        self.dvPt2                     = []
        self.dvPt3                     = []
        self.dvPt4                     = []
        
        
        self.idOfInterest = 0
        
        if self.gateDiagnostics == True and self.nmdaDiagnostics == True:
            raise ValueError('You cannot have NMDA and Gate Diagnostics at the same time')

    
        
    def getNetworkID(self):
        return self.networkID
    def setConnections(self,listOfConnections):
        self.listOfConnections = listOfConnections

    def getFirstSpikeTime(self):
        if len(self.listofSpikeTimes) == 0:
            return -1
        else:
            return self.listofSpikeTimes[0]

    def getSpikeTimes(self):
        return self.listofSpikeTimes

    def getDesignation(self):
        return self.designation
    
    def getID(self):
        return self.idNumber+1
    
    def reset(self):
        ''' 
        this method is out of date. It will not do the job properly
        '''
        
        self.voltage = 0.0
        self.excitatoryConductance = 0
        self.inhibitoryConductance = 0
        self.excitatoryCurrent = 0
        self.inhibitoryCurrent = 0        
        self.refract = False
        self.refractTimer = 0
        self.totalTime = 0.0
        self.excitatoryConductanceDecayTime = 0.0
        self.inhibitoryConductanceDecayTime = 0.0
        self.dv = 0
        self.delayTimer = 0
        self.resolutionCounter = 0
        self.numberOfSpikes = 0
        self.sendSpikeQueue = []
        self.recievedSpikes = deque()
       
        self.listofVoltages      = []
        self.listofInhibCurrent  = []
        self.listofExciteCurrent = []
        self.listofTime          = []
        self.listofSpikeTimes    = []        
    
    
    
    def getFrequency(self):
        return (self.numberOfSpikes/self.totalTime)*1000
       
       
    def getNumberOfSpikes(self):
        return self.numberOfSpikes
    
    
    def spikeMailHandler(self,mailbox,outNetworkSpikeParameters):
        '''
        # spikesInQueue is basically a queue to send AP's to other neurons. Its in a queue because we need
        # to introduce a latency rather than just sending the AP right away. Spike queue looks like:
        # [(Delay time, neuron object), (Delay time, neuron object),...]
        
        This method basically looks through the send spike que, and will send out the mail when the timer is 0.
        It will then delete the spike from the queue. If the time isn't right, then it just decrements the timer
        
        This method will return the updated mailbox
        '''
        for spikesInQueue in self.sendSpikeQueue:
            if spikesInQueue[0] == 0:
                # Differentiates between what type of neuron we are sending
                # current to in order to pick the correct amplitude
                if spikesInQueue[1].getDesignation() == 0:
                    if spikesInQueue[1].getNetworkID() == self.networkID: # neuron is of the same network:
                        mailbox.append([self.exciteSpikeAmplitude,self.designation, self.totalTime,spikesInQueue[1],self.networkID,True])
                    else:  # is of a different network
                        if self.designation == 0:
                            mailbox.append([outNetworkSpikeParameters[spikesInQueue[1].getNetworkID()][3],self.designation, self.totalTime,spikesInQueue[1],self.networkID,True])
                        else:
                            mailbox.append([outNetworkSpikeParameters[spikesInQueue[1].getNetworkID()][1],self.designation, self.totalTime,spikesInQueue[1],self.networkID,True])
                            
                        
                else:
                    if spikesInQueue[1].getNetworkID() == self.networkID: # neuron is of the same network:
                        mailbox.append([self.inhibSpikeAmplitude,self.designation, self.totalTime,spikesInQueue[1],self.networkID,True])
                    else: # is of a different network
                        if self.designation == 0:
                            mailbox.append([outNetworkSpikeParameters[spikesInQueue[1].getNetworkID()][2],self.designation, self.totalTime,spikesInQueue[1],self.networkID,True])
                        else:
                            mailbox.append([outNetworkSpikeParameters[spikesInQueue[1].getNetworkID()][0],self.designation, self.totalTime,spikesInQueue[1],self.networkID,True])
                        
                self.sendSpikeQueue.remove(spikesInQueue)
            else:
                spikesInQueue[0] -= 1        
        
    def calculateExternalCurrent(self,outNetworkSpikeParameters):
        '''
        This method is used to calculate the inhibitory, excitatory, background, 
        and NMDA current (if applicable) for the neuron. It returns the final
        result as the sum of all these currents. This method also cleans out
        the recieved spike data as necessary.
        
        How it works: first it cleans the recieved spike structure by iterating 
        through and decrementing what needs to be decremented, and deleting spikes
        who have reached the time limit defined by the user.
        
        Next, it increments through the recieved spikes and calculates the
        conductance changes for inhib and excite, and it figures other vars
        needed for NMDA. 
        
        Using the stuff calculated in the above step, it will calculate 
        NMDA, inhib, and excite current. It will then randomly generate
        the background current.
        
        Finally it returns the sum of all the currents. 
        '''
        
        # Cleaning out recieved spike list
        while len(self.recievedSpikes) != 0 and self.totalTime - self.recievedSpikes[0][2] > self.cleanTime:
            self.recievedSpikes.popleft()        
        
        # Setup vars for calculating excitatory, inhibitory, and NMDA currents
        self.excitatoryConductance = 0
        self.inhibitoryConductance = 0
            
        
        # Calculating setup vars from recieved spikes
        for recievedSpike in self.recievedSpikes:  # Recieved spike structure: [amplitude, designation, spike time,neuron object,networkID,True]
                                                   # True is used for NMDA calculation. Basically if True, then the NMDA calc will go through, but
                                                   # then it will be switched to false, and 
            
            # NMDA Var
            if self.nmdaON == True:
                # NMDA Addition. We only add once per spike, then not again
                if recievedSpike[1] == 0 and recievedSpike[5] == True:
                    if self.designation == 0:
                        self._nmda_act += self.p_act_step
                        self._nmda_deact += self.p_deact_step
                    else:
                        self._nmda_act += self.b_act_step
                        self._nmda_deact += self.b_deact_step
                    recievedSpike[5] = False

            
            # Excitatory Var
            if recievedSpike[1] == 0:   
                if recievedSpike[4] == self.networkID: # Of the same network
                    self.excitatoryConductance += recievedSpike[0]*self.math.exp(-(self.totalTime-recievedSpike[2])*self.exciteTimeConstant)
                else: # different network
                    if self.designation == 0:
                        self.excitatoryConductance += recievedSpike[0]*self.math.exp(-(self.totalTime-recievedSpike[2])*outNetworkSpikeParameters[recievedSpike[4]][7])
                        
                    else:
                        self.excitatoryConductance += recievedSpike[0]*self.math.exp(-(self.totalTime-recievedSpike[2])*outNetworkSpikeParameters[recievedSpike[4]][5])
                    
            # Inhibitory Var       
            else:                      
                if recievedSpike[4] == self.networkID: # Of the same network
                    self.inhibitoryConductance += recievedSpike[0]*self.math.exp(-(self.totalTime-recievedSpike[2])*self.inhibTimeConstant)
                else: # different network
                    if self.designation == 0:
                        self.inhibitoryConductance += recievedSpike[0]*self.math.exp(-(self.totalTime-recievedSpike[2])*outNetworkSpikeParameters[recievedSpike[3]][6])
                    else:
                        self.inhibitoryConductance += recievedSpike[0]*self.math.exp(-(self.totalTime-recievedSpike[2])*outNetworkSpikeParameters[recievedSpike[3]][4])
                    
            
        
        # Calculating NMDA Current
        if self.nmdaON == True:
            self.calcNMDA()
        else:
            self.nmdaCurrent = 0.0            
        
        # Calculating Excitatory Current
        self.excitatoryCurrent = -self.excitatoryConductance*(self.voltage-self.exciteReversal) # gives micro amps/cm^2
        
        # Calculating Inhibitory Current
        self.inhibitoryCurrent = -self.inhibitoryConductance*(self.voltage-self.inhibReversal)  # gives micro amps/cm^2
        

        #Calculating background current (only for non decoders)
        # MAKE SURE ITS A DECIMAL!
        if self.wantBGCurrent == True: # unit of micro amp
            #self.backGroundCurrent = 18*self.random.random()
            if self.designation == 0:
                if self.numberOfSpikes == 0:
                    self.backGroundCurrent = 20*self.random.random()#
                else:
                    self.backGroundCurrent = 20*self.random.random()# change back to 18 *
            elif self.designation == 1:
                if self.numberOfSpikes == 0:
                    self.backGroundCurrent = 15*self.random.random()
                else:
                    self.backGroundCurrent = 15*self.random.random()
            
        # Final result:
        return (self.excitatoryCurrent + self.inhibitoryCurrent + self.backGroundCurrent + self.nmdaCurrent) # Micro amps/cm^2
    
    def calcNMDA(self):
        # Additional factors for NMDA tuning
        # Taken from "A fast model of voltage-dependent NMDA receptors" - moradi et al
        
        # These parameters are for the voltage dependent tau
        #if self.idNumber == 0 and self._nmda_act != 0:
            #print('hello')
        a_a = 0.0896   # ms
        a_b = 10.0374  # ms
        lam_a = 0.0324 # 1/mV
        lam_b = 0.0239 # 1/mV
        
        mgo = 1    # 1
        ic50 = 4.1 # IC50 at 0 mV
        d = .8     # the electrical distance of the Mg2+ binding site from the outside of the membrane
        F = 96480  # coul
        R = 8.315  # J/degC 
        T = 273.16 # ˚C
        z = 2      # Mg+2 valence
        gVi = 1    # µS Conductance of Voltage Independent component
        tau_g = 7  #
        k = .007   #
        v_0 = -100 #
       
        dactdt = -self._nmda_act / (self.tau_nmda_act+a_a*self.math.exp(-lam_a*self.voltage))
        ddeactdt = -self._nmda_deact / (self.tau_nmda_deact+a_b*(1-self.math.exp(-lam_b*self.voltage)))
        gVd_inf = k*(self.voltage-v_0)
        d_gVd = (gVd_inf - self.gVd)/tau_g
        mgblock = 1/(1+ mgo*(1/ic50)*(self.math.exp(0.001*-z*d*F*self.voltage*(1/R)*(1/T))))
        nmdaConductance = gVi + self.gVd
        
        self.gVd += (d_gVd * self.sampleRate)
        self._nmda_act += (dactdt * self.sampleRate)        
        self._nmda_deact +=  (ddeactdt * self.sampleRate)
        
        
        self.nmdaCurrent = -(self._nmda_deact-self._nmda_act) * nmdaConductance * mgblock * (self.voltage - self.E_nmda)
        
        if self.nmdaDiagnostics == True:
            self.listOfNMDAInfo = []
            self.listOfNMDAInfo.append(nmdaConductance)
            self.listOfNMDAInfo.append(mgblock)
            self.listOfNMDAInfo.append(dactdt)
            self.listOfNMDAInfo.append(ddeactdt)
            
    
    def rungKuttaGate(self, a_n, b_n, a_m, b_m, a_h, b_h):
        '''
        This method is a more accurate way of finding d(gateVar)/dt.
        It will set up m, n, h. Its up for the cycle method to use those
        to actually increment the voltage for the neuron. 
        
        While more accurate, its slower
        '''
        n1 = self.n
        m1 = self.m
        h1 = self.h
        n2 = self.n
        m2 = self.m
        h2 = self.h  
        n3 = self.n
        m3 = self.m
        h3 = self.h              

        
        # Derivative 1 at t(0*)
        dn1 = a_n*(1-self.n) - b_n*self.n
        dm1 = a_m*(1-self.m) - b_m*self.m
        dh1 = a_h*(1-self.h) - b_h*self.h
        
        # Value 1 at y(step/2) using derivative 1
        n1 += dn1*(self.sampleRate/2)
        m1 += dm1*(self.sampleRate/2)
        h1 += dh1*(self.sampleRate/2)
        
        # Derivative 2 at t(step/2) using value 1
        dn2 = a_n*(1-n1) - b_n*n1
        dm2 = a_m*(1-m1) - b_m*m1
        dh2 = a_h*(1-h1) - b_h*h1            
        
        # Value 2 at y(step/2) using derivative 2
        n2 += dn2*(self.sampleRate/2)
        m2 += dm2*(self.sampleRate/2)
        h2 += dh2*(self.sampleRate/2)
        
        # Derivative 3 at t(step/2) using value 2
        dn3 = a_n*(1-n2) - b_n*n2
        dm3 = a_m*(1-m2) - b_m*m2
        dh3 = a_h*(1-h2) - b_h*h2            
        
        # Value 3 at y(step) using derivative 3
        n3 += dn3*(self.sampleRate)
        m3 += dm3*(self.sampleRate)
        h3 += dh3*(self.sampleRate)            
        
        # Derivative 4 at t(step) using value 3
        dn4 = a_n*(1-n3) - b_n*n3
        dm4 = a_m*(1-m3) - b_m*m3
        dh4 = a_h*(1-h3) - b_h*h3            
        
        self.n += (dn1 + 2*dn2 + 2*dn3 + dn4)/6 * self.sampleRate
        self.m += (dm1 + 2*dm2 + 2*dm3 + dm4)/6 * self.sampleRate
        self.h += (dh1 + 2*dh2 + 2*dh3 + dh4)/6 * self.sampleRate
            
    def eulerGate(self,a_n, b_n, a_m, b_m, a_h, b_h):
        '''
        This method is a less accurate way of finding d(gateVar)/dt.
        It will set up m, n, h. Its up for the cycle method to use those
        to actually increment the voltage for the neuron.
        
        While less accurate, its faster
        '''
        dn = a_n*(1-self.n) - b_n*self.n
        dm = a_m*(1-self.m) - b_m*self.m
        dh = a_h*(1-self.h) - b_h*self.h
    
        self.m += dm*self.sampleRate
        self.n += dn*self.sampleRate
        self.h += dh*self.sampleRate
        
    def cycle(self,mailbox,outNetworkSpikeParameters):
        '''
        # This method essentially means to run the neuron through a single time step
        # It will maintain the basics such as decaying current, checking maintaining
        # euler and what not. 
        
        Spike collector (mailbox) will be provided by the neuron network script. It is
        simply an empty list. Every neuron will then dump their spikes into the 
        collective for the neuron network. on the next time step, the neuron
        network will redistribute it to the right place. It essentially a mail service
        and the spikecollector is your mailbox, and neuron network will pick up from there.
        this structure is utilized so that no neuron experiences a time delay due to
        order of processing.
        '''
        
        
        # This step handles the spike queue and mailing out spikes
        self.spikeMailHandler(mailbox,outNetworkSpikeParameters)
        
        # This step calculates the external current
        externalCurrent = self.calculateExternalCurrent(outNetworkSpikeParameters)
        
        
        # Choose your HH eqn
        hh_mode = 3
        if hh_mode == 1:   # HH eqns from: http://neuronaldynamics.epfl.ch/online/Ch2.S2.html
            #if self.voltage > 0:
                #print('hello')
            
            a_n = (0.02*(self.voltage-25)/(1-self.math.exp(-(self.voltage-25)/9)))   
            b_n = (-0.002*(self.voltage-25)/(1-self.math.exp((self.voltage-25)/9)))
        
            a_m = (0.182*(self.voltage+35)/(1-self.math.exp(-(self.voltage+35)/9)))
            b_m = (-0.124*(self.voltage+35)/(1-self.math.exp((self.voltage+35)/9)))
        
            a_h = 0.25*self.math.exp(-(self.voltage+90)/12)
            b_h = 0.25*(self.math.exp((self.voltage+62)/6)/self.math.exp((self.voltage+90)/12))        
        
        elif hh_mode == 2: # HH eqns from: http://icwww.epfl.ch/~gerstner/SPNM/node14.html
            
            a_n = ((0.1-0.01*self.voltage)/(self.math.exp(1-0.01*self.voltage)-1))
            b_n = 0.125*self.math.exp(-self.voltage/80)
            
            a_m = ((2.5-0.1*self.voltage)/(self.math.exp(2.5-0.1*self.voltage)-1))
            b_m = 4*self.math.exp(-self.voltage/18)
            
            a_h = 0.07*self.math.exp(-self.voltage/20)
            b_h = 1/(self.math.exp(3-0.1*self.voltage)+1)
        elif hh_mode == 3: # HH eqns from: 
            a_m = 0.1*(self.voltage+40.0)/(1.0 - self.math.exp(-(self.voltage+40.0) / 10.0))
            b_m = 4.0*self.math.exp(-(self.voltage+65.0) / 18.0)
            
            a_n = 0.01*(self.voltage+55.0)/(1.0 - self.math.exp(-(self.voltage+55.0) / 10.0))
            b_n = 0.125*self.math.exp(-(self.voltage+65) / 80.0)
            
            a_h = 0.07*self.math.exp(-(self.voltage+65.0) / 20.0)
            b_h = 1.0/(1.0 + self.math.exp(-(self.voltage+35.0) / 10.0))
        else:
            raise ValueError('Undefined HH Mode')
        
        # Here we are only incrementing m, n, and h only after the first step
        # because the increment should use the initial values for the firts
        # step
        if self.rungKutta == True:
            if self.totalTime != 0.0:
                self.rungKuttaGate(a_n, b_n, a_m, b_m, a_h, b_h)
        else:
            if self.totalTime != 0.0:
                self.eulerGate(a_n, b_n, a_m, b_m, a_h, b_h)
                
            
            
        dvPt1 = self.NaG * self.math.pow(self.m,3) * self.h * (self.voltage - self.Vna) # micro amps/cm^2
        dvPt2 = self.KG * self.math.pow(self.n,4) * (self.voltage - self.VK)
        dvPt3 = self.LG * (self.voltage - self.VL)
        dvPt4 = externalCurrent

        self.dv = (-dvPt1 - dvPt2 - dvPt3 + dvPt4)/(self.cap)
        self.previousVoltage = self.voltage
        self.voltage += self.dv*self.sampleRate
        

        if self.voltage >= 20 and  self.previousVoltage < 20:
            #if self.idNumber == 0:
                #print('hello') 
            
            if self.idOfInterest == self.idNumber and self.idOfInterest != 0:
                print('Hello')
        
            if len(self.listofSpikeTimes) != 0 and self.totalTime - float(self.listofSpikeTimes[-1]) < .03:
                if self.idOfInterest == 0:
                    self.idOfInterest = self.idNumber
                print('hello') 
            
            self.listofSpikeTimes.append("{:.6}".format(self.totalTime))
            self.numberOfSpikes += 1

            for connections in self.listOfConnections:
                
                if connections.getDesignation() == 0:
                    self.sendSpikeQueue.append([self.delayStepsForExcite,connections])
                else:
                    self.sendSpikeQueue.append([self.delayStepsForInhib,connections])        
        
        
        # The data resolution check. Basically will make sure
        # the data written to files is of the correct resolution
        if self.resolutionCounter % self.resolutionSteps == 0:
            self.listofVoltages.append("{:.6}".format(self.voltage))
            self.listofInhibCurrent.append("{:.6}".format(self.inhibitoryCurrent))
            self.listofExciteCurrent.append("{:.6}".format(self.excitatoryCurrent))
            self.listofTime.append("{:.6}".format(self.totalTime))
            self.listOfNMDACurrent.append("{:.6}".format(self.nmdaCurrent))
            if self.designation != 2 and self.wantBGCurrent == True:
                self.listOfbackGroundCurrent.append("{:.6}".format(self.backGroundCurrent))
            
            if self.nmdaDiagnostics == True:
                self.listOfNMDAConductance.append("{:.6}".format(self.listOfNMDAInfo[0]))
                self.listOfMgBlock.append("{:.6}".format(self.listOfNMDAInfo[1]))
                self.listOfdactdt.append("{:.6}".format(self.listOfNMDAInfo[2]))            
                self.listOfddeactdt.append("{:.6}".format(self.listOfNMDAInfo[3]))          
                self.listOfselfAct.append("{:.6}".format(self._nmda_act))   
                self.listOfselfDeact.append("{:.6}".format(self._nmda_deact))
            
            if self.gateDiagnostics == True:
                self.listOfm.append("{:.6}".format(self.NaG * self.math.pow(self.m,3) * self.h))
                self.listOfh.append("{:.6}".format(self.LG))
                self.listOfn.append("{:.6}".format(self.KG * self.math.pow(self.n,4)))
                self.listOfdv.append("{:.6}".format(self.dv))
                self.dvPt1.append("{:.6}".format(dvPt1))
                self.dvPt2.append("{:.6}".format(dvPt2))
                self.dvPt3.append("{:.6}".format(dvPt3))
                self.dvPt4.append("{:.6}".format(dvPt4))
            
        
        # Basic, Necessary Incrementation
        self.resolutionCounter += 1
        self.totalTime += self.sampleRate
        
       
    
    
    
    def sendConductance(self,amplitude,fromNeuronDesignation,spikeTime):
        
        '''
        # This method will essentially add a provided amplitude to the internal
        # conductance value, and set the countdown timer to 0
        '''
        
        #if fromNeuronDesignation == 0 or fromNeuronDesignation == 2:
            #self.excitatoryConductance += amplitude
            #self.excitatoryConductanceDecayTime = 0
        #else:
            #self.inhibitoryConductance += amplitude
            #self.inhibitoryConductanceDecayTime = 0
    
        self.recievedSpikes.append([amplitude,fromNeuronDesignation,spikeTime])
    
    
    
    def recieveSpike(self,spikeInfo):
        '''
        This method is the way that the neuron network program will
        send the 'mail' back to the appropriate neurons
        '''
        self.recievedSpikes.append(spikeInfo)
    
    def setTime(self, timeStep, duration,resolution):
        '''
        # This method will be used in the higher class to set a universal time
        # and sample rate
        '''
        
        self.sampleRate = timeStep
        self.duration = duration
        self.resolutionSteps = int(resolution/timeStep)
        self.delayStepsForInhib = round(float(self.delay1)/float(self.sampleRate))
        self.delayStepsForExcite = round(float(self.delay2)/float(self.sampleRate))
        
       
        
    
    
   
    def write(self):
        '''
        # This method will first clear the folders once. Then it will create all the 
        # necessary files and write the data to them
        '''
        
        if (len(self.listofSpikeTimes) == 0):
            self.listofSpikeTimes.append(-1)
        
        # Designation for excitatory
        if self.designation == 0:
            voltagedataFile = open(self.networkID+"/ExcitatoryNeuron/Voltage/voltage" + str(self.idNumber+1) + ".txt", "w")
            inhibCurrentdataFile = open(self.networkID+"/ExcitatoryNeuron/InhibitoryCurrent/inhibCurrent" + str(self.idNumber+1) + ".txt", "w")
            exciteCurrentdataFile = open(self.networkID+"/ExcitatoryNeuron/ExcitatoryCurrent/exciteCurrent" + str(self.idNumber+1) + ".txt", "w")
            spikeTimeDataFile = open(self.networkID+"/ExcitatoryNeuron/SpikeTime/spikeTime" + str(self.idNumber+1) + ".txt", "w")
            bgCurrentDataFile = open(self.networkID+"/ExcitatoryNeuron/BackgroundCurrent/bg" + str(self.idNumber+1) + ".txt", "w")
            nmdaCurrentDataFile = open(self.networkID+"/ExcitatoryNeuron/NMDA/nmda" + str(self.idNumber+1) + ".txt", "w")
            if self.nmdaDiagnostics == True:
                nmdaCondcutanceDataFile = open(self.networkID+"/ExcitatoryNeuron/NMDA_Diagnostics/NMDA_Cond/nmda" + str(self.idNumber+1) + ".txt", "w")
                mgBlockDataFile = open(self.networkID+"/ExcitatoryNeuron/NMDA_Diagnostics/mgBlock/nmda" + str(self.idNumber+1) + ".txt", "w")
                dactdtDataFile = open(self.networkID+"/ExcitatoryNeuron/NMDA_Diagnostics/dactdt/nmda" + str(self.idNumber+1) + ".txt", "w")
                ddeactdtDataFile = open(self.networkID+"/ExcitatoryNeuron/NMDA_Diagnostics/ddeactdt/nmda" + str(self.idNumber+1) + ".txt", "w")
                selfDeactDataFile = open(self.networkID+"/ExcitatoryNeuron/NMDA_Diagnostics/self.deact/nmda" + str(self.idNumber+1) + ".txt", "w")
                selfActDataFile = open(self.networkID+"/ExcitatoryNeuron/NMDA_Diagnostics/self.act/nmda" + str(self.idNumber+1) + ".txt", "w")
            if self.gateDiagnostics == True:
                hDataFile   = open(self.networkID+"/ExcitatoryNeuron/gate_diagnostics/h/gate" + str(self.idNumber+1) + ".txt", "w")
                mDataFile   = open(self.networkID+"/ExcitatoryNeuron/gate_diagnostics/m/gate" + str(self.idNumber+1) + ".txt", "w")
                nDataFile   = open(self.networkID+"/ExcitatoryNeuron/gate_diagnostics/n/gate" + str(self.idNumber+1) + ".txt", "w")
                dvDataFile  = open(self.networkID+"/ExcitatoryNeuron/gate_diagnostics/dv/gate" + str(self.idNumber+1) + ".txt", "w")
                dv1DataFile = open(self.networkID+"/ExcitatoryNeuron/gate_diagnostics/dv1/gate" + str(self.idNumber+1) + ".txt", "w")
                dv2DataFile = open(self.networkID+"/ExcitatoryNeuron/gate_diagnostics/dv2/gate" + str(self.idNumber+1) + ".txt", "w")
                dv3DataFile = open(self.networkID+"/ExcitatoryNeuron/gate_diagnostics/dv3/gate" + str(self.idNumber+1) + ".txt", "w")
                dv4DataFile = open(self.networkID+"/ExcitatoryNeuron/gate_diagnostics/dv4/gate" + str(self.idNumber+1) + ".txt", "w")
                
        # Designation for inhibitory
        elif self.designation == 1:
            voltagedataFile = open(self.networkID+"/InhibitoryNeuron/Voltage/voltage" + str(self.idNumber+1) + ".txt", "w")
            inhibCurrentdataFile = open(self.networkID+"/InhibitoryNeuron/InhibitoryCurrent/inhibCurrent" + str(self.idNumber+1) + ".txt", "w")
            exciteCurrentdataFile = open(self.networkID+"/InhibitoryNeuron/ExcitatoryCurrent/exciteCurrent" + str(self.idNumber+1) + ".txt", "w")
            spikeTimeDataFile = open(self.networkID+"/InhibitoryNeuron/SpikeTime/spikeTime" + str(self.idNumber+1) + ".txt", "w")
            bgCurrentDataFile = open(self.networkID+"/InhibitoryNeuron/BackgroundCurrent/bg" + str(self.idNumber+1) + ".txt", "w")
            nmdaCurrentDataFile = open(self.networkID+"/InhibitoryNeuron/NMDA/nmda" + str(self.idNumber+1) + ".txt", "w")
            if self.nmdaDiagnostics == True:
                nmdaCondcutanceDataFile = open(self.networkID+"/InhibitoryNeuron/NMDA_Diagnostics/NMDA_Cond/nmda" + str(self.idNumber+1) + ".txt", "w")
                mgBlockDataFile = open(self.networkID+"/InhibitoryNeuron/NMDA_Diagnostics/mgBlock/nmda" + str(self.idNumber+1) + ".txt", "w")
                dactdtDataFile = open(self.networkID+"/InhibitoryNeuron/NMDA_Diagnostics/dactdt/nmda" + str(self.idNumber+1) + ".txt", "w")
                ddeactdtDataFile = open(self.networkID+"/InhibitoryNeuron/NMDA_Diagnostics/ddeactdt/nmda" + str(self.idNumber+1) + ".txt", "w")
                selfDeactDataFile = open(self.networkID+"/InhibitoryNeuron/NMDA_Diagnostics/self.deact/nmda" + str(self.idNumber+1) + ".txt", "w")
                selfActDataFile = open(self.networkID+"/InhibitoryNeuron/NMDA_Diagnostics/self.act/nmda" + str(self.idNumber+1) + ".txt", "w")  
            if self.gateDiagnostics == True:
                hDataFile  = open(self.networkID+"/InhibitoryNeuron/gate_diagnostics/h/gate" + str(self.idNumber+1) + ".txt", "w")
                mDataFile  = open(self.networkID+"/InhibitoryNeuron/gate_diagnostics/m/gate" + str(self.idNumber+1) + ".txt", "w")
                nDataFile  = open(self.networkID+"/InhibitoryNeuron/gate_diagnostics/n/gate" + str(self.idNumber+1) + ".txt", "w")
                dvDataFile = open(self.networkID+"/InhibitoryNeuron/gate_diagnostics/dv/gate" + str(self.idNumber+1) + ".txt", "w")
                dv1DataFile = open(self.networkID+"/InhibitoryNeuron/gate_diagnostics/dv1/gate" + str(self.idNumber+1) + ".txt", "w")
                dv2DataFile = open(self.networkID+"/InhibitoryNeuron/gate_diagnostics/dv2/gate" + str(self.idNumber+1) + ".txt", "w")
                dv3DataFile = open(self.networkID+"/InhibitoryNeuron/gate_diagnostics/dv3/gate" + str(self.idNumber+1) + ".txt", "w")
                dv4DataFile = open(self.networkID+"/InhibitoryNeuron/gate_diagnostics/dv4/gate" + str(self.idNumber+1) + ".txt", "w")                
        # Designation for decoder
        else:
            voltagedataFile = open(self.networkID+"/DecoderNeuron/Voltage/voltage" + str(self.idNumber+1) + ".txt", "w")
            inhibCurrentdataFile = open(self.networkID+"/DecoderNeuron/InhibitoryCurrent/inhibCurrent" + str(self.idNumber+1) + ".txt", "w")
            exciteCurrentdataFile = open(self.networkID+"/DecoderNeuron/ExcitatoryCurrent/exciteCurrent" + str(self.idNumber+1) + ".txt", "w")
            spikeTimeDataFile = open(self.networkID+"/DecoderNeuron/SpikeTime/spikeTime" + str(self.idNumber+1) + ".txt", "w")
        
        
        voltagedataFile.write(('[%s]' % ', '.join(map(str, self.listofVoltages)))[1:-1])
        inhibCurrentdataFile.write(('[%s]' % ', '.join(map(str, self.listofInhibCurrent)))[1:-1])
        exciteCurrentdataFile.write(('[%s]' % ', '.join(map(str, self.listofExciteCurrent)))[1:-1])
        spikeTimeDataFile.write(('[%s]' % ', '.join(map(str, self.listofSpikeTimes)))[1:-1])
        nmdaCurrentDataFile.write(('[%s]' % ', '.join(map(str, self.listOfNMDACurrent)))[1:-1])
        if self.nmdaDiagnostics == True:
            nmdaCondcutanceDataFile.write(('[%s]' % ', '.join(map(str, self.listOfNMDAConductance)))[1:-1])
            mgBlockDataFile.write(('[%s]' % ', '.join(map(str, self.listOfMgBlock)))[1:-1])
            dactdtDataFile.write(('[%s]' % ', '.join(map(str, self.listOfdactdt)))[1:-1])
            ddeactdtDataFile.write(('[%s]' % ', '.join(map(str, self.listOfddeactdt)))[1:-1])
            selfDeactDataFile.write(('[%s]' % ', '.join(map(str, self.listOfselfDeact)))[1:-1])
            selfActDataFile.write(('[%s]' % ', '.join(map(str, self.listOfselfAct)))[1:-1])           
        if self.designation == 1 or self.designation == 0:
            bgCurrentDataFile.write(('[%s]' % ', '.join(map(str, self.listOfbackGroundCurrent)))[1:-1])
        if self.gateDiagnostics == True:
            hDataFile.write(('[%s]' % ', '.join(map(str, self.listOfh)))[1:-1])
            mDataFile.write(('[%s]' % ', '.join(map(str, self.listOfm)))[1:-1])
            nDataFile.write(('[%s]' % ', '.join(map(str, self.listOfn)))[1:-1])
            dvDataFile.write(('[%s]' % ', '.join(map(str, self.listOfdv)))[1:-1])
            dv1DataFile.write(('[%s]' % ', '.join(map(str, self.dvPt1)))[1:-1])
            dv2DataFile.write(('[%s]' % ', '.join(map(str, self.dvPt2)))[1:-1])
            dv3DataFile.write(('[%s]' % ', '.join(map(str, self.dvPt3)))[1:-1])
            dv4DataFile.write(('[%s]' % ', '.join(map(str, self.dvPt4)))[1:-1])
            
        
        voltagedataFile.close()
        inhibCurrentdataFile.close()
        exciteCurrentdataFile.close()
        spikeTimeDataFile.close()
        nmdaCurrentDataFile.close()
        if self.nmdaDiagnostics == True:
            nmdaCondcutanceDataFile.close()
            mgBlockDataFile.close()
            dactdtDataFile.close()
            ddeactdtDataFile.close()
            selfDeactDataFile.close()
            selfActDataFile.close()        
        if self.designation == 1 or self.designation == 0:
            bgCurrentDataFile.close()
        if self.gateDiagnostics == True:
            hDataFile.close()
            mDataFile.close()
            nDataFile.close()
            dvDataFile.close()
            dv1DataFile.close()
            dv2DataFile.close()
            dv3DataFile.close()
            dv4DataFile.close()
        
        # We only need to do this one once since all the time data files are the exact same
        if (self.idNumber == 0):
            timedataFile = open(self.networkID+"/Time/time" + str(self.idNumber+1) + ".txt", "w")
            timedataFile.write(('[%s]' % ', '.join(map(str, self.listofTime)))[1:-1])        
            timedataFile.close()
        
        
    
        
    def isNMDA_on(self):
        return self.nmdaON
    def isGateDiagnosticsOn(self):
        return self.gateDiagnostics
    def getVoltageList(self):
        return self.listofVoltages
    def getTimeList(self):
        return self.listofTime
    def setNMDAparam(self,b_act_step, b_deact_step,p_act_step, p_deact_step, tau_nmda_act_1, tau_nmda_deact_1, tau_nmda_act_2, tau_nmda_deact_2):
        self.b_act_step, self.b_deact_step = b_act_step,b_deact_step          
        self.p_act_step, self.p_deact_step = p_act_step, p_deact_step
        
        if self.designation == 0:   # excite
            self.tau_nmda_act = tau_nmda_act_1 
            self.tau_nmda_deact = tau_nmda_deact_1
        else:                       # inhib
            self.tau_nmda_act = tau_nmda_act_2
            self.tau_nmda_deact = tau_nmda_deact_2     
        
        pass
if __name__ == '__main__':
    
    #spikeAmplitude = 2
    #exciteTimeConstant = 1
    #inhibTimeConstant = 1
    #threshold = 1
    #refractDuration = 2
    #designation = 0
    #idNumber = 0
    #delay = 0
    
    #spikeAmplitude2 = 2
    #exciteTimeConstant2 = 1
    #inhibTimeConstant2 = 1
    #threshold2 = 1
    #refractDuration2 = 2
    #designation2 = 0
    #idNumber2 = 1
    #delay2 = 0
    
    
    #a = Neuron(spikeAmplitude,exciteTimeConstant,inhibTimeConstant,threshold,refractDuration,designation,idNumber,delay)
    #b = Neuron(spikeAmplitude2,exciteTimeConstant2,inhibTimeConstant2,threshold2,refractDuration2,designation2,idNumber2,delay2)
    
    #a.setTime(.1,10,.1)
    #b.setTime(.1,10,.1)
    
    #asConnections = [b]
    #bsConnections = [a]
    
    #steps = 0
    #while (steps != 100):
        #a.cycle(asConnections)
        #b.cycle(bsConnections)
        
        #if steps%10 == 0:
            #a.sendConductance(2)
        #steps += 1
    
    #a.write()
    #b.write()
    
    print("\n Wrong File!!")