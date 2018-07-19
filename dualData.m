if exist('multiRun','var') == false
    multiRun = false;
end

if multiRun == false
    x = input('WARNING: Are you sure you want to import data for both neuron networks? (y/n) ');

    if x == 'y'
        disp('Okay, running')
    else
        disp('Process Aborted')
        return 
    end
end 

close all
if multiRun == true
    clearvars -except multiRunInfo how_Many_Runs_Do_You_Have multiRun numberOfRuns run1 run2
else
    clearvars -except multiRun
end


%%%% Switches
CA1_ON = false;
NMDA_Diagnostics = false;
Gate_Diagnostics = false;
Auto_guide = true;  % Allows program to smart detect best option for switches
dataIgnore = 1;

%%%% Some preparation

if Auto_guide == true
    if length(dir('NeuronNetworkClass/CA1/InhibitoryNeuron/Voltage'))-2 == 0
        CA1_ON = false;
    elseif multiRun == true
        CA1_ON = false;
    else
        CA1_ON = true;
    end
    
    if length(dir('NeuronNetworkClass/CA3/InhibitoryNeuron/NMDA_Diagnostics/NMDA_Cond'))-2 ~= 0
        NMDA_Diagnostics = true;
    end
    
    if length(dir('NeuronNetworkClass/CA3/InhibitoryNeuron/gate_diagnostics/m'))-2 ~= 0
        Gate_Diagnostics = true;
    end
end

% Gathering the time data file
timeData = csvread('NeuronNetworkClass/CA3/Time/time1.txt');

if Gate_Diagnostics == true
    dataLength = length(timeData)-1;
else
    dataLength = length(timeData);
end

timeData = timeData(dataIgnore:dataLength);
totalTimeinMs = timeData(length(timeData));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gather Inhibitory Neuron Data
ca3_numberOfInhibitoryNeuronFiles = length(dir('NeuronNetworkClass/CA3/InhibitoryNeuron/Voltage'))-2;

ca3_inhibitoryNeuronVoltageData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
for k = 1: ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/Voltage/voltage%d.txt', k);
    ca3_inhibitoryNeuronVoltageData{k} = csvread(myfilename);
    ca3_inhibitoryNeuronVoltageData{k} = ca3_inhibitoryNeuronVoltageData{k}(dataIgnore:dataLength);
end

ca3_inhibitoryNeuronSpikeTimeData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
for k = 1: ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/SpikeTime/spikeTime%d.txt', k);
    ca3_inhibitoryNeuronSpikeTimeData{k} = csvread(myfilename);
    ca3_inhibitoryNeuronSpikeTimeData{k} = ca3_inhibitoryNeuronSpikeTimeData{k}(ca3_inhibitoryNeuronSpikeTimeData{k}>timeData(1));
end

ca3_inhibitoryNeuronInhibitoryCurrentData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
for k = 1: ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/InhibitoryCurrent/inhibCurrent%d.txt', k);
    ca3_inhibitoryNeuronInhibitoryCurrentData{k} = csvread(myfilename);
    ca3_inhibitoryNeuronInhibitoryCurrentData{k}= ca3_inhibitoryNeuronInhibitoryCurrentData{k}(dataIgnore:dataLength);
end

ca3_inhibitoryNeuronExcitatoryCurrentData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
for k = 1: ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/ExcitatoryCurrent/exciteCurrent%d.txt', k);
    ca3_inhibitoryNeuronExcitatoryCurrentData{k} = csvread(myfilename);
    ca3_inhibitoryNeuronExcitatoryCurrentData{k}=ca3_inhibitoryNeuronExcitatoryCurrentData{k}(dataIgnore:dataLength);
end

ca3_inhibitoryNeuronBackgroundCurrentData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
for k = 1: ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/BackgroundCurrent/bg%d.txt', k);
    ca3_inhibitoryNeuronBackgroundCurrentData{k} = csvread(myfilename);
    ca3_inhibitoryNeuronBackgroundCurrentData{k}=ca3_inhibitoryNeuronBackgroundCurrentData{k}(dataIgnore:dataLength);
end

ca3_inhibitoryNeuronNMDACurrentData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
for k = 1: ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/NMDA/nmda%d.txt', k);
    ca3_inhibitoryNeuronNMDACurrentData{k} = csvread(myfilename);
    ca3_inhibitoryNeuronNMDACurrentData{k}=ca3_inhibitoryNeuronNMDACurrentData{k}(dataIgnore:dataLength);
end

if NMDA_Diagnostics == true

    ca3_inhibitoryNeuronNMDAConductanceData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/NMDA_Diagnostics/NMDA_Cond/nmda%d.txt', k);
        ca3_inhibitoryNeuronNMDAConductanceData{k} = csvread(myfilename);
        ca3_inhibitoryNeuronNMDAConductanceData{k}=ca3_inhibitoryNeuronNMDAConductanceData{k}(dataIgnore:dataLength);
    end
    
    ca3_inhibitoryNeuronNMDAmgBlockData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/NMDA_Diagnostics/mgBlock/nmda%d.txt', k);
        ca3_inhibitoryNeuronNMDAmgBlockData{k} = csvread(myfilename);
        ca3_inhibitoryNeuronNMDAmgBlockData{k}=ca3_inhibitoryNeuronNMDAmgBlockData{k}(dataIgnore:dataLength);
    end
    
    ca3_inhibitoryNeuronNMDAdactdtData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/NMDA_Diagnostics/dactdt/nmda%d.txt', k);
        ca3_inhibitoryNeuronNMDAdactdtData{k} = csvread(myfilename);
        ca3_inhibitoryNeuronNMDAdactdtData{k}=ca3_inhibitoryNeuronNMDAdactdtData{k}(dataIgnore:dataLength);
    end
    
    ca3_inhibitoryNeuronNMDAddeactdtData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/NMDA_Diagnostics/ddeactdt/nmda%d.txt', k);
        ca3_inhibitoryNeuronNMDAddeactdtData{k} = csvread(myfilename);
        ca3_inhibitoryNeuronNMDAddeactdtData{k}=ca3_inhibitoryNeuronNMDAddeactdtData{k}(dataIgnore:dataLength);
    end
    
    ca3_inhibitoryNeuronNMDAselfDeactData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/NMDA_Diagnostics/self.deact/nmda%d.txt', k);
        ca3_inhibitoryNeuronNMDAselfDeactData{k} = csvread(myfilename);
        ca3_inhibitoryNeuronNMDAselfDeactData{k}=ca3_inhibitoryNeuronNMDAselfDeactData{k}(dataIgnore:dataLength);
    end
    
    ca3_inhibitoryNeuronNMDAselfActData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/NMDA_Diagnostics/self.act/nmda%d.txt', k);
        ca3_inhibitoryNeuronNMDAselfActData{k} = csvread(myfilename);
        ca3_inhibitoryNeuronNMDAselfActData{k}=ca3_inhibitoryNeuronNMDAselfActData{k}(dataIgnore:dataLength);
    end
    
end

if Gate_Diagnostics == true
    
    ca3_inhibitoryNeuronHData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/gate_diagnostics/h/gate%d.txt', k);
        ca3_inhibitoryNeuronHData{k} = csvread(myfilename);
        ca3_inhibitoryNeuronHData{k}=ca3_inhibitoryNeuronHData{k}(dataIgnore:dataLength);
    end
    
    ca3_inhibitoryNeuronMData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/gate_diagnostics/m/gate%d.txt', k);
        ca3_inhibitoryNeuronMData{k} = csvread(myfilename);
        ca3_inhibitoryNeuronMData{k}=ca3_inhibitoryNeuronMData{k}(dataIgnore:dataLength);
    end
    
    ca3_inhibitoryNeuronNData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/gate_diagnostics/n/gate%d.txt', k);
        ca3_inhibitoryNeuronNData{k} = csvread(myfilename);
        ca3_inhibitoryNeuronNData{k}=ca3_inhibitoryNeuronNData{k}(dataIgnore:dataLength);
    end
    
    ca3_inhibitoryNeuronDvData = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/gate_diagnostics/dv/gate%d.txt', k);
        ca3_inhibitoryNeuronDvData{k} = csvread(myfilename);
        ca3_inhibitoryNeuronDvData{k}=ca3_inhibitoryNeuronDvData{k}(dataIgnore:dataLength);
    end
    
    ca3_inhibitoryNeuronDv1Data = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/gate_diagnostics/dv1/gate%d.txt', k);
        ca3_inhibitoryNeuronDv1Data{k} = csvread(myfilename);
        ca3_inhibitoryNeuronDv1Data{k}=ca3_inhibitoryNeuronDv1Data{k}(dataIgnore:dataLength);
    end
    
    ca3_inhibitoryNeuronDv2Data = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/gate_diagnostics/dv2/gate%d.txt', k);
        ca3_inhibitoryNeuronDv2Data{k} = csvread(myfilename);
        ca3_inhibitoryNeuronDv2Data{k}=ca3_inhibitoryNeuronDv2Data{k}(dataIgnore:dataLength);
    end

    ca3_inhibitoryNeuronDv3Data = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/gate_diagnostics/dv3/gate%d.txt', k);
        ca3_inhibitoryNeuronDv3Data{k} = csvread(myfilename);
        ca3_inhibitoryNeuronDv3Data{k}=ca3_inhibitoryNeuronDv3Data{k}(dataIgnore:dataLength);
    end
    
    ca3_inhibitoryNeuronDv4Data = cell(0,ca3_numberOfInhibitoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/InhibitoryNeuron/gate_diagnostics/dv4/gate%d.txt', k);
        ca3_inhibitoryNeuronDv4Data{k} = csvread(myfilename);
        ca3_inhibitoryNeuronDv4Data{k}=ca3_inhibitoryNeuronDv4Data{k}(dataIgnore:dataLength);
    end
end
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gather Excitatory Neuron Data
ca3_numberOfExcitatoryNeuronFiles = length(dir('NeuronNetworkClass/CA3/ExcitatoryNeuron/Voltage'))-2;

ca3_excitatoryNeuronVoltageData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/Voltage/voltage%d.txt', k);
    ca3_excitatoryNeuronVoltageData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
    ca3_excitatoryNeuronVoltageData{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronVoltageData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
end

ca3_excitatoryNeuronSpikeTimeData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/SpikeTime/spikeTime%d.txt', k);
    ca3_excitatoryNeuronSpikeTimeData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
    ca3_excitatoryNeuronSpikeTimeData{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronSpikeTimeData{k - ca3_numberOfInhibitoryNeuronFiles}(ca3_excitatoryNeuronSpikeTimeData{k - ca3_numberOfInhibitoryNeuronFiles}>timeData(1));
end

ca3_excitatoryNeuronInhibitoryCurrentData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/InhibitoryCurrent/inhibCurrent%d.txt', k);
    ca3_excitatoryNeuronInhibitoryCurrentData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
    ca3_excitatoryNeuronInhibitoryCurrentData{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronInhibitoryCurrentData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
end

ca3_excitatoryNeuronExcitatoryCurrentData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/ExcitatoryCurrent/exciteCurrent%d.txt', k);
    ca3_excitatoryNeuronExcitatoryCurrentData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
    ca3_excitatoryNeuronExcitatoryCurrentData{k - ca3_numberOfInhibitoryNeuronFiles}=ca3_excitatoryNeuronExcitatoryCurrentData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
end

ca3_excitatoryNeuronBackgroundCurrentData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/BackgroundCurrent/bg%d.txt', k);
    ca3_excitatoryNeuronBackgroundCurrentData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
    ca3_excitatoryNeuronBackgroundCurrentData{k - ca3_numberOfInhibitoryNeuronFiles}=ca3_excitatoryNeuronBackgroundCurrentData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
end

ca3_excitatoryNeuronNMDACurrentData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
    myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/NMDA/nmda%d.txt', k);
    ca3_excitatoryNeuronNMDACurrentData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
    ca3_excitatoryNeuronNMDACurrentData{k - ca3_numberOfInhibitoryNeuronFiles}=ca3_excitatoryNeuronNMDACurrentData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
end


if NMDA_Diagnostics == true

    ca3_excitatoryNeuronNMDAConductanceData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/NMDA_Diagnostics/NMDA_Cond/nmda%d.txt', k);
        ca3_excitatoryNeuronNMDAConductanceData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronNMDAConductanceData{k - ca3_numberOfInhibitoryNeuronFiles}=ca3_excitatoryNeuronNMDAConductanceData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronNMDAmgBlockData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/NMDA_Diagnostics/mgBlock/nmda%d.txt', k);
        ca3_excitatoryNeuronNMDAmgBlockData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronNMDAmgBlockData{k - ca3_numberOfInhibitoryNeuronFiles}=ca3_excitatoryNeuronNMDAmgBlockData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronNMDAdactdtData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/NMDA_Diagnostics/dactdt/nmda%d.txt', k);
        ca3_excitatoryNeuronNMDAdactdtData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronNMDAdactdtData{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronNMDAdactdtData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronNMDAselfDeactData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/NMDA_Diagnostics/self.deact/nmda%d.txt', k);
        ca3_excitatoryNeuronNMDAselfDeactData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronNMDAselfDeactData{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronNMDAselfDeactData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronNMDAddeactdtData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/NMDA_Diagnostics/ddeactdt/nmda%d.txt', k);
        ca3_excitatoryNeuronNMDAddeactdtData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronNMDAddeactdtData{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronNMDAddeactdtData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronNMDAselfActData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/NMDA_Diagnostics/self.act/nmda%d.txt', k);
        ca3_excitatoryNeuronNMDAselfActData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronNMDAselfActData{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronNMDAselfActData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    

end

if Gate_Diagnostics == true
    ca3_excitatoryNeuronHData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/gate_diagnostics/h/gate%d.txt', k);
        ca3_excitatoryNeuronHData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronHData{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronHData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronMData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/gate_diagnostics/m/gate%d.txt', k);
        ca3_excitatoryNeuronMData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronMData{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronMData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronNData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/gate_diagnostics/n/gate%d.txt', k);
        ca3_excitatoryNeuronNData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronNData{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronNData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronDvData = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/gate_diagnostics/dv/gate%d.txt', k);
        ca3_excitatoryNeuronDvData{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronDvData{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronDvData{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronDv1Data = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/gate_diagnostics/dv1/gate%d.txt', k);
        ca3_excitatoryNeuronDv1Data{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronDv1Data{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronDv1Data{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronDv2Data = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/gate_diagnostics/dv2/gate%d.txt', k);
        ca3_excitatoryNeuronDv2Data{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronDv2Data{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronDv2Data{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronDv3Data = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/gate_diagnostics/dv3/gate%d.txt', k);
        ca3_excitatoryNeuronDv3Data{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronDv3Data{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronDv3Data{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
    
    ca3_excitatoryNeuronDv4Data = cell(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca3_numberOfInhibitoryNeuronFiles: ca3_numberOfExcitatoryNeuronFiles + ca3_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA3/ExcitatoryNeuron/gate_diagnostics/dv4/gate%d.txt', k);
        ca3_excitatoryNeuronDv4Data{k - ca3_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca3_excitatoryNeuronDv4Data{k - ca3_numberOfInhibitoryNeuronFiles} = ca3_excitatoryNeuronDv4Data{k - ca3_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%                                   CA1



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if CA1_ON == true
    % Gather Inhibitory Neuron Data
    ca1_numberOfInhibitoryNeuronFiles = length(dir('NeuronNetworkClass/CA1/InhibitoryNeuron/Voltage'))-2;

    ca1_inhibitoryNeuronVoltageData = cell(0,ca1_numberOfInhibitoryNeuronFiles);
    for k = 1: ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/InhibitoryNeuron/Voltage/voltage%d.txt', k);
        ca1_inhibitoryNeuronVoltageData{k} = csvread(myfilename);
        ca1_inhibitoryNeuronVoltageData{k} = ca1_inhibitoryNeuronVoltageData{k}(dataIgnore:dataLength);
    end

    ca1_inhibitoryNeuronSpikeTimeData = cell(0,ca1_numberOfInhibitoryNeuronFiles);
    for k = 1: ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/InhibitoryNeuron/SpikeTime/spikeTime%d.txt', k);
        ca1_inhibitoryNeuronSpikeTimeData{k} = csvread(myfilename);
        ca1_inhibitoryNeuronSpikeTimeData{k} = ca1_inhibitoryNeuronSpikeTimeData{k}(ca1_inhibitoryNeuronSpikeTimeData{k}>timeData(1));
    end

    ca1_inhibitoryNeuronInhibitoryCurrentData = cell(0,ca1_numberOfInhibitoryNeuronFiles);
    for k = 1: ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/InhibitoryNeuron/InhibitoryCurrent/inhibCurrent%d.txt', k);
        ca1_inhibitoryNeuronInhibitoryCurrentData{k} = csvread(myfilename);
        ca1_inhibitoryNeuronInhibitoryCurrentData{k} = ca1_inhibitoryNeuronInhibitoryCurrentData{k}(dataIgnore:dataLength);
    end

    ca1_inhibitoryNeuronExcitatoryCurrentData = cell(0,ca1_numberOfInhibitoryNeuronFiles);
    for k = 1: ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/InhibitoryNeuron/ExcitatoryCurrent/exciteCurrent%d.txt', k);
        ca1_inhibitoryNeuronExcitatoryCurrentData{k} = csvread(myfilename);
        ca1_inhibitoryNeuronExcitatoryCurrentData{k} = ca1_inhibitoryNeuronExcitatoryCurrentData{k}(dataIgnore:dataLength);
    end

    ca1_inhibitoryNeuronBackgroundCurrentData = cell(0,ca1_numberOfInhibitoryNeuronFiles);
    for k = 1: ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/InhibitoryNeuron/BackgroundCurrent/bg%d.txt', k);
        ca1_inhibitoryNeuronBackgroundCurrentData{k} = csvread(myfilename);
        ca1_inhibitoryNeuronBackgroundCurrentData{k} = ca1_inhibitoryNeuronBackgroundCurrentData{k}(dataIgnore:dataLength);
    end

    ca1_inhibitoryNeuronReverseConnectionData = cell(0,ca1_numberOfInhibitoryNeuronFiles);
    for k = 1: ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/InhibitoryNeuron/ReverseConnections/neuron%d.txt', k);
        ca1_inhibitoryNeuronReverseConnectionData{k} = csvread(myfilename);

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Gather Excitatory Neuron Data
    ca1_numberOfExcitatoryNeuronFiles = length(dir('NeuronNetworkClass/CA1/ExcitatoryNeuron/Voltage'))-2;

    ca1_excitatoryNeuronVoltageData = cell(0,ca1_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca1_numberOfInhibitoryNeuronFiles: ca1_numberOfExcitatoryNeuronFiles + ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/ExcitatoryNeuron/Voltage/voltage%d.txt', k);
        ca1_excitatoryNeuronVoltageData{k - ca1_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca1_excitatoryNeuronVoltageData{k - ca1_numberOfInhibitoryNeuronFiles} = ca1_excitatoryNeuronVoltageData{k - ca1_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end

    ca1_excitatoryNeuronSpikeTimeData = cell(0,ca1_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca1_numberOfInhibitoryNeuronFiles: ca1_numberOfExcitatoryNeuronFiles + ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/ExcitatoryNeuron/SpikeTime/spikeTime%d.txt', k);
        ca1_excitatoryNeuronSpikeTimeData{k - ca1_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca1_excitatoryNeuronSpikeTimeData{k - ca1_numberOfInhibitoryNeuronFiles} = ca1_excitatoryNeuronSpikeTimeData{k - ca1_numberOfInhibitoryNeuronFiles}(ca1_excitatoryNeuronSpikeTimeData{k - ca1_numberOfInhibitoryNeuronFiles}>timeData(1));
    end

    ca1_excitatoryNeuronInhibitoryCurrentData = cell(0,ca1_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca1_numberOfInhibitoryNeuronFiles: ca1_numberOfExcitatoryNeuronFiles + ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/ExcitatoryNeuron/InhibitoryCurrent/inhibCurrent%d.txt', k);
        ca1_excitatoryNeuronInhibitoryCurrentData{k - ca1_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca1_excitatoryNeuronInhibitoryCurrentData{k - ca1_numberOfInhibitoryNeuronFiles} = ca1_excitatoryNeuronInhibitoryCurrentData{k - ca1_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end

    ca1_excitatoryNeuronExcitatoryCurrentData = cell(0,ca1_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca1_numberOfInhibitoryNeuronFiles: ca1_numberOfExcitatoryNeuronFiles + ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/ExcitatoryNeuron/ExcitatoryCurrent/exciteCurrent%d.txt', k);
        ca1_excitatoryNeuronExcitatoryCurrentData{k - ca1_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca1_excitatoryNeuronExcitatoryCurrentData{k - ca1_numberOfInhibitoryNeuronFiles} = ca1_excitatoryNeuronExcitatoryCurrentData{k - ca1_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end

    ca1_excitatoryNeuronBackgroundCurrentData = cell(0,ca1_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca1_numberOfInhibitoryNeuronFiles: ca1_numberOfExcitatoryNeuronFiles + ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/ExcitatoryNeuron/BackgroundCurrent/bg%d.txt', k);
        ca1_excitatoryNeuronBackgroundCurrentData{k - ca1_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
        ca1_excitatoryNeuronBackgroundCurrentData{k - ca1_numberOfInhibitoryNeuronFiles} = ca1_excitatoryNeuronBackgroundCurrentData{k - ca1_numberOfInhibitoryNeuronFiles}(dataIgnore:dataLength);
    end

    ca1_excitatoryNeuronReverseConnectionData = cell(0,ca1_numberOfExcitatoryNeuronFiles);
    for k = 1 + ca1_numberOfInhibitoryNeuronFiles: ca1_numberOfExcitatoryNeuronFiles + ca1_numberOfInhibitoryNeuronFiles
        myfilename = sprintf('NeuronNetworkClass/CA1/ExcitatoryNeuron/ReverseConnections/neuron%d.txt', k);
        ca1_excitatoryNeuronReverseConnectionData{k - ca1_numberOfInhibitoryNeuronFiles} = csvread(myfilename);
    end
end

%Cleaning
clear dataIgnore
clear dataLength
clear myfilename
clear k