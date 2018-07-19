%%%%%% Acidental Run Protection %%%%%%%%%%%%%%%
if multiRun == false
    x = input('WARNING: Are you sure you want to process data for both neuron networks? (y/n) ');

    if x == 'y'
        disp('Okay, running')
        clear x
    else
        clear x
        disp('Process Aborted')
        return
    end
end


%%% Shut off sdr processing
sdrOff = true;


%%%%%%%%%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Processing for Total Field Potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%  Calculating CA3 Field Potentials  %%%%%

% Lists containing averages
ca3_PNFieldPotential = zeros(1,length(timeData));
ca3_BasketFieldPotential = zeros(1,length(timeData));
ca3_TotalFieldPotential = zeros(1,length(timeData));

% Calculating the total average
for z = 1: length(timeData)
    pnNeuronVoltageOfASingleDt = zeros(1,length(ca3_excitatoryNeuronVoltageData));
    for k = 1: length(ca3_excitatoryNeuronVoltageData)
        pnNeuronVoltageOfASingleDt(k) = ca3_excitatoryNeuronVoltageData{k}(z);
    end

    bCellVoltageOfASingleDt = zeros(1,length(ca3_inhibitoryNeuronVoltageData));
    for l = 1: length(ca3_inhibitoryNeuronVoltageData)
        bCellVoltageOfASingleDt(l) = ca3_inhibitoryNeuronVoltageData{l}(z);
    end
    
    ca3_PNFieldPotential(z) = mean(pnNeuronVoltageOfASingleDt);
    ca3_BasketFieldPotential(z) = mean(bCellVoltageOfASingleDt);
    
    allNeuronVoltageOfASingleDt = [pnNeuronVoltageOfASingleDt bCellVoltageOfASingleDt];
    ca3_TotalFieldPotential(z) = mean(allNeuronVoltageOfASingleDt);
end


%%%% Calculating CA1 Field Potentials %%%%%%%%%%%%

if CA1_ON == true

    % Lists containing averages
    ca1_PNFieldPotential = zeros(1,length(timeData));
    ca1_BasketFieldPotential = zeros(1,length(timeData));
    ca1_TotalFieldPotential = zeros(1,length(timeData));

    % Calculating the total average
    for z = 1: length(timeData)
        pnNeuronVoltageOfASingleDt = zeros(1,length(ca1_excitatoryNeuronVoltageData));
        for k = 1: length(ca1_excitatoryNeuronVoltageData)
            pnNeuronVoltageOfASingleDt(k) = ca1_excitatoryNeuronVoltageData{k}(z);
        end

        bCellVoltageOfASingleDt = zeros(1,length(ca1_inhibitoryNeuronVoltageData));
        for l = 1: length(ca1_inhibitoryNeuronVoltageData)
            bCellVoltageOfASingleDt(l) = ca1_inhibitoryNeuronVoltageData{l}(z);
        end

        ca1_PNFieldPotential(z) = mean(pnNeuronVoltageOfASingleDt);
        ca1_BasketFieldPotential(z) = mean(bCellVoltageOfASingleDt);

        allNeuronVoltageOfASingleDt = [pnNeuronVoltageOfASingleDt bCellVoltageOfASingleDt];
        ca1_TotalFieldPotential(z) = mean(allNeuronVoltageOfASingleDt);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Processing for Signal Iron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CA3 Signal Iron

% This smoothens the total field potential and graphs it
ca3_smootheningResolution = 50;  %40 for ca1 build
ca3_yValue = -50;
ca3_timeIntervalFilter = 10;
ca3_heightFilter = -62.5;  %-0.1 for ca1 build;
ca3_coeff24hMA = ones(1, ca3_smootheningResolution)/ca3_smootheningResolution;
ca3_fDelay = (length(ca3_coeff24hMA)-1)/2;
ca3_avgSmoothFieldPotential = filter(ca3_coeff24hMA, 1, ca3_TotalFieldPotential);

% This finds the peak times in the smoothened version
ca3_listOfPeakTimes = zeros(0,0);
ca3_listOfPeakTimesYVector = zeros(0,0);
for k = 2: length(ca3_avgSmoothFieldPotential)
    if (k ~= length(ca3_avgSmoothFieldPotential)) % Prevent to high an index
        if isempty(ca3_listOfPeakTimes) % block to check if empty to run checks without looking at past peak time
            if (ca3_avgSmoothFieldPotential(k-1) <= ca3_avgSmoothFieldPotential(k)) && (ca3_avgSmoothFieldPotential(k+1) < ca3_avgSmoothFieldPotential(k)) && (ca3_avgSmoothFieldPotential(k) > ca3_heightFilter)
                peakTime = timeData(k)-ca3_fDelay/ca3_smootheningResolution;
                ca3_listOfPeakTimes = [ca3_listOfPeakTimes, zeros(0,1)];
                ca3_listOfPeakTimes(length(ca3_listOfPeakTimes)+1) = peakTime;

                ca3_listOfPeakTimesYVector = [ca3_listOfPeakTimesYVector, zeros(0,1)];
                ca3_listOfPeakTimesYVector(length(ca3_listOfPeakTimesYVector)+1) = ca3_yValue;  
            end
        else % checks that include looking at past peak time
            if (ca3_avgSmoothFieldPotential(k-1) <= ca3_avgSmoothFieldPotential(k)) && (ca3_avgSmoothFieldPotential(k+1) < ca3_avgSmoothFieldPotential(k)) && (ca3_avgSmoothFieldPotential(k) > ca3_heightFilter) && ((timeData(k)-ca3_fDelay/ca3_smootheningResolution) - ca3_listOfPeakTimes(length(ca3_listOfPeakTimes)) > ca3_timeIntervalFilter)
                peakTime = timeData(k)-ca3_fDelay/ca3_smootheningResolution;
                ca3_listOfPeakTimes = [ca3_listOfPeakTimes, zeros(0,1)];
                ca3_listOfPeakTimes(length(ca3_listOfPeakTimes)+1) = peakTime;

                ca3_listOfPeakTimesYVector = [ca3_listOfPeakTimesYVector, zeros(0,1)];
                ca3_listOfPeakTimesYVector(length(ca3_listOfPeakTimesYVector)+1) = ca3_yValue;
            end
        end
    end
end

if length(ca3_listOfPeakTimes) < 1
    warning('You Have 0 or 1 CA3 cycle(s). Many Calculations may fail as a result.')
end

%%% CA1 Signal Iron
if CA1_ON == true

    % This smoothens the total field potential and graphs it
    ca1_smootheningResolution = 40;
    ca1_yValue = .4;
    ca1_timeIntervalFilter = 10;
    ca1_heightFilter = -0.3;
    ca1_coeff24hMA = ones(1, ca1_smootheningResolution)/ca1_smootheningResolution;
    ca1_fDelay = (length(ca1_coeff24hMA)-1)/2;
    ca1_avgSmoothFieldPotential = filter(ca1_coeff24hMA, 1, ca1_TotalFieldPotential);

    % This finds the peak times in the smoothened version
    ca1_listOfPeakTimes = zeros(0,0);
    ca1_listOfPeakTimesYVector = zeros(0,0);
    for k = 2: length(ca1_avgSmoothFieldPotential)
        if (k ~= length(ca1_avgSmoothFieldPotential)) % Prevent to high an index
            if isempty(ca1_listOfPeakTimes) % block to check if empty to run checks without looking at past peak time
                if (ca1_avgSmoothFieldPotential(k-1) <= ca1_avgSmoothFieldPotential(k)) && (ca1_avgSmoothFieldPotential(k+1) < ca1_avgSmoothFieldPotential(k)) && (ca1_avgSmoothFieldPotential(k) > ca1_heightFilter)
                    peakTime = timeData(k)-ca1_fDelay/ca1_smootheningResolution;
                    ca1_listOfPeakTimes = [ca1_listOfPeakTimes, zeros(0,1)];
                    ca1_listOfPeakTimes(length(ca1_listOfPeakTimes)+1) = peakTime;

                    ca1_listOfPeakTimesYVector = [ca1_listOfPeakTimesYVector, zeros(0,1)];
                    ca1_listOfPeakTimesYVector(length(ca1_listOfPeakTimesYVector)+1) = ca1_yValue;  
                end
            else % checks that include looking at past peak time
                if (ca1_avgSmoothFieldPotential(k-1) <= ca1_avgSmoothFieldPotential(k)) && (ca1_avgSmoothFieldPotential(k+1) < ca1_avgSmoothFieldPotential(k)) && (ca1_avgSmoothFieldPotential(k) > ca1_heightFilter) && ((timeData(k)-ca1_fDelay/ca1_smootheningResolution) - ca1_listOfPeakTimes(length(ca1_listOfPeakTimes)) > ca1_timeIntervalFilter)
                    peakTime = timeData(k)-ca1_fDelay/ca1_smootheningResolution;
                    ca1_listOfPeakTimes = [ca1_listOfPeakTimes, zeros(0,1)];
                    ca1_listOfPeakTimes(length(ca1_listOfPeakTimes)+1) = peakTime;

                    ca1_listOfPeakTimesYVector = [ca1_listOfPeakTimesYVector, zeros(0,1)];
                    ca1_listOfPeakTimesYVector(length(ca1_listOfPeakTimesYVector)+1) = ca1_yValue;
                end
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Fourrier Tansform %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CA3 Fourrier

Fs = 1000/(timeData(2)-timeData(1));   % Sampling frequency                    
%T  = 1/Fs;                             % Sampling period       
L  = length(ca3_TotalFieldPotential);  % Length of signal
%t  = (0:L-1)*T;                        % Time vector

X = ca3_TotalFieldPotential - mean(ca3_TotalFieldPotential); % x = x-mean(x)
ca3_fourrierTransformOfFieldPotential = fft(X);

ca3_P2 = abs(ca3_fourrierTransformOfFieldPotential/L);
ca3_P1 = ca3_P2(1:L/2+1);
ca3_P1(2:end-1) = 2*ca3_P1(2:end-1);
ca3_frequencyVectorForTransform = Fs*(0:(L/2))/L;


%%% CA1 Fourrier

if CA1_ON == true

    Fs = 1000/(timeData(2)-timeData(1));   % Sampling frequency                    
    T  = 1/Fs;                             % Sampling period       
    L  = length(ca1_TotalFieldPotential);  % Length of signal
    %t  = (0:L-1)*T;                        % Time vector

    X = ca1_TotalFieldPotential - mean(ca1_TotalFieldPotential); % x = x-mean(x)
    ca1_fourrierTransformOfFieldPotential = fft(X);

    ca1_P2 = abs(ca1_fourrierTransformOfFieldPotential/L);
    ca1_P1 = ca1_P2(1:L/2+1);
    ca1_P1(2:end-1) = 2*ca1_P1(2:end-1);
    ca1_frequencyVectorForTransform = Fs*(0:(L/2))/L;

end


%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Cycle Processor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Cycle Processor finds which neurons spiked in a cycle and when they
% spiked. It does this for every cycle. The resulting variable is cycle
% spike data and it has n cells for n cycles. In each cell there is a 2xn
% matrix. The first row is the id of the neuron and the second row is the
% exact time of the spike 
%
% Note: pnCycleSpike is a variant of cycle spike data, but here there are
% only PNs and no INs. This will be used in jitter calculations


%%% CA3 Cycle Processor
spikeTimeData = horzcat(ca3_inhibitoryNeuronSpikeTimeData,ca3_excitatoryNeuronSpikeTimeData);

ca3_longestList = 0;
for k = 1:length(spikeTimeData)
    if length(spikeTimeData{k}) > ca3_longestList
        ca3_longestList = length(spikeTimeData{k});
    end
end

ca3_pnLongestList = 0;
for k = 1:length(ca3_excitatoryNeuronSpikeTimeData)
    if length(ca3_excitatoryNeuronSpikeTimeData{k}) > ca3_pnLongestList
        ca3_pnLongestList = length(ca3_excitatoryNeuronSpikeTimeData{k});
    end
end

ca3_cycleSpikeData = cell(0,length(ca3_listOfPeakTimes)-1);
for k = 1: length(ca3_listOfPeakTimes)-1 % scan for each cycle
    singleCycleData = zeros(2,0);
    zeroLocation = 1; % counter so we can insert the data into right slot for single cycleData
    x = 1;
    while (x <= ca3_longestList) % scan through each neuron's data

        z = 1;
        while (z <= length(spikeTimeData))  % scan for each neuron

            % pre allocated for 10 spikes:
            if x <= length(spikeTimeData{z}) % if data point iterator is not longer than the list
            
                
                if (spikeTimeData{z}(x) > ca3_listOfPeakTimes(k)) && (spikeTimeData{z}(x) <= ca3_listOfPeakTimes(k+1)) % if the data point fits the criterion

                    singleCycleData(1,zeroLocation) = spikeTimeData{z}(x);
                    singleCycleData(2,zeroLocation) = z;
                    zeroLocation = zeroLocation + 1;
     
                end
       
            end
        
        z = z + 1;    
        end
        
        
    x = x + 1;    
    end
    
ca3_cycleSpikeData{k} = singleCycleData;
end

ca3_PNcycleSpikeData = cell(0,length(ca3_listOfPeakTimes)-1);
for k = 1: length(ca3_listOfPeakTimes)-1 % scan for each cycle
    singleCycleData = zeros(2,0);
    zeroLocation = 1; % counter so we can insert the data into right slot for single cycleData
    x = 1;
    while (x <= ca3_pnLongestList) % scan through each neuron's data

        z = 1;
        while (z <= length(ca3_excitatoryNeuronSpikeTimeData))  % scan for each neuron

            % pre allocated for 10 spikes:
            if x <= length(ca3_excitatoryNeuronSpikeTimeData{z}) % if data point iterator is not longer than the list
            
                
                if (ca3_excitatoryNeuronSpikeTimeData{z}(x) > ca3_listOfPeakTimes(k)) && (ca3_excitatoryNeuronSpikeTimeData{z}(x) <= ca3_listOfPeakTimes(k+1)) % if the data point fits the criterion

                    singleCycleData(1,zeroLocation) = ca3_excitatoryNeuronSpikeTimeData{z}(x);
                    singleCycleData(2,zeroLocation) = z;
                    zeroLocation = zeroLocation + 1;
     
                end
       
            end
        
        z = z + 1;    
        end
        
        
    x = x + 1;    
    end
    
ca3_PNcycleSpikeData{k} = singleCycleData;
end


%%% CA1 Cycle Processor

if CA1_ON == true

    spikeTimeData = horzcat(ca1_inhibitoryNeuronSpikeTimeData,ca1_excitatoryNeuronSpikeTimeData);

    ca1_longestList = 0;
    for k = 1:length(spikeTimeData)
        if length(spikeTimeData{k}) > ca1_longestList
            ca1_longestList = length(spikeTimeData{k});
        end
    end

    ca1_pnLongestList = 0;
    for k = 1:length(ca1_excitatoryNeuronSpikeTimeData)
        if length(ca1_excitatoryNeuronSpikeTimeData{k}) > ca1_pnLongestList
            ca1_pnLongestList = length(ca1_excitatoryNeuronSpikeTimeData{k});
        end
    end

    ca1_cycleSpikeData = cell(0,length(ca1_listOfPeakTimes)-1);
    for k = 1: length(ca1_listOfPeakTimes)-1 % scan for each cycle
        singleCycleData = zeros(2,0);
        zeroLocation = 1; % counter so we can insert the data into right slot for single cycleData
        x = 1;
        while (x <= ca1_longestList) % scan through each neuron's data

            z = 1;
            while (z <= length(spikeTimeData))  % scan for each neuron

                % pre allocated for 10 spikes:
                if x <= length(spikeTimeData{z}) % if data point iterator is not longer than the list


                    if (spikeTimeData{z}(x) > ca1_listOfPeakTimes(k)) && (spikeTimeData{z}(x) <= ca1_listOfPeakTimes(k+1)) % if the data point fits the criterion

                        singleCycleData(1,zeroLocation) = spikeTimeData{z}(x);
                        singleCycleData(2,zeroLocation) = z;
                        zeroLocation = zeroLocation + 1;

                    end

                end

            z = z + 1;    
            end


        x = x + 1;    
        end

    ca1_cycleSpikeData{k} = singleCycleData;

    end

    ca1_PNcycleSpikeData = cell(0,length(ca1_listOfPeakTimes)-1);
    for k = 1: length(ca1_listOfPeakTimes)-1 % scan for each cycle
        singleCycleData = zeros(2,0);
        zeroLocation = 1; % counter so we can insert the data into right slot for single cycleData
        x = 1;
        while (x <= ca1_pnLongestList) % scan through each neuron's data

            z = 1;
            while (z <= length(ca1_excitatoryNeuronSpikeTimeData))  % scan for each neuron

                % pre allocated for 10 spikes:
                if x <= length(ca1_excitatoryNeuronSpikeTimeData{z}) % if data point iterator is not longer than the list


                    if (ca1_excitatoryNeuronSpikeTimeData{z}(x) > ca1_listOfPeakTimes(k)) && (ca1_excitatoryNeuronSpikeTimeData{z}(x) <= ca1_listOfPeakTimes(k+1)) % if the data point fits the criterion

                        singleCycleData(1,zeroLocation) = ca1_excitatoryNeuronSpikeTimeData{z}(x);
                        singleCycleData(2,zeroLocation) = z;
                        zeroLocation = zeroLocation + 1;

                    end

                end

            z = z + 1;    
            end


        x = x + 1;    
        end

    ca1_PNcycleSpikeData{k} = singleCycleData;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Spikes Per Cycle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Spikes Per cycle shows each neurons data for how many times it spiked for
% every cycle. Its pretty similar to cycle spike data. The difference is
% basically how they're organized. If looking for data specifically for a
% cycle, cycle spike data is the best bet for finding that info. But lets
% say you want spike info on a particular neuron, then spikes per cycle is
% the best bet. So the data is organized like this: there are cells for
% every neuron, and in each cell there is an array of size n for n cycles
% and in the array, a tally for the amount of times the neuron spiked for
% that cycle

%%% CA3 Spikes Per Cycle Data

% Excitatory Neurons
ca3_listOfNumberOfExcitatoryNeuronSpikesPerCycle = cell(0,ca3_numberOfExcitatoryNeuronFiles);
for k = 1: ca3_numberOfExcitatoryNeuronFiles
    ca3_listOfNumberOfExcitatoryNeuronSpikesPerCycle{k} = zeros(0,length(ca3_listOfPeakTimes)-1);
end
for k = 2: length(ca3_listOfPeakTimes)
    for y = 1:ca3_numberOfExcitatoryNeuronFiles
        z = 1;
        numberOfSpikes = 0;
        while (z <= length(ca3_excitatoryNeuronSpikeTimeData{y})) && (ca3_excitatoryNeuronSpikeTimeData{y}(z) < ca3_listOfPeakTimes(k))
            
            if ca3_listOfPeakTimes(k-1) <= ca3_excitatoryNeuronSpikeTimeData{y}(z)
                numberOfSpikes = numberOfSpikes + 1;
            end
            z = z+1;
        end
        ca3_listOfNumberOfExcitatoryNeuronSpikesPerCycle{y}(k-1) = numberOfSpikes;
    end   
end

% Inhibitory Neurons
ca3_listOfNumberOfInhibitoryNeuronSpikesPerCycle = cell(0,ca3_numberOfInhibitoryNeuronFiles);
for k = 1: ca3_numberOfInhibitoryNeuronFiles
    ca3_listOfNumberOfInhibitoryNeuronSpikesPerCycle{k} = zeros(0,length(ca3_listOfPeakTimes)-1);
end
for k = 2: length(ca3_listOfPeakTimes)
    for y = 1:ca3_numberOfInhibitoryNeuronFiles
        z = 1;
        numberOfSpikes = 0;
        while (z <= length(ca3_inhibitoryNeuronSpikeTimeData{y})) && (ca3_inhibitoryNeuronSpikeTimeData{y}(z) < ca3_listOfPeakTimes(k))
            
            if ca3_listOfPeakTimes(k-1) <= ca3_inhibitoryNeuronSpikeTimeData{y}(z)
                numberOfSpikes = numberOfSpikes + 1;
            end
            z = z+1;
        end
        ca3_listOfNumberOfInhibitoryNeuronSpikesPerCycle{y}(k-1) = numberOfSpikes;
    end   
end


%%% CA1 Spikes Per Cycle Data

if CA1_ON == true

    % Excitatory Neurons
    ca1_listOfNumberOfExcitatoryNeuronSpikesPerCycle = cell(0,ca1_numberOfExcitatoryNeuronFiles);
    for k = 1: ca1_numberOfExcitatoryNeuronFiles
        ca1_listOfNumberOfExcitatoryNeuronSpikesPerCycle{k} = zeros(0,length(ca1_listOfPeakTimes)-1);
    end
    for k = 2: length(ca1_listOfPeakTimes)
        for y = 1:ca1_numberOfExcitatoryNeuronFiles
            z = 1;
            numberOfSpikes = 0;
            while (z <= length(ca1_excitatoryNeuronSpikeTimeData{y})) && (ca1_excitatoryNeuronSpikeTimeData{y}(z) < ca1_listOfPeakTimes(k))

                if ca1_listOfPeakTimes(k-1) <= ca1_excitatoryNeuronSpikeTimeData{y}(z)
                    numberOfSpikes = numberOfSpikes + 1;
                end
                z = z+1;
            end
            ca1_listOfNumberOfExcitatoryNeuronSpikesPerCycle{y}(k-1) = numberOfSpikes;
        end   
    end     

    % Inhibitory Neurons
    ca1_listOfNumberOfInhibitoryNeuronSpikesPerCycle = cell(0,ca1_numberOfInhibitoryNeuronFiles);
    for k = 1: ca1_numberOfInhibitoryNeuronFiles
        ca1_listOfNumberOfInhibitoryNeuronSpikesPerCycle{k} = zeros(0,length(ca1_listOfPeakTimes)-1);
    end
    for k = 2: length(ca1_listOfPeakTimes)
        for y = 1:ca1_numberOfInhibitoryNeuronFiles
            z = 1;
            numberOfSpikes = 0;
            while (z <= length(ca1_inhibitoryNeuronSpikeTimeData{y})) && (ca1_inhibitoryNeuronSpikeTimeData{y}(z) < ca1_listOfPeakTimes(k))

                if ca1_listOfPeakTimes(k-1) <= ca1_inhibitoryNeuronSpikeTimeData{y}(z)
                    numberOfSpikes = numberOfSpikes + 1;
                end
                z = z+1;
            end
            ca1_listOfNumberOfInhibitoryNeuronSpikesPerCycle{y}(k-1) = numberOfSpikes;
        end   
    end

end

%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Hz Per Cycle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CA3 Hz Per Cycle

% Calculating the average firing rate per cycle of the excitatory neuron
ca3_listOfENeuronSpikesPerCycle = zeros(0,length(ca3_listOfPeakTimes)-1);
for z = 1: length(ca3_listOfPeakTimes)-1
    pnNeuronSpikesOfASingleCycle = zeros(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1: ca3_numberOfExcitatoryNeuronFiles
        pnNeuronSpikesOfASingleCycle(k) = ca3_listOfNumberOfExcitatoryNeuronSpikesPerCycle{k}(z);
    end   
    ca3_listOfENeuronSpikesPerCycle(z) = sum(pnNeuronSpikesOfASingleCycle);
end 

% Calculating the average firing rate per cycle of the inhibitory neuron
ca3_listOfINeuronSpikesPerCycle = zeros(0,length(ca3_listOfPeakTimes)-1);
for z = 1: length(ca3_listOfPeakTimes)-1
    InNeuronSpikesOfASingleCycle = zeros(0,ca3_numberOfExcitatoryNeuronFiles);
    for k = 1: ca3_numberOfInhibitoryNeuronFiles
        InNeuronSpikesOfASingleCycle(k) = ca3_listOfNumberOfInhibitoryNeuronSpikesPerCycle{k}(z);
    end   
    ca3_listOfINeuronSpikesPerCycle(z) = sum(InNeuronSpikesOfASingleCycle);
end


%%% CA1 Hz Per Cycle

if CA1_ON == true

    % Calculating the average firing rate per cycle of the excitatory neuron
    ca1_listOfENeuronSpikesPerCycle = zeros(0,length(ca1_listOfPeakTimes)-1);
    for z = 1: length(ca1_listOfPeakTimes)-1
        pnNeuronSpikesOfASingleCycle = zeros(0,ca1_numberOfExcitatoryNeuronFiles);
        for k = 1: ca1_numberOfExcitatoryNeuronFiles
            pnNeuronSpikesOfASingleCycle(k) = ca1_listOfNumberOfExcitatoryNeuronSpikesPerCycle{k}(z);
        end   
        ca1_listOfENeuronSpikesPerCycle(z) = sum(pnNeuronSpikesOfASingleCycle);
    end 

    % Calculating the average firing rate per cycle of the inhibitory neuron
    ca1_listOfINeuronSpikesPerCycle = zeros(0,length(ca1_listOfPeakTimes)-1);
    for z = 1: length(ca1_listOfPeakTimes)-1
        InNeuronSpikesOfASingleCycle = zeros(0,ca1_numberOfExcitatoryNeuronFiles);
        for k = 1: ca1_numberOfInhibitoryNeuronFiles
            InNeuronSpikesOfASingleCycle(k) = ca1_listOfNumberOfInhibitoryNeuronSpikesPerCycle{k}(z);
        end   
        ca1_listOfINeuronSpikesPerCycle(z) = sum(InNeuronSpikesOfASingleCycle);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% SDR Calculator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sdrOff == false
    %%% CA3 SDR 
    try
        ca3_SDRforEachDistance = cell(2,length(ca3_cycleSpikeData)-1);
        for distance = 1:length(ca3_cycleSpikeData) - 1
            ca3_SDRforEachDistance{distance} = zeros(2,length(ca3_cycleSpikeData) - distance);

            for q = 1:length(ca3_cycleSpikeData) - distance
                x = q;
                y = q+distance;

                n = unique(ca3_cycleSpikeData{x}(2,:));
                n = n(n>50);
                k = unique(ca3_cycleSpikeData{y}(2,:));
                k = k(k>50);
                r = length(union(n,k));
                s = length(intersect(n,k));
                n = length(n);
                k = length(k);

                sdr = (r-s)/(n+k);
                sdr_adjusted = (r-s)/(n+k) - abs(n-k)/(n+k);
                if isnan(sdr)
                    sdr = 1;
                    sdr_adjusted = 1;
                end
                ca3_SDRforEachDistance{distance}(1,q) = sdr;
                ca3_SDRforEachDistance{distance}(2,q) = sdr_adjusted;

            end


        end
        ca3_sdrBar = mean(ca3_SDRforEachDistance{1}(1,:));
        ca3_weightSDRBAr = mean(ca3_SDRforEachDistance{1}(2,:));
        ca3_sdrBarErr = std(ca3_SDRforEachDistance{1}(1,:));
        ca3_weightSdrBarErr = std(ca3_SDRforEachDistance{1}(2,:));
        ca3_aggregateSDR = zeros(0, length(ca3_SDRforEachDistance));
        ca3_weightAggregate = zeros(0, length(ca3_SDRforEachDistance));
        ca3_aggregateSDRError = zeros(0, length(ca3_SDRforEachDistance));
        ca3_aggregateWeightedSDRError =zeros(0, length(ca3_SDRforEachDistance));
        for k = 1: length(ca3_SDRforEachDistance)
            ca3_aggregateSDR(k) = mean(ca3_SDRforEachDistance{k}(1,:));
            ca3_weightAggregate(k) = mean(ca3_SDRforEachDistance{k}(2,:));
            ca3_aggregateSDRError(k) = std(ca3_SDRforEachDistance{k}(1,:));
            ca3_aggregateWeightedSDRError(k) = std(ca3_SDRforEachDistance{k}(2,:));
        end
    catch
        warning('Not enough cycles for SDR data')
    end


    %%% CA1 SDR

    if CA1_ON == true
        ca1_SDRforEachDistance = cell(2,length(ca1_cycleSpikeData)-1);
        if length(ca1_cycleSpikeData) >= 2
            for distance = 1:length(ca1_cycleSpikeData) - 1
                ca1_SDRforEachDistance{distance} = zeros(2,length(ca1_cycleSpikeData) - distance);

                for q = 1:length(ca1_cycleSpikeData) - distance
                    x = q;
                    y = q+distance;

                    n = unique(ca1_cycleSpikeData{x}(2,:));
                    n = n(n>50);
                    k = unique(ca1_cycleSpikeData{y}(2,:));
                    k = k(k>50);
                    r = length(union(n,k));
                    s = length(intersect(n,k));
                    n = length(n);
                    k = length(k);

                    sdr = (r-s)/(n+k);
                    sdr_adjusted = (r-s)/(n+k) - abs(n-k)/(n+k);
                    if isnan(sdr)
                        sdr = 1;
                        sdr_adjusted = 1;
                    end
                    ca1_SDRforEachDistance{distance}(1,q) = sdr;
                    ca1_SDRforEachDistance{distance}(2,q) = sdr_adjusted;
                end


            end
            ca1_sdrBar = mean(ca1_SDRforEachDistance{1}(1,:));
            ca1_weightSDRBAr = mean(ca1_SDRforEachDistance{1}(2,:));
            ca1_sdrBarErr = std(ca1_SDRforEachDistance{1}(1,:));
            ca1_weightSdrBarErr = std(ca1_SDRforEachDistance{1}(2,:));
            ca1_aggregateSDR = zeros(0, length(ca1_SDRforEachDistance));
            ca1_weightAggregate = zeros(0, length(ca1_SDRforEachDistance));
            ca1_aggregateSDRError = zeros(0, length(ca1_SDRforEachDistance));
            ca1_aggregateWeightedSDRError = zeros(0, length(ca1_SDRforEachDistance));
            for k = 1: length(ca1_SDRforEachDistance)
                ca1_aggregateSDR(k) = mean(ca1_SDRforEachDistance{k}(1,:));
                ca1_weightAggregate(k) = mean(ca1_SDRforEachDistance{k}(2,:));
                ca1_aggregateSDRError(k) = std(ca1_SDRforEachDistance{k}(1,:));
                ca1_aggregateWeightedSDRError(k) = std(ca1_SDRforEachDistance{k}(2,:));
            end
        else
            disp('Warning: Not enough cycles for SDR data')
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Cross Network Cycle Lag Calculator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if CA1_ON == true

    timeDistance = 15;  % Define how far is considered lag and not just a different cycle

    afterCycles = zeros(0,length(ca3_listOfPeakTimes)); % = list of CA1 cycles that are subsequent ones
    timeLag = zeros(0,length(ca3_listOfPeakTimes));
    cycleLength = zeros(0,length(ca3_listOfPeakTimes)-1);
    phaseLag = zeros(0,length(ca3_listOfPeakTimes));
    comparativeSDR = zeros(2,length(ca3_listOfPeakTimes));
    ca1Index = 1;
    lastCA1Index = ca1Index;
    CA1ConcurrentPeaks = zeros(0,length(ca1_listOfPeakTimes)); % In these two we will store cycles
    CA3ConcurrentPeaks = zeros(0,length(ca3_listOfPeakTimes)); % indexes that are concurrent
    for k = 1: length(ca3_listOfPeakTimes)
        while ca1Index <= length(ca1_listOfPeakTimes) && ca1_listOfPeakTimes(ca1Index) <= ca3_listOfPeakTimes(k)
            ca1Index = ca1Index + 1;
        end
        if ca1Index ~= lastCA1Index % Prevents the ability for two subsequent cycles to belong to one
                                    % precursor cycle

            if k ~= length(ca3_listOfPeakTimes)
                cycleLength(k) = ca3_listOfPeakTimes(k+1)-ca3_listOfPeakTimes(k);
            end

            if ca1Index <= length(ca1_listOfPeakTimes)
                if ca1_listOfPeakTimes(ca1Index) - ca3_listOfPeakTimes(k) <= timeDistance
                    afterCycles(k) = ca1_listOfPeakTimes(ca1Index);
                    if lastCA1Index == ca1Index
                        disp('repeat')
                        disp(ca1Index)
                    end
                    lastCA1Index = ca1Index;
                    CA1ConcurrentPeaks(ca1Index) = ca1Index;
                    CA3ConcurrentPeaks(k) = k;
                    if ca1_listOfPeakTimes(ca1Index) - ca3_listOfPeakTimes(k) == 0
                        timeLag(k) = -1;
                    else
                        timeLag(k) = ca1_listOfPeakTimes(ca1Index) - ca3_listOfPeakTimes(k);
                    end
                end
            end
        end
    end

    % Processing to remove place holding zeros
    afterCycles = afterCycles(afterCycles~=0);
    CA1ConcurrentPeaks = CA1ConcurrentPeaks(CA1ConcurrentPeaks~=0);
    CA3ConcurrentPeaks = CA3ConcurrentPeaks(CA3ConcurrentPeaks~=0);

    % Processing to remove place holding -1
    timeLag = timeLag(timeLag~=0);
    timeLag(timeLag == -1) = 0;

    avgTimeLag = mean(timeLag);
    stdTimeLag = std(timeLag);
    avgCycleLen = mean(cycleLength);

    for k = 1:length(timeLag)
        phaseLag(k) = (2 * pi * timeLag(k)) / avgCycleLen;
    end

    avgPhaseLag = mean(phaseLag);
    stdPhaseLag = std(phaseLag);
end


%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Cross SDR Calculator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sdrOff == false;
    
    if CA1_ON == true
        % Here we will get data to plot all the different possible SDR's. A total
        % of n choose 2. 

        x = combnk(1:length(CA3ConcurrentPeaks),2);
        ca3SDRList = zeros(0,length(x));
        ca1SDRList = zeros(0,length(x));
        for k = 1:length(x)

            index1 = x(k,1);  
            index2 = x(k,2);  % cycle indexes to extract sdr

            a = CA3ConcurrentPeaks(index1);
            b = CA3ConcurrentPeaks(index2);
            c = CA1ConcurrentPeaks(index1);
            d = CA1ConcurrentPeaks(index2);



            %CA3 SDR Extraction
            distance = abs(a - b); 
            if a < b
                ca3SDR = ca3_SDRforEachDistance{distance}(a);
            else
                ca3SDR = ca3_SDRforEachDistance{distance}(b);
            end

            ca3SDRList(k) = ca3SDR;


            %CA1 SDR Extraction
            distance = abs(c - d); 
            if c < d
                ca1SDR = ca1_SDRforEachDistance{distance}(c);
            else
                ca1SDR = ca1_SDRforEachDistance{distance}(d);
            end

            ca1SDRList(k) = ca1SDR;

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Frequency Calculator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CA3 Frequency Calculation

pnHzDataList = zeros(1,ca3_numberOfExcitatoryNeuronFiles);
ca3_pnNumberOfSpikes = zeros(1,ca3_numberOfExcitatoryNeuronFiles);
for k = 1:ca3_numberOfExcitatoryNeuronFiles
    pnHzDataList(k) = length(ca3_excitatoryNeuronSpikeTimeData{k})/(totalTimeinMs/1000);
    ca3_pnNumberOfSpikes(k) = length(ca3_excitatoryNeuronSpikeTimeData{k});
end
ca3_pnAvgHz = mean(pnHzDataList);
ca3_pnStDev = std2(pnHzDataList);

inHzDataList = zeros(1,ca3_numberOfInhibitoryNeuronFiles);
ca3_inNumberOfSpikes = zeros(1,ca3_numberOfInhibitoryNeuronFiles);
for k = 1:ca3_numberOfInhibitoryNeuronFiles
    inHzDataList(k) = length(ca3_inhibitoryNeuronSpikeTimeData{k})/(totalTimeinMs/1000);
    ca3_inNumberOfSpikes(k) = length(ca3_inhibitoryNeuronSpikeTimeData{k});
end
ca3_inAvgHz= mean(inHzDataList);
ca3_inStDev = std2(inHzDataList);


%%% CA1 Frequency Calculation

if CA1_ON == true
    pnHzDataList = zeros(1,ca1_numberOfExcitatoryNeuronFiles);
    for k = 1:ca1_numberOfExcitatoryNeuronFiles
        pnHzDataList(k) = length(ca1_excitatoryNeuronSpikeTimeData{k})/(totalTimeinMs/1000);
    end
    ca1_pnAvgHz = mean(pnHzDataList);
    ca1_pnStDev = std2(pnHzDataList);

    inHzDataList = zeros(1,ca1_numberOfInhibitoryNeuronFiles);
    for k = 1:ca1_numberOfInhibitoryNeuronFiles
        inHzDataList(k) = length(ca1_inhibitoryNeuronSpikeTimeData{k})/(totalTimeinMs/1000);
    end
    ca1_inAvgHz = mean(inHzDataList);
    ca1_inStDev = std2(inHzDataList);
end


%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Jitter Calculator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CA3 Jitter Calculations

ca3_jitterPerCycle = zeros(0,length(ca3_PNcycleSpikeData));
for k = 1: length(ca3_PNcycleSpikeData)
    n = std(ca3_PNcycleSpikeData{k},0,2); % this command makes n the std for every row
    ca3_jitterPerCycle(k) = n(1,1);
end

ca3_jitterOfJitter = mean(ca3_jitterPerCycle);
ca3_jOjError = std(ca3_jitterPerCycle);


% CA1 Jitter Calculations

if CA1_ON == true
    ca1_jitterPerCycle = zeros(0,length(ca1_PNcycleSpikeData));
    for k = 1: length(ca1_PNcycleSpikeData)
        n = std(ca1_PNcycleSpikeData{k},0,2);
        ca1_jitterPerCycle(k) = n(1,1);
    end

    ca1_jitterPerCycle(isnan(ca1_jitterPerCycle)) = 0;  % Code to remove NaN
    ca1_jitterOfJitter = mean(ca1_jitterPerCycle);
    ca1_jOjError = std(ca1_jitterPerCycle);

    % Calculations to subtract jitters for the correct cycles
    ca3_jitterDifference = zeros(0,length(CA3ConcurrentPeaks));
    ca1_jitterDifference = zeros(0,length(CA3ConcurrentPeaks));
    for k = 1:length(CA3ConcurrentPeaks)
        if CA3ConcurrentPeaks(k) ~= 1 && CA1ConcurrentPeaks(k) ~= 1

            ca3_cycle = CA3ConcurrentPeaks(k) - 1;
            ca1_cycle = CA1ConcurrentPeaks(k) - 1;

            if ca3_jitterPerCycle(ca3_cycle) ~= 0
                ca3_jitterDifference(k) = ca3_jitterPerCycle(ca3_cycle);
            else
                ca3_jitterDifference(k) = -1;
            end


            if ca1_jitterPerCycle(ca1_cycle) ~= 0
                ca1_jitterDifference(k) = ca1_jitterPerCycle(ca1_cycle);
            else
                ca1_jitterDifference(k) = -1;
            end
        end
    end
end


if CA1_ON == true
    % Processing to remove place holding zeros
    ca3_jitterDifference = ca3_jitterDifference(ca3_jitterDifference~=0);
    ca1_jitterDifference = ca1_jitterDifference(ca1_jitterDifference~=0);

    % Processing to remove place holding -1
    ca3_jitterDifference(ca3_jitterDifference == -1) = 0;
    ca1_jitterDifference(ca1_jitterDifference == -1) = 0;

    jitterDifference = ca3_jitterDifference - ca1_jitterDifference;

    % Used for plotting CA1 vs CA3 Jitter
    jitterCombo = vertcat(ca3_jitterDifference,ca1_jitterDifference);
    [~, order] = sort(jitterCombo(1,:));
    jitterCombo = jitterCombo(:,order);
end


%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PNs VS Cycle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CA3 Number of PN's Per Cycle 
try
    ca3_numberOfPNsPerCycleData = zeros(3,length(ca3_cycleSpikeData));
    for k = 1: length(ca3_cycleSpikeData)
        ca3_numberOfPNsPerCycleData(1,k) = length(unique(ca3_cycleSpikeData{k}(2,:)));
        ca3_numberOfPNsPerCycleData(2,k) = ca3_jitterPerCycle(k);
    end
    [~, order] = sort(ca3_numberOfPNsPerCycleData(1,:));
    ca3_numberOfPNsPerCycleData = ca3_numberOfPNsPerCycleData(:,order);

    for k = 1:length(ca3_numberOfPNsPerCycleData)
        z = 1;
        tempListOfElements = zeros(1,1);
        tempListOfElements(1,1) = ca3_numberOfPNsPerCycleData(2,k);
        while (k+z <= length(ca3_numberOfPNsPerCycleData)) && (ca3_numberOfPNsPerCycleData(1,k) == ca3_numberOfPNsPerCycleData(1,k+z)) && (ca3_numberOfPNsPerCycleData(1,k) ~= 0)
            tempListOfElements(1,z+1) = ca3_numberOfPNsPerCycleData(2,k+z);
            ca3_numberOfPNsPerCycleData(1,k+z) = 0;
            ca3_numberOfPNsPerCycleData(2,k+z) = 0;

            z = z + 1;
        end
        if ca3_numberOfPNsPerCycleData(1,k) ~= 0
            ca3_numberOfPNsPerCycleData(2,k) = mean(tempListOfElements);
            ca3_numberOfPNsPerCycleData(3,k) = std(tempListOfElements);

        end


    end

    k = 1;
    while k <= length(ca3_numberOfPNsPerCycleData)
        if ca3_numberOfPNsPerCycleData(1,k) == 0
            ca3_numberOfPNsPerCycleData(:,k) = [];
        else
            k = k + 1;
        end
    end

catch
    warning('Not Enough Data To Calc PNs vs Cycle');
end
% CA1 Number of PN's Per Cycle 
if CA1_ON == true
    
    ca1_numberOfPNsPerCycleData = zeros(3,length(ca1_cycleSpikeData));
    for k = 1: length(ca1_cycleSpikeData)
        ca1_numberOfPNsPerCycleData(1,k) = length(unique(ca1_cycleSpikeData{k}(2,:)));
        ca1_numberOfPNsPerCycleData(2,k) = ca1_jitterPerCycle(k);
    end
    [~, order] = sort(ca1_numberOfPNsPerCycleData(1,:));
    ca1_numberOfPNsPerCycleData = ca1_numberOfPNsPerCycleData(:,order);

    for k = 1:length(ca1_numberOfPNsPerCycleData)
        z = 1;
        tempListOfElements = zeros(1,1);
        tempListOfElements(1,1) = ca1_numberOfPNsPerCycleData(2,k);
        while (k+z <= length(ca1_numberOfPNsPerCycleData)) && (ca1_numberOfPNsPerCycleData(1,k) == ca1_numberOfPNsPerCycleData(1,k+z)) && (ca1_numberOfPNsPerCycleData(1,k) ~= 0)
            tempListOfElements(1,z+1) = ca1_numberOfPNsPerCycleData(2,k+z);
            ca1_numberOfPNsPerCycleData(1,k+z) = 0;
            ca1_numberOfPNsPerCycleData(2,k+z) = 0;

            z = z + 1;
        end
        if ca1_numberOfPNsPerCycleData(1,k) ~= 0
            ca1_numberOfPNsPerCycleData(2,k) = mean(tempListOfElements);
            ca1_numberOfPNsPerCycleData(3,k) = std(tempListOfElements);

        end


    end

    k = 1;
    while k <= length(ca1_numberOfPNsPerCycleData)
        if ca1_numberOfPNsPerCycleData(1,k) == 0
            ca1_numberOfPNsPerCycleData(:,k) = [];
        else
            k = k + 1;
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Presynaptic Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if CA1_ON == true
    
    preSynapticCycleData = cell(1,length(CA3ConcurrentPeaks));
    synapticXVector = zeros(1,length(CA3ConcurrentPeaks));
    for k = 1:length(preSynapticCycleData)
        preSynapticCycleData{k} = cell(1,ca1_numberOfExcitatoryNeuronFiles);
    end

    for k = 1:length(preSynapticCycleData) % For every cycle
        ca3_cycle = CA3ConcurrentPeaks(k);
        ca1_cycle = CA1ConcurrentPeaks(k);

        if ca3_cycle ~= 1 && ca1_cycle ~= 1
            ca3_cycle = ca3_cycle - 1;
            ca1_cycle = ca1_cycle - 1;
            listOfCA3Neurons = unique(ca3_cycleSpikeData{ca3_cycle}(2,:));
            listOfCA3Neurons = listOfCA3Neurons(listOfCA3Neurons>50);
            listOfCA1Neurons = unique(ca1_cycleSpikeData{ca1_cycle}(2,:));
            listOfCA1Neurons = listOfCA1Neurons(listOfCA1Neurons>50);
            synapticXVector(k) = ca1_cycle;


            for z = 1: length(listOfCA1Neurons)

                ca1Neuron = listOfCA1Neurons(z);

                presynaptics = ca1_excitatoryNeuronReverseConnectionData{ca1Neuron-50}; 
                presynaptics = intersect(presynaptics,listOfCA3Neurons);

                preSynapticCycleData{k}{ca1Neuron-50} = presynaptics;

            end

        end
    end

    avgOverPN = zeros(2,200);

    for k = 1: 200
        totalLength = zeros(1,length(CA3ConcurrentPeaks));
        for l=1:length(CA3ConcurrentPeaks)
        totalLength(l) = length(preSynapticCycleData{l}{k});
        end
        avgOverPN(1,k) = mean(totalLength);
        avgOverPN(2,k) = std(totalLength);
    end

    avgOverCycle = zeros(2,length(CA3ConcurrentPeaks));
    for l = 1:length(CA3ConcurrentPeaks)
        totalLength = zeros(1,200);
        for k = 1:200
            totalLength(k) = length(preSynapticCycleData{l}{k});
        end
        avgOverCycle(1,l) = mean(totalLength);
        avgOverCycle(2,l) = std(totalLength);
    end

end



%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Instant Cycle Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


binSize = 10;
maxHz = 200;
minHz = 0;

ca3_listOfCyclesPerHz = zeros(1,maxHz-minHz-binSize+1);
tickString = cell(1,length(ca3_listOfCyclesPerHz)); 
m = 1;
while minHz + binSize <= maxHz
    
    numCyclesForThisRange = 0;
    for k = 1:length(ca3_listOfPeakTimes)-1  % -1 because period is defined
                                             % as the time between the
                                             % spikes so you need to have
                                             % one less iteration
        period = ca3_listOfPeakTimes(k+1) - ca3_listOfPeakTimes(k);
        frequency = inv(period/1000);
        if frequency >= minHz && frequency < minHz+binSize
            numCyclesForThisRange = numCyclesForThisRange + 1;
        end
        
        
    end
    tickString{m} = strcat(num2str(minHz),'-',num2str(minHz+binSize));
    ca3_listOfCyclesPerHz(m) = numCyclesForThisRange;
    
    minHz = minHz + 1;
    m = m + 1;
end



ca3_individualCycleHz = zeros(1,length(ca3_listOfPeakTimes)-1);
for k = 1:length(ca3_listOfPeakTimes)-1
    period = ca3_listOfPeakTimes(k+1) - ca3_listOfPeakTimes(k);
    frequency = inv(period/1000);
    ca3_individualCycleHz(k) = frequency;
end
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
if multiRun == true
    doj = getSection(run1,'Design Data');
    doj2 = getSection(run2,'Design Data');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if CA1_ON == true
    ca1_pnSpikePerCycleMean = zeros(1,length(ca1_listOfNumberOfExcitatoryNeuronSpikesPerCycle));
    ca1_inSpikePerCycleMean = zeros(1,length(ca1_listOfNumberOfInhibitoryNeuronSpikesPerCycle));
end
ca3_pnSpikePerCycleMean = zeros(1,length(ca3_listOfNumberOfExcitatoryNeuronSpikesPerCycle));
ca3_inSpikePerCycleMean = zeros(1,length(ca3_listOfNumberOfInhibitoryNeuronSpikesPerCycle));

for k = 1:200
    if CA1_ON == true
        ca1_pnSpikePerCycleMean(k) = mean(ca1_listOfNumberOfExcitatoryNeuronSpikesPerCycle{k});
    end
    ca3_pnSpikePerCycleMean(k) = mean(ca3_listOfNumberOfExcitatoryNeuronSpikesPerCycle{k});
end

for k = 1:50
    if CA1_ON == true
        ca1_inSpikePerCycleMean(k) = mean(ca1_listOfNumberOfInhibitoryNeuronSpikesPerCycle{k});
    end
    ca3_inSpikePerCycleMean(k) = mean(ca3_listOfNumberOfInhibitoryNeuronSpikesPerCycle{k});
end



%X(isnan(X)) = []
% Cleaning
clear k
clear z
clear l
clear x
clear q
clear y
clear n
clear r
clear s
clear t
clear T
clear X
clear L
clear a
clear b
clear c
clear d
clear m

clear ca3_yValue
clear ca3_yValue
clear ca3_timeIntervalFilter
clear ca3_heightFilter
clear ca3_coeff24hMA
clear ca3_cycle
clear ca3SDR

clear ca1_yValue
clear ca1_yValue
clear ca1_timeIntervalFilter
clear ca1_heightFilter
clear ca1_coeff24hMA
clear ca1_longestList
clear ca1_cycle
clear ca1SDR

clear peakTime
clear pnNeuronVoltageOfASingleDt
clear bCellVoltageOfASingleDt
clear allNeuronVoltageOfASingleDt
clear numberOfSpikes
clear pnNeuronSpikesOfASingleCycle
clear InNeuronSpikesOfASingleCycle
clear singleCycleData
clear zeroLocation
clear longestList
clear distance
clear sdr
clear sdr_adjusted
clear inHzDataList
clear pnHzDataList
clear temp
clear tempListOfElements

clear timeDistance
clear ca1Index
clear lastCA1Index
clear index1
clear index2
clear dataIgnore
clear dataLength
clear Fs
clear myfilename
clear order

clear period
clear frequency
clear numCyclesForThisRange