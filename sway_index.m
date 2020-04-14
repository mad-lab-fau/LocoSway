function estimated_sway_index = sway_index(File_name,paw_variation)
% This code estimates sway index the mouse paw position.
% The input for this is code the CatWalk run file.
% Which paws are included in the calculation can be decided by the
% paw_variation: 'all_paws', 'forepaws', or 'hindpaws'.
%
% This code is related to the publication Timotius et.al, 
% "Dynamic footprint based locomotion sway assessment in alpha-synucleinopathic
% mice using Fast Fourier Transform and Low Pass Filter", Journal of Neuroscience Methods, 2018.
%
% Example: sway_index([pwd ,'\3628_Run008.xlsx'],'all_paws')
%
% Ivanna K. Timotius (2018)

        % Sampling Frequency of the CatWalk video
        Fs = 100;   

        % Parameter of the Hamming Filter
        Ny = 20;         % Order in y-axis
        Nx = 3;          % Order in x-axis

        Fc   = 1;               % Cutoff Frequency in y-axis
        Fcx = 20;               % Cutoff Frequency in x-axis
        offset = floor(Ny/2);   % To exclude the first and last half-part of the window 

        flag = 'scale';         % Sampling Flag

        % Create the window vector for the design algorithm.
        win_y = hamming(Ny+1);
        win_x = hamming(Nx+1);
        % Calculate the coefficients using the FIR1 function.
        by  = fir1(Ny, Fc/(Fs/2), 'low', win_y, flag);  % filter design  
        bx  = fir1(Nx, Fcx/(Fs/2), 'low', win_x, flag); 

        % Read the file
        B = xlsread(File_name);

        % Estimating Sway Index:
        % Reading x- & y-position and intensity of paws for STFFT-based algorithm
        PositionyRF = B(:,4);
        PositionyRH = B(:,12);
        PositionyLF = B(:,20);
        PositionyLH = B(:,28);
        PositionxRF = B(:,3);
        PositionxRH = B(:,11);
        PositionxLF = B(:,19);
        PositionxLH = B(:,27);
        
        % Calculating the center point
        if paw_variation == 'all_paws',
            % 4 paws:
            PositionxMean = nanmean([PositionxRF,PositionxRH,PositionxLF,PositionxLH],2);
            PositionyMean = nanmean([PositionyRF,PositionyRH,PositionyLF,PositionyLH],2);
        elseif paw_variation == 'forepaws',
            % Front paws:
            PositionxMean = nanmean([PositionxRF,PositionxLF],2);
            PositionyMean = nanmean([PositionyRF,PositionyLF],2);
        elseif paw_variation == 'hindpaws',
            % Hind paws:
            PositionxMean = nanmean([PositionxRH,PositionxLH],2);
            PositionyMean = nanmean([PositionyRH,PositionyLH],2);
        end
        
        PositionNan = isnan(PositionyMean);
        PositionyPaw = PositionyMean;
        PositionxPaw = PositionxMean;
        PositionyPaw(PositionNan == 1) = [];
        PositionxPaw(PositionNan == 1) = [];
        PositionyPaw = PositionyPaw - mean(PositionyPaw);
        
        % Hamming filters in both axis
        PositionyFilt = filtfilt(by,1,PositionyPaw);
        PositionxFilt = filtfilt(bx,1,PositionxPaw);
        
        PositionyFiltoffset = PositionyFilt;
        PositionyFiltoffset(1:offset) = NaN;
        PositionyFiltoffset(end-offset+1:end) = NaN;
        
        PositionxFiltoffset = PositionxFilt;
        PositionxFiltoffset(1:offset) = NaN;
        PositionxFiltoffset(end-offset+1:end) = NaN;
        
        % Path-distance calculation
        x_j = PositionxFilt(offset+1:end-offset)-PositionxFilt(offset:end-offset-1);
        y_j = PositionyFilt(offset+1:end-offset)-PositionyFilt(offset:end-offset-1);
        path_distance = sum(sqrt(x_j.^2 + y_j.^2));
        
        % Horizontal displacement
        x_min = min(PositionxFilt(offset+1:end-offset));
        x_max = max(PositionxFilt(offset+1:end-offset));
        displacement = x_max - x_min;
        
        % Sway Index
        estimated_sway_index = path_distance./displacement;    