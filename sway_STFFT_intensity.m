function  EPosIntNorm2 = sway_STFFT_intensity(File_name,paw_variation)
% This code assess sway based the mouse paw position and intensity using FFT.
% The input for this is code the CatWalk run file.
% Which paws are included in the calculation can be decided by the
% paw_variation: 'all_paws', 'forepaws', or 'hindpaws'.
%
% This code is related to the publication Timotius et.al, 
% "Dynamic footprint based locomotion sway assessment in alpha-synucleinopathic
% mice using Fast Fourier Transform and Low Pass Filter", Journal of Neuroscience Methods, 2018.
%
% Example: sway_STFFT_intensity([pwd ,'\3628_Run008.xlsx'],'all_paws')
%
% Ivanna K. Timotius (2018)

        % STFFT parameter
        size_fft_window = 64; 
        size_fft_overlap = 5;

        % Read the file
        B = xlsread(File_name);
                
        % Reading y-position and intensity of paws for STFFT-based algorithm
        PositionRF = B(:,4);
        PositionRH = B(:,12);
        PositionLF = B(:,20);
        PositionLH = B(:,28);
        
        IntensityRF = B(:,10)-B(:,8);           % MeanIntensity - MinIntensity (for normalisation)
        IntensityRH = B(:,18)-B(:,16);          % MinIntensity is different due to green intensity threshold
        IntensityLF = B(:,26)-B(:,24);
        IntensityLH = B(:,34)-B(:,32);
        
        PositionIntensityRF = PositionRF.*IntensityRF;
        PositionIntensityRH = PositionRH.*IntensityRH;
        PositionIntensityLF = PositionLF.*IntensityLF;
        PositionIntensityLH = PositionLH.*IntensityLH;
        
        % Calculating the center point
        if paw_variation == 'all_paws',
            % 4 paws:
            PositionIntensityMean = nanmean([PositionIntensityRF,PositionIntensityRH,PositionIntensityLF,PositionIntensityLH],2);
        elseif paw_variation == 'forepaws',
            % Front paws:
            PositionIntensityMean = nanmean([PositionIntensityRF,PositionIntensityLF],2);
        elseif paw_variation == 'hindpaws',
            % Hind paws:
            PositionIntensityMean = nanmean([PositionIntensityRH,PositionIntensityLH],2);
        end
                
       % Freq analysis from Position*Intensity:
        PositionIntensityNan = isnan(PositionIntensityMean);
        PositionIntensityPaw = PositionIntensityMean;
        PositionIntensityPaw(PositionIntensityNan == 1) = [];
        % 1. Substracting by its mean value to eliminate the effect of
        % mid-path
        PositionIntensityPaw = PositionIntensityPaw - mean(PositionIntensityPaw);
        length_PositionIntensityPaw = length(PositionIntensityPaw);
        % 2. Absolute value of STFFT
        number_window_PI = 1 + floor((length_PositionIntensityPaw-size_fft_window)./size_fft_overlap);
        for j = 1:number_window_PI,
            win_start_PI = (j-1)*size_fft_overlap+1;
            win_signal_PI = PositionIntensityPaw(win_start_PI:win_start_PI+size_fft_window-1);
            win_fft_PositionIntensityPaw(j,:) = abs(fft(win_signal_PI));
        end
        % 3. Mean from all the window
        fft_PositionIntensityPaw_ave = mean(win_fft_PositionIntensityPaw,1);
        % 4. Calculating the Energy
        EPosInt2 = sum(fft_PositionIntensityPaw_ave(1:2).^2);   %E_s
        EPosInt35 = sum(fft_PositionIntensityPaw_ave(3:4).^2);  %E_n
        EPosIntNorm2 = EPosInt2./EPosInt35;                     %S_i