function EPosNorm2 = sway_STFFT_position(File_name,paw_variation)
% This code assess sway based the mouse paw position using FFT.
% The input for this is code the CatWalk run file.
% Which paws are included in the calculation can be decided by the
% paw_variation: 'all_paws', 'forepaws', or 'hindpaws'.
%
% This code is related to the publication Timotius et.al, 
% "Dynamic footprint based locomotion sway assessment in alpha-synucleinopathic
% mice using Fast Fourier Transform and Low Pass Filter", Journal of Neuroscience Methods, 2018.
%
% Example: sway_STFFT_position([pwd ,'\3628_Run008.xlsx'],'all_paws')
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
                
        % Calculating the center point
        if paw_variation == 'all_paws',
            % 4 paws:
            PositionMean = nanmean([PositionRF,PositionRH,PositionLF,PositionLH],2);
        elseif paw_variation == 'forepaws',
            % Front paws:
            PositionMean = nanmean([PositionRF,PositionLF],2);
        elseif paw_variation == 'hindpaws',
            % Hind paws:
            PositionMean = nanmean([PositionRH,PositionLH],2);
        end
        
        % Frequency analysis from center position:
        PositionNan = isnan(PositionMean);
        PositionPaw = PositionMean;
        PositionPaw(PositionNan == 1) = [];
        % 1. Substracting by its mean value to eliminate the effect of
        % mid-path
        PositionPaw = PositionPaw - mean(PositionPaw);
        length_PositionPaw = length(PositionPaw);
        % 2. Absolute value of STFFT
        number_window_P = 1 + floor((length_PositionPaw-size_fft_window)./size_fft_overlap);
        for j = 1:number_window_P,
            win_start_P = (j-1)*size_fft_overlap+1;
            win_signal_P = PositionPaw(win_start_P:win_start_P+size_fft_window-1);
            win_fft_PositionPaw(j,:) = abs(fft(win_signal_P));
        end
        % 3. Mean from all the window
        fft_PositionPaw_ave = mean(win_fft_PositionPaw,1);
        % 4. Calculating the Energy
        EPos2 = sum(fft_PositionPaw_ave(1:2).^2);       %E_s
        EPos35 = sum(fft_PositionPaw_ave(3:4).^2);      %E_n
        EPosNorm2 = EPos2./EPos35;                   %S_p