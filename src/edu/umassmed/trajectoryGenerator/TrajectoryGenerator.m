% Vanni Galli, 2012
% Alex Rigano, 2013
% ErrorEstimation v1.0, 2013-5-31 
%
% This code calls a modified version of the 'artificialTrajectories2' function by Ivo Sbalzarini.
% Changes in 'artificialTrajectories' are marked with 'VG'.
%
% Vanni modifies
% 2012-10-16, v0.2: fixed some parts of the code
% 2012-10-17, v0.2: new fixes on random (!)
% 2012-10-18, v0.3: general refactoring, removed useless array passed to artificialTrajectories
% 2012-10-23, v0.4: changed the way we are adding noise, as per Ivo's e-mail on 2012-10-23

function NoisyTrajectoryGenerator(snrFrom, snrTo)
    % this function calculates the uncertainty for D and SMSS for each combination of SNR, L, SMSS, D
    % then aggregates the results by SNR and L
    
    % note: in order to change the number and the 'types' of combinations, change the i, ii, iii and iiii %intervals
    % in the OMEGA's road-map is listed: 12SNR*20L*10SMSS*10D = 24'000 types
    % remember that for EACH type, Ntracks trajectories will be generated
	
    print_tracks = true;
	pathSep = '/';
	
	%folder = strcat('F:', pathSep, '2014-10-06_TrajectoryGeneratorValidation_NoNoise');
	%distriFolder = strcat(folder, pathSep, '2014-01-24_Image for bias and sigma estimation');
	folder = strcat(pathSep, 'storage', pathSep, 'raid10', pathSep, 'Matlab', pathSep, '2015-07-29_TrajectoryFullSimulation');
	distriFolder = strcat(folder, pathSep, '2014-01-24_Image for bias and sigma estimation');
    
    % specify here the L, SMSS and D values to use
	%SNR_values = [1.291059];
    %L_values    = [20, 30, 38, 49, 63, 81, 103, 132, 169, 216, 277, 354, 453, 580, 743, 1000];
    %SMSS_values  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    %D_values     = [0.00005, 0.0002, 0.0008, 0.003, 0.01, 0.05, 0.2, 0.8, 3.2, 13.1, 52.4];
    %SNR_values   = [1.291059,   3.494379,   11.632132,  31.306549];
    %signal_values = [15.00,    28.73,  154.70,  1000.00];  
    lengthDivisor = 3;
    
    % for a "complete" simulation use
    L_values     = [10, 30, 38, 49, 63, 81, 103, 132, 169, 216, 277, 354, 453, 580, 743, 1000];
    SMSS_values  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    D_values     = [0.00005, 0.0002, 0.0008, 0.003, 0.01, 0.05, 0.2, 0.8, 3.2, 13.1, 52.4];
	%D_values      = [0.001, 0.005,  0.0075, 0.01,   0.02,   0.05,   0.1,    0.5,    0.75,   1,  2,  5,  20];
    % figure 4, Sbalzarini and Koumotsakos, 2005
    SNR_values   = [1.291059, 1.990510, 2.846111, 3.494379, 4.556798, 6.516668, 8.832892, 11.632132, 15.067460, 19.326731, 24.642859, 31.306549];
    signal_values = [15.00,   18.58,  23.90,   28.73,  38.10,   60.80,   97.00, 154.70,  246.60,  393.30,  627.10,  1000.00];
    %Bias_values  = [0.6135    0.3918    0.2131    0.1435    0.0833    0.0417    0.0267    0.0182     0.0130     0.0100     0.0084     0.0071];
    %Sigma_values = [0.4275    0.3484    0.2449    0.2043    0.1677    0.1121    0.0887    0.0681     0.0545     0.0458     0.0383     0.0345];
    
    % number of tracks to be generated in each call to 'artificialTrajectories'
    Ntracks = 1000;
        
    % loop through SNR
    for i = snrFrom:snrTo
        SNR   = SNR_values(i);
        signal = signal_values(i);
        distriSignalFolder = strcat(distriFolder, pathSep, sprintf('Special_100_%4.2f_10.0_16_Bits_Grey_1_null_test', signal));
        distriFile = strcat(distriSignalFolder, pathSep, 'P2PDistanceCalculator_Results_distrib.txt');
        
        % loop through L
        l_L = length(L_values);
        for ii = 1:l_L
            L = L_values(ii);
            
            % we need to sort files by name in Java
            i_placeholder  = '';
            ii_placeholder = '';
            if(i < 10)
                i_placeholder  = '0';
            end
            if(ii < 10)
                ii_placeholder = '0';
            end 
            
            % open files
            S_folderName = strcat(folder, pathSep, 'SMSS');
            if(exist(S_folderName, 'dir') == 0)
                mkdir(S_folderName)
            end
            D_folderName = strcat(folder, pathSep, 'D');
            if(exist(D_folderName, 'dir') == 0)
                mkdir(D_folderName)
            end
            
            S_fname = strcat(S_folderName, pathSep, 'SMSS_values_SNR_', i_placeholder, int2str(i), '_L_', ii_placeholder, int2str(ii), '.csv');
            D_fname = strcat(D_folderName, pathSep, 'D_values_SNR_',       i_placeholder, int2str(i), '_L_', ii_placeholder, int2str(ii), '.csv');
            %S_fname = strcat('C:\Users\galliva\Desktop\tests\SMSS_values_SNR_', int2str(i), '_L_', int2str(ii), '.csv');
            %D_fname = strcat('C:\Users\galliva\Desktop\tests\D_values_SNR_',    int2str(i), '_L_', int2str(ii), '.csv');
            
            [fidS,msg1] = fopen(S_fname, 'a')
            [fidD,msg2] = fopen(D_fname, 'a')
            
            % just the index to be saved in the output files
            SNR_index = 1;

            % loop through SMSS
            l_SMSS = length(SMSS_values);
            for iii = 1:l_SMSS
                SMSS = SMSS_values(iii);
                
                % loop through D
                l_D = length(D_values);
                for iiii = 1:l_D
                    D = D_values(iiii);
                    
                    % SMSS array used in polyval later
                    array_SMSS(1:Ntracks)= SMSS;
                    
                    % only for debugging reasons
                    str = sprintf('index SNR: %d of %d, L: %d of %d, SMSS: %d of %d, D: %d of %d.', i, snrTo, ii, l_L, iii, l_SMSS, iiii, l_D);
                    disp(str);
                    
                    t1 = tic;  
                    [Ds, SMSSs] = artificialTrajectories(folder, SNR, Ntracks, array_SMSS, D, L, lengthDivisor, distriFile, pathSep, print_tracks);
                    toc(t1)
                    
                    % save to files
                    fprintf(fidS,'%d %f %f %f %f;***', [SNR_index; SNR; L; SMSS; D]);
                    fprintf(fidD,'%d %f %f %f %f;***', [SNR_index; SNR; L; SMSS; D]);
                    
                    for k = 1:Ntracks
                        fprintf(fidS,';%f', SMSSs{k});
                        fprintf(fidD,';%f', Ds{k});
                    end
                    
                    fprintf(fidS,'\n');
                    fprintf(fidD,'\n');
                    
                    SNR_index = SNR_index + 1;
                end
            end
            
            % close files
            fclose(fidS);
            fclose(fidD);
        end
    end
end

% VG: added bias and sigma as input, Ds and SMSSs as output
% VG: removed tracks as output
function [Ds, SMSSs] = artificialTrajectories(folder, SNR, Ntracks, mssSlope, D2, L, lengthDivisor, distriFile, pathSep, print_tracks)
% Produces artificial trajectories with arbitrary diffusion-type which is
% adjusted by alpha. Tracks are saved in textfiles

    % some params
    micronPerPixel = 1;
    secPerFrame    = 1;

    % VG: removed check input
    
    if(print_tracks == true)
        snrString = sprintf('%f', SNR);
        snrString = strrep(snrString, '.', '-');
        folderName = strcat(folder, pathSep, 'tracks_', snrString);
        if(exist(folderName, 'dir') == 0)
            mkdir(folderName);
        end
        dString = sprintf('%f', D2);
        dString = strrep(dString, '.', '-');
        mssString = sprintf('%f', mssSlope(1));
        mssString = strrep(mssString, '.', '-');
        subFolderName = strcat(folderName, pathSep, 'L_' , int2str(L), '_SMSS_', mssString, '_D_', dString);
        if(exist(subFolderName, 'dir') == 0)
            mkdir(subFolderName);
        end
        filename = strcat(subFolderName, pathSep, 'track_');
    end
    
    % use LUT instead
    % these values are the result of a test
    P = [8.22753559071985 -12.81933944852967 6.18021769632377 1.27901146576493 -0.00191148801039];
    alpha = polyval(P,mssSlope);

    gamma = sqrt(pi);

    tracks   = cell(Ntracks,1);
    foos     = cell(Ntracks,1);
    frameNos = cell(Ntracks,1);

    % VG: the noisy track
    tracks_new = cell(Ntracks,1);
    
    % VG: removed unused cntFile = 0;
    
    % VG: in order to store D and SMSS
    Ds    = cell(Ntracks,1);
    SMSSs = cell(Ntracks,1);
    
    distri_data = importdata(distriFile, '\t');
    data_row = size(distri_data, 1);
    
    for itrack = 1:Ntracks
        theta_x = rand(57,1)*2*pi;
        theta_y = rand(57,1)*2*pi;

        H = alpha(itrack)/2;

        t_max   = L;
        frameNo = [1:1:t_max];
        Wx = zeros(1,t_max);
        Wy = zeros(1,t_max);
        % VG: the noisy points
        Rx = zeros(1,t_max);
        Ry = zeros(1,t_max);
        
        foo = zeros(1,t_max);

        foos{itrack} = foo;
        frameNos{itrack} = frameNo;

        for t = 1:t_max
            Wx(t) = 0;
            Wy(t) = 0;

            for n = -8:1:48
                Wx(t) = Wx(t) + ...
                    (cos(theta_x(n+9))-cos((gamma^n)*2*pi*t/t_max +...
                    theta_x(n+9)))/gamma^(n*H);
                Wy(t) = Wy(t) + ...
                    (cos(theta_y(n+9))-cos((gamma^n)*2*pi*t/t_max +...
                    theta_y(n+9)))/gamma^(n*H);
            end
        end

        
        % VG: 2) re-scale trajectory with ground-truth D
        clear points;
        points(:,1) = Wx;
        points(:,2) = Wy;
        
        % check D and scale
        [D_sim(itrack,:) ,mss_sim(itrack,:),mssSlope_sim(itrack)] = getMss(points, lengthDivisor, micronPerPixel, secPerFrame);
        
        D2_sim(itrack) = D_sim(itrack,3);
        points = points * sqrt(D2/D2_sim(itrack));
        
        % store track
        tracks{itrack} = points;
   
        
        % VG: 3) add noise
        clear noisy_points
        
        for t = 1:t_max
            % VG: generate values from a normal distribution with mean (point_position + bias) and standard deviation sigma
            index = unidrnd(data_row);
            noisy_points(t,1) = points(t,1) + distri_data(index, 1); 
            noisy_points(t,2) = points(t,2) + distri_data(index, 2); 
        end
       
        
        % VG: 4) analize resulting scaled, noisy trajectory with getMss
        %[D_sim_new, ~, mssSlope_sim_new] = getMss(noisy_points, lengthDivisor, micronPerPixel,secPerFrame);
        [D_sim_new, ~, mssSlope_sim_new] = getMss(points, lengthDivisor, micronPerPixel,secPerFrame);

        % store track
        %tracks_new{itrack} = noisy_points;
        tracks_new{itrack} = points;
        
             
        % VG: store D and SMSS of the noisy points
        Ds{itrack}    = D_sim_new(3);
        SMSSs{itrack} = mssSlope_sim_new;
    end
    
    % VG: changed the file saving method
    if(print_tracks == true)
        % save track
        for itrack = 1:Ntracks 
            fname = strcat(filename, int2str(itrack), '.out');
            fid   = fopen(fname, 'w');

            t  = tracks{itrack};
            fr = frameNos{itrack};

            fprintf(fid,'%d\t%4.8f\t%4.8f\n', [fr; t(:,1)'; t(:,2)']);

            fclose(fid);
        end
        % save the noisy track
        for itrack = 1:Ntracks 
            fname = strcat(filename, int2str(itrack), '_noise.out');
            fid   = fopen(fname, 'w');

            t  = tracks_new{itrack};
            fr = frameNos{itrack};

            fprintf(fid,'%d\t%4.8f\t%4.8f\n', [fr; t(:,1)'; t(:,2)']);

            fclose(fid);
        end
    end
end