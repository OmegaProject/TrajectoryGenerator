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

function NoisyTrajectoryGenerator()
    % this function calculates the uncertainty for D and SMSS for each combination of SNR, L, SMSS, D
    % then aggregates the results by SNR and L
    
    % note: in order to change the number and the 'types' of combinations, change the i, ii, iii and iiii intervals
    % in the OMEGA's road-map is listed: 12SNR*20L*10SMSS*10D = 24'000 types
    % remember that for EACH type, Ntracks trajectories will be generated
    
    print_tracks = true;
    pathSep = '/';
    
    folder = strcat('F:', pathSep, '2014-10-23_TrajectoryGenerator_Brownian_Larry');
    
    % specify here the L, SMSS and D values to use
    %L_values    = [20, 30, 38, 49, 63, 81, 103, 132, 169, 216, 277, 354, 453, 580, 743, 1000, 10000];
	L_values    = [10000, 50000];
    %SMSS_values  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    %D_values     = [0.00005, 0.0002, 0.0008, 0.003, 0.01, 0.05, 0.2, 0.8, 3.2, 13.1, 52.4];
	D_values     = [0.01, 0.1];
    
    lengthDivisor = 10;
    
    % for a "complete" simulation use
    %L_values     = [10, 30, 38, 49, 63, 81, 103, 132, 169, 216, 277, 354, 453, 580, 743, 1000];
    %SMSS_values  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    %D_values     = [0.00005, 0.0002, 0.0008, 0.003, 0.01, 0.05, 0.2, 0.8, 3.2, 13.1, 52.4];
    
    % number of tracks to be generated in each call to 'artificialTrajectories'
    Ntracks = 50;
    
    % loop through L
    l_L = length(L_values);
    for i = 1:l_L
        L = L_values(i);
        
        % we need to sort files by name in Java
        i_placeholder  = '';
        if(i < 10)
            i_placeholder = '0';
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
        
        S_fname = strcat(S_folderName, pathSep, 'SMSS_values_L_', i_placeholder, int2str(i), '.csv');
        D_fname = strcat(D_folderName, pathSep, 'D_values_L_', i_placeholder, int2str(i), '.csv');
        
        fidS    = fopen(S_fname, 'a');
        fidD    = fopen(D_fname, 'a');
        
        % just the index to be saved in the output files
        index = 1;
        
        
        % loop through D
        l_D = length(D_values);
        for ii = 1:l_D
            D = D_values(ii);
            
            % only for debugging reasons
            str = sprintf('index L: %d of %d, D: %d of %d.', i, l_L, ii, l_D);
            disp(str);
            
            t1 = tic;
            [Ds, SMSSs] = artificialTrajectories(folder, Ntracks, D, L, lengthDivisor, pathSep, print_tracks);
            toc(t1)
            
            % save to files
            fprintf(fidS,'%d %f %f;***', [index; L; D]);
            fprintf(fidD,'%d %f %f;***', [index; L; D]);
            
            for k = 1:Ntracks
                fprintf(fidS,';%f', SMSSs{k});
                fprintf(fidD,';%f', Ds{k});
            end
            
            fprintf(fidS,'\n');
            fprintf(fidD,'\n');
            
            index = index + 1;
        end
    end
    
    % close files
    fclose(fidS);
    fclose(fidD);
end

% VG: added bias and sigma as input, Ds and SMSSs as output
% VG: removed tracks as output
function [Ds, SMSSs] = artificialTrajectories(folder, Ntracks, D, L, lengthDivisor, pathSep, print_tracks)
% Produces artificial trajectories with arbitrary diffusion-type which is
% adjusted by alpha. Tracks are saved in textfiles

    % some params
    micronPerPixel = 1;
    secPerFrame    = 1;

    % VG: removed check input
    
    if(print_tracks == true)
        snrString = sprintf('%f', L);
        folderName = strcat(folder, pathSep, 'L_', snrString);
        if(exist(folderName, 'dir') == 0)
            mkdir(folderName);
        end
        dString = sprintf('%f', D);
        subFolderName = strcat(folderName, pathSep, 'D_', dString);
        if(exist(subFolderName, 'dir') == 0)
            mkdir(subFolderName);
        end
        filename = strcat(subFolderName, pathSep, 'track_');
    end

    tracks   = cell(Ntracks,1);
    frameNos = cell(Ntracks,1);
    
    %in order to store D and SMSS
    Ds    = cell(Ntracks,1);
    SMSSs = cell(Ntracks,1);
    
    for itrack = 1:Ntracks
        t_max = L;
        frameNo = [1:1:t_max];
        Wx = zeros(1,t_max);
        Wy = zeros(1,t_max);
        frameNos{itrack} = frameNo;

        stDev = sqrt(4 * D * secPerFrame);
        for t = 2:t_max
            angle = pi * rand;
            step = stDev * randn;
            k = t - 1;
            Wx(t) = Wx(k) + step * cos(angle);
            Wy(t) = Wy(k) + step * sin(angle);
        end
        
        clear points;
        points(:,1) = Wx;
        points(:,2) = Wy;
        
        % store track
        tracks{itrack} = points;
        
        % analize resulting trajectories with getMss
        [D_sim_new, ~, mssSlope_sim_new] = getMss(points, lengthDivisor, micronPerPixel, secPerFrame);
             
        % store D and SMSS of the points
        Ds{itrack}    = D_sim_new(3);
        SMSSs{itrack} = mssSlope_sim_new;
    end
    
    % file saving method
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
    end
end