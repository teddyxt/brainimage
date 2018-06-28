%% This script is a rough sketch of how I would like to batch process images from the macroscope
addpath(genpath('~/Documents/MATLAB/'))

clear all
clc

% find all files for 1803140
[mouseDirs, imFiles]=findImDirsFiles('/media/westonlab/3C0A970335974DE3/imagingDataWillieTobin/',...
    1803140);

% remove regIms from the file list
imFiles(find(contains(imFiles,'regIms')))=[];


%% pull in a recent red stack, identify pixel locations (sensor space) over the
%  infected patch

% Calculate a mask based a recent red image stack
redStackPath='/media/westonlab/3C0A970335974DE3/imagingDataWillieTobin/09-Apr-2018/1803140/bestRed/Acq_1/Acq_1_MMStack_Pos0.ome.tif';
redStackParms=load('/media/westonlab/3C0A970335974DE3/imagingDataWillieTobin/09-Apr-2018/1803140/bestRed/Acq_1/acqParms.mat');
redStackInfo=imfinfo(redStackPath);

tic
% reading and filtering
parfor i = 1:length(redStackInfo)
    redStack(:,:,i) = imgaussfilt(double(imread(redStackPath,'Index',i,'Info',redStackInfo)),2);
end
toc

%take the median of the stack
redStackMed=median(redStack,3);


%turn it into a mask
BW1=imbinarize(redStackMed);
% b=strel('disk',5);
% BW2=imdilate(BW1,b);


%Put coords in terms of sensor pixel locations
[maskXCoords, maskYCoords]=find(BW1);
maskXCoords=maskXCoords+round(redStackParms.acqParms.roiY)-1;
maskYCoords=maskYCoords+round(redStackParms.acqParms.roiX)-1;

%% Build filters
Fs=40;
highpassCutoff = 0.01; % Hz
heartbeatBandStop = [9 14];
[b100s, a100s] = butter(2, highpassCutoff/(Fs/2), 'high');
[bHeart, aHeart] = butter(2, heartbeatBandStop/(Fs/2), 'stop');

%% Ananlyze all sessions for this animal

%group files by session
filesBySession=groupFilesBySession(imFiles);

%structure to hold products of session analysis
%sessionSummaries=struct();
load('/home/westonlab/Documents/MATLAB/willieT/imageAnalysisCode/macroscopicImages/1803140/pipelineScratch/sessionSummaries.mat')

%Add loop over sessions here
fprintf(1, 'Loop over sessions\n');
for sess=[length(sessionSummaries):length(filesBySession)]
    fprintf(1, 'The %d session: \n', sess);
    ssCatIms=[];
    counter=1;
    
    %loop over files in session
    for file=1:numel(filesBySession{sess})
        
        %Eventually I will loop over files here
        curFile=filesBySession{sess}{file};
        
        info = imfinfo(curFile); % open file info
        
        F = length(info); % number of frame
        
        I = zeros(info(1).Height,info(1).Width,F,'uint16'); % datamatrix
        tic
        % reading and filtering
        parfor i = 1:F
            I(:,:,i) = imgaussfilt(double(imread(curFile,'Index',i,'Info',info)),1);
        end
        toc
        
        tic
        for j=1:size(I,3)
            %Subsample the data
            ssCatIms(:,:,counter)=I(1:5:end,1:5:end,j);
            counter=counter+1;
        end
        toc
        
        clear I
    end
    
    %Calculate dF/F
    w = 50; % sliding window length
    ssCatIms = single(imstack(:,:,1:100));
    dFoFStack=dFoF(ssCatIms,w,size(ssCatIms,3));
    dFoFStack2=dFoFWT(ssCatIms,w,size(ssCatIms,3));
    
    %get slashes for curFil
    cFSlashes=regexp(curFile,'/');
    
    %load parms for the last file in the session to get roi location on sensor
    curFileParms=load([curFile(1:cFSlashes(end)),'acqParms.mat']);
    
    curMaskXs=maskXCoords-round(curFileParms.acqParms.roiY);
    curMaskYs=maskYCoords-round(curFileParms.acqParms.roiX);
    
    curMaskXs=round(curMaskXs/5);
    curMaskYs=round(curMaskYs/5);
    
    seedPixel=[round(mean(curMaskXs)), round(mean(curMaskYs))];
    seedPixelInds=sub2ind([size(dFoFStack,1), size(dFoFStack,2)],...
        repmat([seedPixel(1)-1:seedPixel(1)+1], [3, 1]) ,...
        repmat([seedPixel(2)-1:seedPixel(2)+1], [3 1])');
    
    dFoFVector=reshape(dFoFStack,size(dFoFStack,1)*size(dFoFStack,2),size(dFoFStack,3));
    
    %Filter the dF stack
    filtDf=zeros(size(dFoFVector));
    
    parfor p=1:size(dFoFVector,1)
        
        % dP=detrend(dFoFVector(p,:),'linear');
        %fHeart=filtfilt(bHeart,aHeart,dP);%,[],2);
        
        fHeart=filtfilt(bHeart,aHeart,dFoFVector(p,:));
        filtDf(p,:)=detrend(filtfilt(b100s,a100s,fHeart),'linear');%,[],2);
       
    end
    
    % correlation map
    CM = corr(squeeze(mean(filtDf(seedPixelInds,:)))',filtDf');
    CM = reshape(CM,size(dFoFStack,1),size(dFoFStack,2));
    
    sessionSummaries(sess).SPCM=CM;
    sessionSummaries(sess).maxDFoF=max(dFoFStack,[],3);
    sessionSummaries(sess).stdDFoF=std(dFoFStack,[],3);
    
    %the differential of dFoF stack
    dFilt=reshape(diff(filtDf,[],2),...
        size(dFoFStack,1),size(dFoFStack,2),size(dFoFStack,3)-1);
    
    sessionSummaries(sess).diffMax=max(dFilt,[],3);
    
    meanPkProm=zeros(1,size(filtDf,1));
    stdPkProm=zeros(1,size(filtDf,1));
    meanPkWidth=zeros(1,size(filtDf,1));
    stdPkWidth=zeros(1,size(filtDf,1));
    pkNum=zeros(1,size(filtDf,1));
    
    parfor p=1:size(filtDf,1)
        
        [pkVals, pkLocs, pkWidths, pkProms]=findpeaks(filtDf(p,:),'MinPeakProminence',10);
        
        meanPkProm(p)=mean(pkProms);
        stdPkProm(p)=std(pkProms);
        meanPkWidth(p)=mean(pkWidths);
        stdPkWidth(p)=std(pkWidths);
        pkNum(p)=numel(pkLocs);
        
    end
    sessionSummaries(sess).meanPkProm=meanPkProm;
    sessionSummaries(sess).stdPkProm=stdPkProm;
    sessionSummaries(sess).meanPkWidth=meanPkWidth;
    sessionSummaries(sess).stdPkWidth=stdPkWidth;
    sessionSummaries(sess).pkNum=pkNum;
    
    % apply pca to filtDf
    % Navgframes: number of consecutive frames taken into avg
    % NSVDcomp: number of svd comp want to keep
    Navgframes = 5;
    NSVDcomp = 100;
    ops = setSVDParams(filtDf,Navgframes,NSVDcomp);
    [U, Sv, V, totalVar] = get_svdcomps(ops,filtDf);
    sessionSummaries(sess).U=U;
    sessionSummaries(sess).Sv=Sv;
    sessionSummaries(sess).V=V;
    sessionSummaries(sees).totalVar= totalVar;
    
    %% save
    save('/home/westonlab/Documents/MATLAB/willieT/imageAnalysisCode/macroscopicImages/1803140/pipelineScratch/sessionSummaries.mat',...
        'sessionSummaries', '-v7.3')
    
    clear dFoFStack dFoFVector dFilt filtDf U Sv V ops
    
end

%% svd reconstruction

% View components
% extract U, Sv, V, totalVar
U = sessionSummaries(sess).U;
V = sessionSummaries(sess).V;
Sv = sessionSummaries(sess).Sv;
totalVar = sessionSummaries(sess).totalVar;
Fs = 40;
svdViewer(U, Sv, V, Fs, totalVar);

U0 = U(:,:,1:ncomps); 
V0 = V(1:ncomps,:); 
frameofRecon = 201:300; % define the frames you want to reconstruct
frameRecon = double(svdFrameReconstruct(U0,V0(:,frameofRecon)));
stackRecon = bsxfun(@plus, frameRecon, single(ops.mimg));
stackRecon = uint16(stackRecon);

