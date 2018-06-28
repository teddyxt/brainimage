
function [prod] = dFoF(stack,w,nframe)
% w: sliding window width; nframe: total frames of the input video
% Function accepts an image file and frame #s that bound a baseline period
% and returns a stack of dF/F images.
% The ?F/F of individual ROI traces was calculated by subtracting the mean
% of the lowest 25% of values within a w-frame sliding window from each 
% value at the center of this window, divided by the mean of the lowest 
% 25% of the window.

        % dFoF interval
        t_0 = floor(w/2);
        t_n = nframe - (w - t_0) + 1;
        fprintf(1,'dFoF starting frame %d\n', t_0);
        fprintf(1,'dFoF end frame %d\n', t_n);
        % initiate prod
        prod = zeros(size(stack,1)*size(stack,2),t_n);
        % number of lowest 25% of values
        n_25 = floor(w * 0.25);
        % batch
        batchsize = 1000;
        nbatch = floor((t_n - t_0 + 1)/batchsize);
        % convert video from 3d to 2d
        stackflat = single(reshape(stack,[],size(stack,3)));
        % calculate dF/F
        fprintf(1,'calculating dF/F...\n');
        
        tic
        for b = 1:nbatch
            fprintf(1,'    frame %d\n', b*batchsize);
            % starting/end frame in each batch
            fstart = (b-1)*batchsize+t_0;
            fend = b*batchsize+t_0-1;
            parfor f= fstart:fend
                % lowest 25% window
                substack = stackflat(:,(f-t_0+1):(f+w-t_0));
                %mean_lowest = zeros(size(substack,1),1);
                % find the lowest 25% for each pixel 
                mink = sort(substack,2);
                mean_lowest=single(mean(mink(:,1:n_25),2))
                % subtract mean over mean
                prod(:,f) = (stackflat(:,f) - mean_lowest)./mean_lowest*100;
            end
        end
        % do the same for the last batch
        fprintf(1,'    frame %d\n', t_n);
        parfor f = (nbatch*batchsize+t_0-1):t_n
            % lowest 25% window
            substack = stackflat(:,(f-t_0+1):(f+w-t_0));
            %mean_lowest = zeros(size(substack,1),1);
            % find the lowest 25% for each pixel 
            mink = sort(substack,2);
            mean_lowest=single(mean(mink(:,1:n_25),2))
            % subtract mean over mean
            prod(:,f) = (stackflat(:,f) - mean_lowest)./mean_lowest*100;
        end
        toc
        clear stackflat substack result mean_lowest mink
        % discard 1:t_0 (empty) and convert to 3d
        prod = reshape(prod(:,t_0:t_n),size(stack,1),size(stack,2));
        fprintf(1,'dF/F is completed');
end
