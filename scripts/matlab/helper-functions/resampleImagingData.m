% -------------------------------------------------------------------------
% This function takes calcium imaging data and resamples it based on the
% reference (final) frame rate, which typically is the lowest frame rate
% acquired from a group of animals/sessions
%
% For data with >1 dimensions, only the 1st dimension (rows) are resampled
% Make sure your time series is located in the 1st dimension
% -------------------------------------------------------------------------

function ResampledData = resampleImagingData(data,FrameRate,finalFrameRate)
    if FrameRate == finalFrameRate
        ResampledData = data;
    else
        switch length(size(data))
            case [1,2]
                tempData = resample(data,finalFrameRate,FrameRate);
                % store final data
                ResampledData = tempData(1:end-1,:);
            case 3
                tData = resample(data(:,:,1),finalFrameRate,FrameRate);
                newData = zeros(size(tData,1),size(data,2),size(data,3));
                for i = 1:size(data,3)
                    newData(:,:,i) = resample(data(:,:,i),...
                        finalFrameRate,FrameRate);
                end
                ResampledData = newData(1:end-1,:,:);
            case 4
                tData = resample(data(:,:,1,1),finalFrameRate,FrameRate);
                newData = zeros(size(tData,1),size(data,2),...
                    size(data,3),size(data,4));
                for j = 1:size(data,4)
                    for i = 1:size(data,3)
                        newData(:,:,i,j) = resample(data(:,:,i,j),...
                            finalFrameRate,FrameRate);
                    end
                end
                ResampledData = newData(1:end-1,:,:,:);
            otherwise
                warning(['cannot resample data with >4 dimensions, ',...
                    'returning original data'])
                ResampledData = data;
        end
    end
end