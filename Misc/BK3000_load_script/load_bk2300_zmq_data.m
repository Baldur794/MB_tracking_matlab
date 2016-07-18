% LOAD_BK2300_ZMQ_DATA - load data from the BK2300 scanner stored through
% ZMQ
%
% This function loads a binary file stored by the ZmqReceiver application.
% It is up to the user to specify the image dimensions, sample type, etc.
%
% Usage:
%   img = load_bk2300_zmq_data(filename, n_lines, n_samples,
%                              sample_type, iq)
%
% Inputs:
%   filename - string, the filename to load
%   n_lines - scalar, number of lines in the image
%   n_samples - scalar, number of samples along a line
%   sample_type - string, one of [u]int{8,16,32,64}
%   iq - logical, flag indicating IQ samples (1) or simple sampling (0)
%
% Outputs:
%   img - matrix [n_samples n_lines], the data in the file, converted to
%         double precision floating point
%
% History:
%   v1.0.0, 2016-03-31, mbst, Initial version

function img=load_bk2300_zmq_data(filename, n_lines, n_samples, sample_type, iq)

    % Open the file
    fid = fopen(filename);
    if fid==-1
        error('File not found: %s', filename);
    end

    % Calculate bytes per sample
    if strcmp(sample_type, 'int8') || strcmp(sample_type, 'uint8')
        bytes_per_sample = 1;
    elseif strcmp(sample_type, 'int16') || strcmp(sample_type, 'uint16')
        bytes_per_sample = 2;
    elseif strcmp(sample_type, 'int32') || strcmp(sample_type, 'uint32')
        bytes_per_sample = 4;
    elseif strcmp(sample_type, 'int64') || strcmp(sample_type, 'uint64')
        bytes_per_sample = 8;
    else
        error('Unsupported sample type %s', sample_type);
    end

    % Read the contents
    if iq
        n_raw_samples=2*n_samples;
    else
        n_raw_samples=n_samples;
    end
    raw = double(fread(fid, [n_raw_samples n_lines], sample_type));
    
    % Close the file
    fclose(fid);
    
    % Make sure the amount of data matches the expected amount of data
    data_expected = n_lines*n_samples*bytes_per_sample;
    if iq
        data_expected = 2*data_expected;
    end
    bytes_read=numel(raw)*bytes_per_sample;
    if bytes_read ~= data_expected
        error('File contains %d bytes, expected %d bytes', bytes_read, data_expected);
    end

    % Produce the output image
    if iq
        img = raw(1:2:end,:)+1i*raw(2:2:end,:);
    else
        img = raw;
    end
    
