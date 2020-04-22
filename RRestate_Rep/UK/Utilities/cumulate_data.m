function arrData2 = cumulate_data(arrData, intDimCum, intDim2, vecCum)
% Turns array of one-period data (observations or forecasts) into array with
% mixed one-period and cumulative data. Most of the columns (j) are unchanged: 
%       [x(1,j,:)   x(2,j,:)          ...     x(T,j,:)] ,
% but those identified by vecCum (j*) are cumulated along the 1st
% dimension:
%       [x(1,j*,:)   x(1,j*,:)+x(2,j*,:) ...     x(1,j*,:)+...+x(T,j*,:)]
%
% Works on 2 and 3 dimensional arrays, input and output have the same
% dimensions. Dimension 1, along which we cumulate, is supposed to be time
% (T).

% Number of dimensions:
numDims = length(size(arrData));

% Dimensions:
[d1 d2 d3] = size(arrData);
% Where d1 is time and eiter d2 or d3 is the dimension along which the data
% must be cumulated.

% Size of the array along the cum dimension:
intSize = size(arrData, intDimCum);

% Create new array
arrData2 = arrData;


if intDimCum == 1 & intDim2 == 2
    arrData2(:, vecCum, :) = cumsum(arrData(:, vecCum, :), 1);

    % Adjust annualization
    arrData2(:, vecCum, :) = arrData2(:, vecCum, :) ./ repmat([1:d1]', [1 length(vecCum)]);
    
elseif intDimCum == 2 & intDim2 == 3
    arrData2(:, :, vecCum) = cumsum(arrData(:, :, vecCum), 2);
    
	% Adjust annualization:
    arrData2(:, :, vecCum) = arrData2(:, :, vecCum) ./ repmat([1:d2], [d1 1 length(vecCum)]) ;
    
else
    ('Error: ...')
end

