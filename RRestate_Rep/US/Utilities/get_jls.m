function jointLogScore = get_jls(A,B,horizons)
% Inputs: A is 3D matrix where dimension 3 indexes variables
%         B is matrix of observed values at horizon h where columns
%         indicate variables
%         horizons is vector of horizons to store values

% Output: jointLogScore is vector of joint log scores for data 
%%
% initialise vector of log-scores
jointLogScore=NaN(1,length(horizons));



% loop over horizons to get log-score for each predictive density

for i = 1:length(horizons)
    
   ii=horizons(i);
   
   pathsA=squeeze(A(:,ii,:));
   ObsB=single(B(ii,:));
   
   % create grid
   [~,d]=size(pathsA); ng=25;
  
%  for kk=1:d
%     grid(:,kk)=linspace(min(pathsA(:,kk))-5*std(pathsA(:,kk)),max(pathsA(:,kk))+5*std(pathsA(:,kk)),10000)'; %grid for kernel evaluating
%  end
   MAX=max(pathsA,[],1); MIN=min(pathsA,[],1); scaling=MAX-MIN;
  % create meshgrid in 3-dimensions
  [X1,X2,X3,X4]=ndgrid(MIN(1):scaling(1)/(ng-1):MAX(1),...
      MIN(2):scaling(2)/(ng-1):MAX(2),MIN(3):scaling(3)/(ng-1):MAX(3),MIN(4):scaling(4)/(ng-1):MAX(4));
  grid=reshape([X1(:), X2(:), X3(:), X4(:)],ng^d,d); % create points for plotting

   
   % Estimate joint pdf on grid
   [pdf]=akde(pathsA,grid);
   % -------------------------------------------------------------------
	% Replace tiny negatives (due to numerical issues) with eps
    if min(min(pdf))<-0.00001
        warning('kde2d returns negative pdf')
    end
    pdf_ind = (pdf<0);
    pdf(pdf_ind) = eps; 
    
    pdf1=reshape(pdf,size(X1));
    
    pdf0=interpn(pdf1,ObsB(:,1),ObsB(:,2),ObsB(:,3),ObsB(:,4));
    % -------------------------------------------------------------------
	% Replace NaN due to the observation being outside the domain of the
	% pdf with eps:
    %if isnan(pdf0)
    %    if obsB(:,1)<min(min(ObsB)) || obsB(:,1)>max(max(X)) | obsY<min(min(Y)) | obsY>max(max(Y)) 
    %        pdf0 = eps;
    %    else
    %        error('pdf interpolation returns NaN')
    %    end
    % end
	% -------------------------------------------------------------------

    jointLogScore(i) = log(pdf0);
    
    
end



end