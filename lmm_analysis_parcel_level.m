function [p, r2_adj, BIC, r] = lmm_analysis_parcel_level(ConMats, SubjRecs, RecAges, LMM_order)
%LMM_ANALYSIS_PARCEL_LEVEL Parcel-level statistical test with linear mixed effects model (LMM)

% INPUT ARGUMENTS
%   ConMats: cell array {1 x N freq.} of connectivity matrices for all 
%   recordings [N parcels x N parcels x N subj.]. Connectivity matrices must 
%   be symmetric square matrices and ordered so that all recordings of a subject are
%   consecutively.
%   SubjRecs: double vector [N recs x 1] which specifies recordings of each
%   subj as ascending integers from 1 to N subj with the same integer for all
%   recordings of the same subject. Order of recordings should be consistent
%   with ConMats.
%   RecAges: double vector [N recs x 1] with subject ages at each recording.
%   Order of recordings should be consistent with ConMats and SubjRecs.
%   LMM_order (optional): integer value 1 or 2, where 1 = linear (first order) and 2 =
%   quadratic (second order) LMM. Default is 1.

% OUTPUT ARGUMENTS
%   p: cell array {1 x N freq.} of p-value matrices [1 x N parcels]
%   r2_adj: cell array {1 x N freq.} of adjusted R squared  matrices [1 x N parcels]
%   BIC: cell array {1 x N freq.} of BIC value matrices [1 x N parcels]
%   r: cell array {1 x N freq.} of r-value matrices [1 x N parcels]

% Get parameters from Input arguments
N_parcels = size(ConMats{1,1},1);    % Number of parcels
N_Fc = size(ConMats,2);              % Number of frequency bands

% Check arguments
if ~exist('LMM_order', 'var')
    LMM_order = 1;
end

% Initialize variables
p{1,N_Fc} = [];
r{1,N_Fc} = [];
r2_adj{1,N_Fc} = [];
BIC{1,N_Fc} = [];
K = zeros(2,N_Fc);
   
    for f = 1:N_Fc
        
        % Get connectivity matrices
        connectivity_matrices = ConMats{1,f};

        % Parcel level
        for p1 = 1:N_parcels            

            % Get connectivity vector
            Connectivity = squeeze(mean(connectivity_matrices(p1,:,:), 'omitnan')); 

            % Prepare table for LMME
            Table = table(SubjRecs, RecAges, Connectivity);

            % Perform first or second order LMM
            if LMM_order == 1
                % First order LMM
                lmm_buf = fitlme(Table, 'Connectivity ~ 1 + RecAges + (1 + RecAges | SubjRecs)');
                p_buf = lmm_buf.Coefficients.pValue(2);  
                r_buf = lmm_buf.Coefficients.Estimate(2);
                r2_adj_buf = lmm_buf.Rsquared.Adjusted;  
                BIC_buf = lmm_buf.ModelCriterion.BIC;    
            elseif LMM_order == 2
                % Second order LMM
                lmm_buf = fitlme(Table, 'Connectivity ~ 1 + RecAges + RecAges^2 + (1 + RecAges | SubjRecs)');
                p_buf = lmm_buf.Coefficients.pValue(3);   
                r2_adj_buf = lmm_buf.Rsquared.Adjusted;
                BIC_buf = lmm_buf.ModelCriterion.BIC;
            else
                disp(['Unaccepted value for LMM_order.'])
                break
            end

            p{1,f}(p1) = p_buf;
            r2_adj{1,f}(p1) = r2_adj_buf;             
            BIC{1,f}(p1) = BIC_buf;

            if LMM_order == 1                    
                r{1,f}(p1) = r_buf;
            end

            clear Table
            clear Connectivity
            clear lmm_buf p_buf r2_adj_buf BIC_buf

        end    

    end
    
end


