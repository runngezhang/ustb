function [iqGainedA, iqGainedB, gcurveA, gcurveB] = gainCompensation(iqA, iqB, method,varargin)
    % carries out gain compensation of two images
    % [iqGainedA, iqGainedB, gcurveA, gcurveB] = gainCompensation(iqA, iqB, method,varargin)
    % both the gained IQ signal and the gain curves are returned. 
    % pseudo code: iqGainedA = iqA.*gcurveA; iqGainedB = iqB.*gcurveB;
    %    gcourve(A|B) will be a vector if it is constant for all lines)
    switch method
        case 'latavg'
            % varargin: beams, start_depth, poly_order
            % takes the lateral average over all of the beams or the beams 
            % given in parameter beams of the image A and B.
            % values up to start_depth are ignored.
            % Then computes the gain relation: A/B, approximates it with a 
            % polynome of order poly_order. 
            % This function is then used to make B equally gained as A:
            % iqGainedB = iqB.*(poly(A/B)) (pseudocode)
            % iqGainedA = iqA
            beamdim = 2;
            if nargin < 4
                beams = 1:size(iqA,beamdim);
            else
                beams = varargin{1};
            end
            if nargin < 5
                start_depth = 1; %samples                
            else
                start_depth = varargin{2};
            end
            if nargin < 6
                poly_order = 3;
            else
                poly_order = varargin{3};
            end;
            computeSum = @(x)sum(abs(x(:,beams,1)), beamdim);
            iqAsum = computeSum(iqA);
            iqBsum = computeSum(iqB);
            
            gain_curve = iqAsum./iqBsum;
            gain_curve(isinf(gain_curve)) = 0;
            % drop values at the start
            gain_curve(1:start_depth-1) = gain_curve(start_depth);
            
            % approximate with a polynome
            p = polyfit(1:length(gain_curve),gain_curve',poly_order);
            gcurveB = polyval(p,1:length(gain_curve));
            
            % carry out gaining
            iqGainedA = iqA;
            iqGainedB = iqB.*repmat(gcurveB.', [1, size(iqB,2) size(iqB,3)]);
            gcurveA = ones(size(gcurveB));
        case 'apriori'
            gain_curve = load_gain_curve(varargin{1},SURF_TCI,varargin{2},varargin{3},size(iqA,1));
            gain_curve = gain_curve(:,ones(1,size(iqA,2)),ones(1,size(iqA,3)));
            iqGainedA = iqA;
            iqGainedB = gain_curve.*iqB;
            
            gcurveA = ones(size(gain_curve(:,1)));
            gcurveB = gain_curve(:,1);
        case 'parametric'
            % varargin: 
            %   probe id
            %   application (SURF_TCI)
            %   c
            %   fs
            %   N number of samples
            % will get a gaincurve and then multiply the gaincurve with iqB
            probeid = varargin{1};
            application = varargin{2};
            c = varargin{3};
            fs = varargin{4};
            N = size(iqA,1);
            gain_curve = load_gain_curve(varargin{1},application,c,fs,N);
            gain_curve = repmat(gain_curve,[1,size(iqB,2),size(iqB,3)]);
            iqGainedA = iqA;
            iqGainedB = gain_curve.*iqB;
            gcurveA = ones(size(gain_curve(:,1)));
            gcurveB = gain_curve(:,1);
        case 'parametric-filt-plus'
            % Hack Svein-Erik 
            % varargin: 
            %   probe id
            %   application (SURF_TCI)
            %   c
            %   fs
            %   N number of samples
            % will get a gaincurve and then multiply the gaincurve with iqB
            probeid = varargin{1};
            application = varargin{2};
            c = varargin{3};
            fs = varargin{4};
            N = size(iqA,1);
            gain_curve = load_gain_curve(varargin{1},application,c,fs,N);
            gain_curve_raw = gain_curve;
            thrs = 3.8;
            gain_curve_raw(gain_curve_raw>thrs) = thrs;
            smoothx = 100;
            gain_curve = filtfilt(ones(1,smoothx)/smoothx,1,gain_curve_raw);
            gain_curve = gain_curve*10^(3.5/20);
            gain_curve = repmat(gain_curve,[1,size(iqB,2),size(iqB,3)]);
            iqGainedA = iqA;
            iqGainedB = gain_curve.*iqB;
            gcurveA = ones(size(gain_curve(:,1)));
            gcurveB = gain_curve(:,1);
        case 'parametric-filt-half'
            % Hack Svein-Erik 
            % varargin: 
            %   probe id
            %   application (SURF_TCI)
            %   c
            %   fs
            %   N number of samples
            % will get a gaincurve and then multiply the gaincurve with iqB
            probeid = varargin{1};
            application = varargin{2};
            c = varargin{3};
            fs = varargin{4};
            N = size(iqA,1);
            gain_curve = load_gain_curve(varargin{1},application,c,fs,N);
            gain_curve_raw = gain_curve;
            thrs = 3.8;
            gain_curve_raw(gain_curve_raw>thrs) = thrs;
            smoothx = 100;
            gain_curve = filtfilt(ones(1,smoothx)/smoothx,1,gain_curve_raw);
%             gain_curve = gain_curve*10^(3.5/20);
            gain_curve = gain_curve*10^(-3.5/20);
            gain_curve = repmat(gain_curve,[1,size(iqB,2),size(iqB,3)]);
            iqGainedA = iqA;
            iqGainedB = gain_curve.*iqB;
            gcurveA = ones(size(gain_curve(:,1)));
            gcurveB = gain_curve(:,1);
        otherwise
            error('Method %s unknown',method);
    end
end
