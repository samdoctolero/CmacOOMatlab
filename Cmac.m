classdef Cmac 
    %CMAC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access=protected)
        mNumLayers
        mNumQ
        mHashtable
        mOffset
        mMemSize
        mMinState
        mMaxState
        mGammas
        mLocations
        mMem
        mWeights
        mNumInputs
        mNumOutputs
        mCurrentOutput
    end
    
    methods
        function obj = Cmac(numLayers, numQ ...
                , memsize, numinputs ...
                , numoutputs, minstate, maxstate)
            %CMAC Construct an instance of this class
            obj.mMem = zeros(memsize,numoutputs);
            
            obj.mHashtable = zeros(1,numinputs*numLayers*numQ + numQ*numLayers + numQ);
            for i=1:numinputs*numLayers*numQ + numQ*numLayers + numQ
                hashnum = rand();
                hashnum = hashnum*(10000);
                obj.mHashtable(i) =hashnum;
            end
            
            obj.mOffset = zeros(numinputs,numLayers);
            for i=1:numinputs
                for j=1:numLayers
                    obj.mOffset(i,j)=rand();
                end
            end
            
            obj.mLocations = ones(numLayers,1);
            obj.mGammas = zeros(1,numLayers);
            obj.mWeights = zeros(numLayers, numoutputs);
            obj.mCurrentOutput = zeros(numoutputs,1);
            
            % set the rest of the properties
            obj.mNumLayers = int32(numLayers);
            obj.mNumQ = int32(numQ);
            obj.mMemSize = int32(memsize);
            obj.mMinState = minstate;
            obj.mMaxState = maxstate;
            obj.mNumInputs = int32(numinputs);
            obj.mNumOutputs = int32(numoutputs);
            
        end
        
        function output = GetOutput(obj, originalinputs)
            %GETOUTPUT Summary of this method goes here
            %   Normalizes the inputs, finds the basis function
            %   values, finds the active weights based on the inputs,
            %   then the basis functions are multiplied by the 
            %   active weights to get the output.
            
            % normalize inputs
            originalinputs = originalinputs';
            denom = obj.mMaxState(1:obj.mNumInputs)-obj.mMinState(1:obj.mNumInputs);
            input= (originalinputs(1:obj.mNumInputs)-obj.mMinState(1:obj.mNumInputs))./denom;
            input=min(input,ones(obj.mNumInputs,1));
            input=max(input,zeros(obj.mNumInputs,1));
            
            % find the basis function values
            gammasum = 0;
            for j=1:obj.mNumLayers
                
                totalLocations= 0.0 ;
                cell = zeros(obj.mNumInputs,1);
                for i=1:obj.mNumInputs
                    obj.mGammas(j)=1.0;
                    place = input(i)*double(obj.mNumQ -1) + obj.mOffset(i,j);
                    cell(i)=floor(place);      %/*activated cell*/
                    h=place-cell(i);  %/*between 0 and 1*/
                    if(h > 1)
                        fprintf('Possible error');
                    end
                    func=16.0*(h*h-2.0*h*h*h+h*h*h*h);%spline
                    
                    
                    obj.mGammas(j)=obj.mGammas(j)*func;
                    var = cell(i)+1 + obj.mNumQ*(j-1) + obj.mNumQ*obj.mNumLayers*(i-1);
                    totalLocations =totalLocations+ double(obj.mHashtable(var));
                    loc = int32(floor( mod(totalLocations,obj.mMemSize-1)) )+1;
                    obj.mLocations(j)=loc;
                    
                    gammasum = gammasum + obj.mGammas(j);
                    
                end
            end
            
            % normalize the basis functions
            obj.mGammas=obj.mGammas./gammasum;
            
            if(obj.mGammas > 1)
                fprintf('Gamma2 is too large');
            end
            
            % find the weights
            for k = 1: obj.mNumOutputs
                obj.mWeights(:,k) = obj.mMem(obj.mLocations,k);
            end
            
            % this can be commented out
            if(obj.mWeights > 5)
                fprintf('Weights possibly too large');
            end
            
            % calculate the output
            for i = 1: obj.mNumOutputs
                obj.mCurrentOutput(i) = obj.mGammas*obj.mWeights(:,i);
            end
            
            output = obj.mCurrentOutput;
        end
        
        function TrainEmod(obj,betadt, nu, z, normZ)
            %TRAINEMOD Summary of this method goes here
            %   Updates the weights based on an
            %   E-modification scheme.
            
            for k = 1 : obj.mNumOutputs
                    a1 = obj.mGammas'*z(k);
                    a2 =nu*normZ*obj.mWeights(:,k);
                    
                    dW = betadt*(a1 - a2);
                    obj.mMem(obj.mLocations,k) = ...
                        obj.mMem(obj.mLocations,k) + dW;
            end
        end
        
    end
end

