classdef SupervisoryCmac < Cmac
    %SUPERVISORYCMAC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = private)
        mSupervisoryValues
    end
    
    methods
        function obj = SupervisoryCmac(numLayers, numQ ...
                , memsize, numinputs ...
                , numoutputs, minstate, maxstate ...
                , supervisoryValues)
            %SUPERVISORYCMAC Construct an instance of this class
            %   Detailed explanation goes here
            
            % call the base class constructor
            obj@Cmac(numLayers, numQ ...
                , memsize, numinputs ...
                , numoutputs, minstate, maxstate);
            
            % set own member values
            obj.mSupervisoryValues = supervisoryValues;
        end
        
        function output = GetOutput(obj, originalinputs)
            % use the same get output
            output = GetOutput@Cmac(obj, originalinputs);
            % bias weigh output with supervisory values
            output = output + obj.mSupervisoryValues;
        end
        
        function TrainEmod(obj,betadt, nu, z, normZ)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            for k = 1 : obj.mNumOutputs
                a1 = obj.mGammas'*z(k);
                a2 =nu*normZ*(obj.mSupervisoryValues(k) - obj.mWeights(:,k));
                
                dW = betadt*(a1 + a2);
                obj.mMem(obj.mLocations,k) = ...
                    obj.mMem(obj.mLocations,k) + dW;
            end
        end
    end
end

