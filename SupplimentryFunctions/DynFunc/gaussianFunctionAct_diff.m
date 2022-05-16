function dfact_dlce_nor = gaussianFunctionAct_diff(lce_nor, b1, b2, b3, b4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the differential function of active muscle force-length relationship 
% gaussian function, from Groote 2016 paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fact = b1*exp(-0.5.*(lce_nor - b2).^2./(b3 + b4.*lce_nor));

dfact_dlce_nor = fact.*(-(lce_nor - b2).*(b3 + b4.*lce_nor) +...
                        b4*0.5.*(lce_nor - b2).^2)./((b3 + b4.*lce_nor).^2);

end