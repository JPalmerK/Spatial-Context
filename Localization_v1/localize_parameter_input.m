function [localize_struct,array_struct]=localize_parameter_input(hydrophone_struct,array_struct,parm);

  for(j=1:length(hydrophone_struct))
        
        depth(j)=hydrophone_struct(j).depth;
    end

ssp_choice = input('Enter (1) for constant SSP, (2) for user provided file');
localize_struct.parm.ssp_choice=ssp_choice;

if(ssp_choice ==2)
    
    load ssp_PMRF_winter   
    zmax=max(ssp(:,1));    
    if(max(depth) > zmax)
        
        display('SSP does not cover full depth');
        localization_struct=[];
        return
    end
    
else
   ssp_C = input('Input constant sound speed value');  
   ssp(:,1)=linspace(0,max(depth),101);
   ssp(:,2)=ssp_C;

end

localize_struct.parm.grid_depth = input('Nominal depth for TDOA calculations');
localize_struct.parm.num_calls = input('number of calls to cross correlate?');
localize_struct.parm.number_maxima = input('number of maxima?');

% this is outdated feature
%localize_struct.parm.phase = input('enter phase range (not in seconds yet?) ');

localize_struct.parm.lsq_cutoff=input('LSQ cutoff (.005 recommended)');
localize_struct.parm.pow=input('Power of FFT amplitude for cross corr (1.5)');

% currently not used
%localize_struct.parm.cc_cutoff=input('cross corr minimum score (0 through 1)');

localize_struct.parm.ssp = ssp;

localize_struct.parm.c_nominal= input('cst speed to use for max travel times');

[array_struct] = setup_TDOA_grid(hydrophone_struct, array_struct, localize_struct,parm); 
