function Plot_update_burn(folder_name,output_name,new_burnin)
% Plot_update_burn(folder_name,output_name,new_burnin)
% note: burn-in in percent of total iterations
% will over-write the figures under "figures"

load(['models/' folder_name '/' output_name '.mat']);

burn = fix(new_burnin*length(ensemble));

PlotEnsembleRef2D_verD
