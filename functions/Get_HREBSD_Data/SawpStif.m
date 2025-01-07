function [DirStiffness,R] = SawpStif(Dirresults,RefID)
%   if ~isfield(Dir.Maps,'Correct_Stiffness_Tensor')
  try
    DirXEBSD = [erase(Dirresults, '.mat') '_XEBSD.mat'];
    load(DirXEBSD,'Map_RefID','C_voight', 'Map_EBSD_MTEX',...
                'MicroscopeData','Data_InputMap','GrainData'); 
  catch err
    try
        load(Dirresults,'Map_RefID','C_voight', 'Map_EBSD_MTEX',...
                'MicroscopeData','Data_InputMap','GrainData');
    catch
    [~,A,~] = fileparts(fileparts(Dirresults));
    fname = fullfile(fileparts(Dirresults), [A '_XEBSD.mat']);
    load(fname,'Map_RefID','C_voight', 'Map_EBSD_MTEX',...
                'MicroscopeData','Data_InputMap','GrainData');
    end
  end
  
if ~exist('RefID','var')
    set(0,'DefaultLineMarkerSize',10)
    % plot to select the boundray
    close all;  	warning off;        s3=subplot(1,1,1);
    imagesc(Data_InputMap.X_axis,Data_InputMap.Y_axis,Map_RefID); hold on
    axis off;           axis image;         axis xy;        s3.YDir='reverse';   
    colormap jet;       s3.XDir='reverse';        colorbar          
    title('Respond in the Command Line')
    if isempty(GrainData.RefPoint)
        [GrainData.RefPoint] = to_label(Data_InputMap, Map_RefID);
    else
        for i=1:length(GrainData.RefPoint.x)       
            GrainData.RefPoint.prop.labels{i} = num2str(i);     
        end
    end
    scatter(GrainData.RefPoint.prop.x,GrainData.RefPoint.prop.y,'k','filled');
    scatter(GrainData.RefPoint.prop.x,GrainData.RefPoint.prop.y,'w');
    labelpoints(GrainData.RefPoint.prop.x,GrainData.RefPoint.prop.y,GrainData.RefPoint.prop.labels);     
    hold off;      set(gcf,'position',[30 50 1300 950])
    if length(GrainData.RefPoint.x) ~=1
        Spec    = input('Stifness tensor of which grain?   '); % select
    else
        Spec = 1;
    end
    close;
    R = Map_EBSD_MTEX(sub2ind([MicroscopeData.NROWS,MicroscopeData.NCOLS],...
                   GrainData.RefPoint.prop.yi(Spec),GrainData.RefPoint.prop.xi(Spec)))...
                  .orientations.matrix;
    DirStiffness = Korsunsky_StiffnessRot(R,C_voight(:,:,Spec)); 
else
    R = Map_EBSD_MTEX(sub2ind([MicroscopeData.NROWS,MicroscopeData.NCOLS],...
                   GrainData.RefPoint.prop.yi(RefID),GrainData.RefPoint.prop.xi(RefID)))...
                  .orientations.matrix;
    DirStiffness = Korsunsky_StiffnessRot(R,C_voight(:,:,RefID)); 
end
end