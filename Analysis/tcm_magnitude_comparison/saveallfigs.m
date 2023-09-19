function saveallfigs(output_path, tikz_flag, pdf_flag)
% SAVEALLFIGS  Saves and exports all open figures to FIG and PNG files.
%   SAVEALLFIGS saves all the figures in the current folder.
%
%   SAVEALLFIGS(OUTPUT_PATH) saves all the figures into the directory OUTPUT_PATH.
%   OUTPUT_PATH must be an absolute path to an existing folder.
%
%   The function works in LiveScript, if it is placed in a separate section.
%   
%   V1.0 01/12/19 by Samuel Londner


%validate/correct output path
if nargin==0 || isempty(output_path) || output_path==""
    output_path = pwd;
end
output_path = char(output_path); %in case string was given
%TODO: support relative/non-existing folder
%if output_path(end) ~= filesep
%    output_path(end+1) = filesep;
%end
%[filepath,name,ext] = fileparts(output_path) ;

all_figures = findobj('Type', 'figure');


dirname = output_path;
%TODO: if not absolute path add pwd in front of it

for jj=1:length(all_figures)
    current_fig = all_figures(jj);
    try
        % attempt to find the title of the figure
        try
            ax = findall(current_fig, 'type', 'axes');
            if isempty(ax) %empty figure with no axes
                continue;
            end
            if length(ax)>1 %if there are several axes (subplot?) choose first
                ax = ax(1);
            end
            figname = ax.Title.String;
            figname = matlab.lang.makeValidName(figname);
        catch %if failed for some reason (no title etc.) assign empty name
            figname = '';
        end
        
        % save visibility status (needed so that LiveScript saved FIGs are openable)
        old_state = get(current_fig, 'Visible');
        set(current_fig, 'Visible', 'on');
        
        %build fig filename
        filename = fullfile(dirname, ['fig' num2str(jj) '_' figname]);
        % export to PNG
        print(current_fig, [filename '.png'], '-dpng', '-r450');
        % save to FIG
        savefig(current_fig, filename);
        if tikz_flag
            % save to tikz
            matlab2tikz('figurehandle',current_fig,'filename',strcat(filename, '.tikz'), 'showInfo', false)
        end
        % export as pdf
        if nargin == 3
            exportgraphics(current_fig,strcat(filename, '.pdf'),'ContentType','image')
        else
		    exportgraphics(current_fig,strcat(filename, '.pdf'),'ContentType','vector')
        end
        
        %restore visibility status
        set(current_fig, 'Visible', old_state);
    catch Exception
        disp(['Figure ' num2str(jj) ' not saved:' getReport(Exception)]);
    end
end

end