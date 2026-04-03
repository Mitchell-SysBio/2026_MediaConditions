function [reClassWell, reClass] = inputDialog
%MVDLG moveable input dialog box.
%  ANSWER = MVDLG(PROMPT,TITLE,POSITION) creates a dialog box that returns
%  user input for a prompt in cell array ANSWER. PROMPT is a string. TITLE
%  is a string that species a title for the dialog box. POSITION is a
%  four-element vector that specifies the size and location on the screen
%  of the dialog window. Specify POSITION in normalized coordinates in the
%  form: [left, bottom, width, height] where left and bottom are the
%  distance from the lower-left corner of the screen to the lower-left
%  corner of the dialog box. width and height are the dimensions.  
%
%
%  Examples:
%
%  prompt='Enter a value here:';
%  name='Data Entry Dialog';
%  position = [.5 .5 .2 .1];
%
%  answer=mvdlg(prompt,name,position);
%
%
%  prompt='Enter a value here:';
%  name='Data Entry Dialog';
%  position = [.7 .2 .2 .3]; 
%  answer=mvdlg(prompt,name,position);

% Set up main window figure
position = [.8 .5 .3 .4];

dialogWind = figure('Units','normalized','Position',position,'MenuBar','none','NumberTitle','off');

% Set up GUI controls and text
prompt1 = 'grey = junk; red = empty/ tiny pellet; blue = marginal; green = full';
prompt = 'Enter wellname in space seperated list';
uicontrol('style','text','String',sprintf('%s\n%s',prompt1, prompt),'Units','normalized', ...
    'Position',[.05,.9,0.9,.1],'HorizontalAlignment','left')
hedit1 = uicontrol('style','edit','Units','normalized','Position', ...
    [.05,.8,1,.1],'HorizontalAlignment','left');
prompt2 = '01_messy, 02_bubbles, 03_empty, 04_pellet_tiny, 05_cloud_small, 06_cloud_white, 07_pellet_cloudy_small, 08_pellet_small, 09_cloud_large, 10_cloud_and_pellet_large, 11_pellet_cloudy_large, 12_pellet_large';
prompt3 = 'Enter new class as number in space seperated list';
uicontrol('style','text','String',sprintf('%s\n%s',prompt2, prompt3),'Units','normalized', ...
    'Position',[.05,.6,0.9,.1],'HorizontalAlignment','left')
hedit2 = uicontrol('style','edit','Units','normalized','Position', ...
    [.05,.5,1,.1],'HorizontalAlignment','left');
okayButton = uicontrol('style','pushbutton','Units','normalized',...
    'position', [.05,.1,.8,.3],'string','OK','callback',@okCallback);

% Initialize ANSWER to empty cell array
reClassWell = {};
reClass = {};

% Wait for user input, and close once a button is pressed 
uiwait(dialogWind);

% Callbacks for 'OK' and 'Cancel' buttons  
function okCallback(hObject,eventdata)
        reClassWell = cell(get(hedit1,{'String'}));
        reClass = cell(get(hedit2,{'String'}));
        uiresume(dialogWind);
        close(dialogWind);
end
end
