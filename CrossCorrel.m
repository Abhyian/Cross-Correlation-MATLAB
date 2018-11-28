function varargout = CrossCorrel(varargin)
% CROSSCORREL MATLAB code for CrossCorrel.fig
%      CROSSCORREL, by itself, creates a new CROSSCORREL or raises the existing
%      singleton*.
%
%      H = CROSSCORREL returns the handle to a new CROSSCORREL or the handle to
%      the existing singleton*.
%
%      CROSSCORREL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CROSSCORREL.M with the given input arguments.
%
%      CROSSCORREL('Property','Value',...) creates a new CROSSCORREL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CrossCorrel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CrossCorrel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CrossCorrel

% Last Modified by GUIDE v2.5 08-Jun-2017 14:16:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CrossCorrel_OpeningFcn, ...
                   'gui_OutputFcn',  @CrossCorrel_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CrossCorrel is made visible.
function CrossCorrel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CrossCorrel (see VARARGIN)

% Choose default command line output for CrossCorrel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CrossCorrel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CrossCorrel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in channel1.
function channel1_Callback(hObject, eventdata, handles)
% hObject    handle to channel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname1,dirname] = uigetfile('*.mat','Select video')
h1 = findobj('Tag','slider1');
h2 = findobj('Tag','slider1');
channel1 = load(fname1,'-mat');
size(channel1.ImageWavelet,3)
set(hObject,'userdata',channel1.ImageWavelet);

 set(h1,'Min',1/size(channel1.ImageWavelet,3));
% 
 set(h1,'Max',1);
set(h1,'value',1/size(channel1.ImageWavelet,3));
set(h1,'sliderstep',[1/size(channel1.ImageWavelet,3),3/size(channel1.ImageWavelet,3)]);

set(h2,'Min',1/size(channel1.ImageWavelet,3));
% 
 set(h2,'Max',1);
set(h2,'value',1/size(channel1.ImageWavelet,3));
set(h2,'sliderstep',[1/size(channel1.ImageWavelet,3),3/size(channel1.ImageWavelet,3)]);




% --- Executes on button press in channel2.
function channel2_Callback(hObject, eventdata, handles)
% hObject    handle to channel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname1,dirname] = uigetfile('*.mat','Select video')
h1 = findobj('Tag','slider1');
channel2 = load(fname1,'-mat');
size(channel2.ImageWavelet,3)
set(hObject,'userdata',channel2.ImageWavelet);
set(h1,'Min',1/size(channel2.ImageWavelet,3));
set(h1,'Max',1);
set(h1,'value',1/size(channel2.ImageWavelet,3));
set(h1,'sliderstep',[1/size(channel2.ImageWavelet,3),3/size(channel2.ImageWavelet,3)]);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h1 = findobj('Tag','channel1');
h2 = findobj('Tag','channel2');
channel1 = get(h1,'UserData');
channel2 = get(h2,'UserData');
k = round(get(hObject,'value')*size(channel1,3))
axes(handles.axes1);
imagesc(channel1(:,:,k));
hold(handles.axes1,'off')
axes(handles.axes2);
imagesc(channel2(:,:,k));
hold(handles.axes2,'off')




% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in crop.
function crop_Callback(hObject, eventdata, handles)
% hObject    handle to crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
hold(handles.axes6,'off')
rect = getrect(handles.axes1);
set(hObject,'UserData',rect);
reg_edge = {[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]};
I1 = getimage(handles.axes1);
Ic1 = imcrop(I1,rect);
axes(handles.axes3);
imagesc(Ic1)
I2 = getimage(handles.axes2);
Ic2 = imcrop(I2,rect);
axes(handles.axes4);
imagesc(Ic2)
axes(handles.axes5);
Icorr= xcorr2((Ic1),(Ic2));
imagesc(Icorr);
hold on
plot(size(Icorr,2)/2,size(Icorr,1)/2,'+r')
hold off


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

h1 = findobj('Tag','channel1');
h2 = findobj('Tag','channel2');
h3 = findobj('Tag','crop');
rect = get(h3,'userdata');
channel1 = get(h1,'UserData');
channel2 = get(h2,'UserData');
k = round(get(hObject,'value')*size(channel1,3))
axes(handles.axes1);
imagesc(channel1(:,:,k));
hold(handles.axes1,'off')
axes(handles.axes2);
imagesc(channel2(:,:,k));
hold(handles.axes2,'off')

axes(handles.axes1);


reg_edge = {[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]};
I1 = getimage(handles.axes1);
Ic1 = imcrop(I1,rect);
axes(handles.axes3);
imagesc(Ic1)
I2 = getimage(handles.axes2);
Ic2 = imcrop(I2,rect);
axes(handles.axes4);
imagesc(Ic2)
axes(handles.axes5);
Icorr= xcorr2((Ic1),(Ic2));
imagesc(Icorr);
hold on
plot(size(Icorr,2)/2,size(Icorr,1)/2,'+r')
size(Icorr)
plot(handles.axes6,k,Icorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2))/max(Icorr(:)),'*');
hold(handles.axes6,'on')
round(size(Icorr,1)/2)
round(size(Icorr,2)/2)
Icorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2))/max(Icorr(:))
hold off

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




h1 = findobj('Tag','channel1');
h2 = findobj('Tag','channel2');
h3 = findobj('Tag','crop');

h4 = findobj('Tag','tl');
h5 = findobj('Tag','th');
h6 = findobj('Tag','ROI');

tmin = str2double(get(h4,'string'));
tmax = str2double(get(h5,'string'));
rect = get(h3,'userdata');
channel1 = get(h1,'UserData');
% save('ch1_cropped.mat','channel1');

max_ch1 = max(channel1(:))
thr_ch1 = 0;
channel2 = get(h2,'UserData');
% save('ch2_cropped.mat','channel1');
max_ch2 = max(channel2(:));
thr_ch2 = 0;

% Initialization
filename = '1_result_testdata.xlsx';
ch1_binary = '2_result_ch1_binary.tif';
ch2_binary = '2_result_ch2_binary.tif';
ch1_or_ch2file = '2_result_ch1orch2.tif';
ch1_and_ch2file = '2_result_ch1andch2.tif';
ch1_correl_ch2file = '2_result_correlation2D.tif';
statistics = '1_result_statisticalAnalysis.xlsx';
data = zeros(size(channel1,3),3);
% Imch1 = zeros(size(channel1(:,:,1)));
% Imch2 = zeros(size(channel1(:,:,1)));
num_channel = zeros(size(channel1,3)+1,5);

num_channel(1,1:11) = zeros(1,11);
Ic1 = zeros(size(channel1(:,:,1)));
Ic2 = zeros(size(channel2(:,:,1)));



IcT1 = imcrop(Ic1,rect);
IcT2 = imcrop(Ic2,rect);
chan1 = zeros(size(IcT1,1),size(IcT1,2),size(channel1,3));
chan2 = zeros(size(IcT2,1),size(IcT2,2),size(channel2,3));


xc_size = size(xcorr2(IcT1,IcT2));
crossCorr = zeros(xc_size(1),xc_size(2),size(channel1,3));
save 1_result_crosscorr.mat crossCorr -v7.3;
save 1_result_chan1_cropped.mat chan1 -v7.3;
save 1_result_chan2_cropped.mat chan2 -v7.3;
% k = round(get(hObject,'value')*size(channel1,3))
hold(handles.axes6,'off');
for k = 1:size(channel1,3)
    
%     Imch1 = zeros(size(channel1(:,:,1)));
%     Imch2 = zeros(size(channel1(:,:,1)));

    axes(handles.axes1);
    imagesc(channel1(:,:,k));

   
    
    
    hold(handles.axes1,'off')
    axes(handles.axes2);
    imagesc(channel2(:,:,k));
    hold(handles.axes2,'off')
    axes(handles.axes1);
    % reg_edge = {[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]};
    I1 = getimage(handles.axes1);
    Ic1 = imcrop(I1,rect);
    axes(handles.axes3);
    imagesc(Ic1)
    I2 = getimage(handles.axes2);
    Ic2 = imcrop(I2,rect);
    
    Imch1 = zeros(size(Ic1(:,:)));
    Imch1(Ic1(:,:)>thr_ch1) = 1  %%binary image threshold channel 1
    
    Imch2 = zeros(size(Ic2(:,:)));
    Imch2(Ic2(:,:)>thr_ch2) = 1; 
    
    axes(handles.axes4);
    imagesc(Ic2)
    axes(handles.axes5);
    Icorr= xcorr2((Ic1),(Ic2));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    normcorr= xcorr2((Ic1),(Ic2));
    imagesc(Icorr);
    example = matfile('1_result_crosscorr.mat','Writable',true);
    example.crossCorr(:,:,k) = normcorr;
    
    example_ch1 = matfile('1_result_chan1_cropped.mat','Writable',true);
    example_ch1.chan1(:,:,k) = Ic1;
    
    example_ch2 = matfile('1_result_chan2_cropped.mat','Writable',true);
    example_ch2.chan2(:,:,k) = Ic2;
    imwrite(Icorr, ch1_correl_ch2file, 'writemode', 'append');
    hold on
    plot(size(Icorr,2)/2,size(Icorr,1)/2,'+r')
    size(Icorr)
    plot(handles.axes6,k,normcorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2)),'*r');
    hold(handles.axes6,'on')
    round(size(Icorr,1)/2)
    round(size(Icorr,2)/2)
    Icorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2));
    data(k,1)=k;
    data(k,2)=Icorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2));
    data(k,3)=normcorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2));
    pause(0.001)
    hold off
    
    
    IcT1 = IcT1+Imch1;
    imwrite(Imch1, ch1_binary, 'writemode', 'append' );
     %%binary image threshold channel 2
    IcT2 = IcT2+Imch2;
    
    imwrite(Imch2, ch2_binary, 'writemode', 'append' );
    
    ch1_or_ch2 = (Imch1+Imch2);
    ch1_or_ch2(ch1_or_ch2>1)=1;
    imwrite(ch1_or_ch2, ch1_or_ch2file, 'writemode', 'append');
    ch1_and_ch2 = Imch1.*Imch2;
    imwrite(ch1_and_ch2, ch1_and_ch2file, 'writemode', 'append');
    
    %Statistical analysis
    CC_ch1 = bwconncomp(Imch1);
    CC_ch2 = bwconncomp(Imch2);
    CC_ch1_or_ch2 = bwconncomp(ch1_or_ch2);
    CC_ch1_and_ch2 = bwconncomp(ch1_and_ch2);
    
   
    num_channel(k+1,1) = CC_ch1.NumObjects;
    num_channel(k+1,2) = CC_ch2.NumObjects;
    num_channel(k+1,3) = CC_ch1_or_ch2.NumObjects;
    num_channel(k+1,4) = CC_ch1_and_ch2.NumObjects;
    num_channel(k+1,5) = num_channel(k+1,4)/num_channel(k+1,3)*100;
    num_channel(k+1,7) = length(find(Imch1==1));
    num_channel(k+1,8) = length(find(Imch2==1));
    num_channel(k+1,9) = length(find(ch1_or_ch2==1));
    num_channel(k+1,10) = length(find(ch1_and_ch2==1));
    num_channel(k+1,11) = num_channel(k+1,10)/num_channel(k+1,9);
     
end
% axes(handles.axes1);
% imagesc(Imch1);
% hold(handles.axes1,'off')
% axes(handles.axes2);
% imagesc(Imch2);
% hold(handles.axes2,'off')
xlswrite(filename,data)
xlswrite(statistics,num_channel);

IcTMax = max(IcT1,IcT2);
maxI = max(IcTMax(:));

axes(handles.axes7);
mesh(IcTMax,'FaceAlpha',0.1)

IcTFilter = IcTMax;
IcTFilter(IcTFilter<tmin) = 0;
% IcTFilter(IcTFilter>tmax) = 0;
IROI = maxI*imregionalmax(IcTFilter,26);
% IcTFilter = imregionalmax(IcTFilter);
hold(handles.axes7,'on')
imagesc(IROI);
IcMxlevel = tmax*ones(size(IcTFilter));
surf(IcMxlevel,'FaceAlpha',0.1)
hold(handles.axes7,'on')
IcMnlevel = tmin*ones(size(IcTFilter));
surf(IcMnlevel,'FaceAlpha',0.1)
hold(handles.axes7,'off')
set(h6,'UserData',IROI);


col_header={'cha1','cha2','OR','AND','Percentage','','Areacha1','Areacha2','AreachaORch2','Areach1ANDch2','PercentageArea'};
xlswrite(statistics,col_header,'Sheet1','A1');     %Write column header
% xlswrite(statistics,cha2,'Sheet1','A1');     %Write column header
% xlswrite(statistics,OR,'Sheet1','A1');     %Write column header
% xlswrite(statistics,AND,'Sheet1','A1');     %Write column header
% xlswrite(statistics,Percentage,'Sheet1','A1');     %Write column header

% --- Executes on button press in ROI.
function ROI_Callback(hObject, eventdata, handles)
% hObject    handle to ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h1 = findobj('Tag','channel1');
h1o = findobj('Tag','chan1');
h2 = findobj('Tag','channel2');
h2o = findobj('Tag','chan2');
h3 = findobj('Tag','crop');
rect = get(h3,'userdata');
channel1 = get(h1,'UserData');
channelori1 =   get(h1o,'UserData');

max_ch1 =max(channel1(:));
thr_ch1 = 0;
channel2 = get(h2,'UserData');
channelori2 =   get(h2o,'UserData');

max_ch2 = max(channel2(:));
thr_ch2 = 0;


Imch1 = zeros(size(channel1(:,:,1)));

Imch2 = zeros(size(channel1(:,:,2)));


% k = round(get(hObject,'value')*size(channel1,3))
hold(handles.axes6,'off');
for k = 1:size(channel1,3)
%     Imch1 = zeros(size(channel1(:,:,1)));
% 
%     Imch2 = zeros(size(channel1(:,:,2)));
%     axes(handles.axes1);
%     imagesc(channel1(:,:,k));
    Imch1(channel1(:,:,k)>thr_ch1) = 1;
%     Imch1(channel1(:,:,k)==0) = 0;
    Imcp1 = imcrop(Imch1,rect);
    Ich1 = imcrop(channel1(:,:,1),rect);
    if (k==1)
            CC = bwconncomp(Ich1)
            sum = 0;
            obj = CC.NumObjects
            for i = 1:obj
                
                    sum = sum + length(CC.PixelIdxList{i});
            end
                
            size_mean = sum/CC.NumObjects;
%             CC = bwconncomp(im_mg)
%             CC.PixelIdxList{2} 
    end
    
%     hold(handles.axes1,'off')
%     axes(handles.axes2);
%     imagesc(channel2(:,:,k));
    Imch2(channel2(:,:,k)>thr_ch2) = 1; 
%     Imch2(channel2(:,:,k)==0) = 0;
    Imcp2 = imcrop(Imch2,rect);
    Ich2 = imcrop(channel2(:,:,1),rect);
%     hold(handles.axes2,'off')
%     axes(handles.axes1);
    % reg_edge = {[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]};
%     I1 = getimage(handles.axes1);
%     Ic1 = imcrop(I1,rect);
%     axes(handles.axes3);
%     imagesc(Ic1)
%     I2 = getimage(handles.axes2);
%     Ic2 = imcrop(I2,rect);
%     axes(handles.axes4);
%     imagesc(Ic2)
%     axes(handles.axes5);
%     Icorr= xcorr2((Ic1),(Ic2));
%     imagesc(Icorr);
%     hold on
%     plot(size(Icorr,2)/2,size(Icorr,1)/2,'+r')
%     size(Icorr)
%     plot(handles.axes6,k,Icorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2))/max(Icorr(:)),'*r');
%     hold(handles.axes6,'on')
%     round(size(Icorr,1)/2)
%     round(size(Icorr,2)/2)
%     Icorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2))/max(Icorr(:))
%     pause(0.001)
%     hold off
    
end
axes(handles.axes3);
imagesc(Imcp1);
hold(handles.axes3,'off')
axes(handles.axes4);
imagesc(Imcp2);
hold(handles.axes4,'off')
% im_mg = or(Imcp1,Imcp2);
im_mg = get(hObject,'UserData');
im_mg(im_mg>0) = 1;
axes(handles.axes5);
imagesc(im_mg);

CC = bwconncomp(im_mg);
% CC.PixelIdxList{2}; 
% ind = CC.PixelIdxList{2} 
% mean(Ich2(ind))
S = regionprops(CC,'BoundingBox' )
map = colormap(jet(CC.NumObjects));
size(map)
for k = 1 : length(S)
  thisBB = 1*S(k).BoundingBox;
  rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],'EdgeColor',map(k,:),'LineWidth',2 )
end


obj1 = zeros(CC.NumObjects,size(channel1,3));
obj2 = zeros(CC.NumObjects,size(channel1,3));
for i = 1:CC.NumObjects
    length(CC.PixelIdxList{i});
%     if (length(CC.PixelIdxList{i})<=size_mean*40)
%         length(CC.PixelIdxList{i})
           i
        for k = 1:size(channel1,3)
            Ich1 = imcrop(channel1(:,:,k),rect);
            Ich2 = imcrop(channel2(:,:,k),rect);

            ind = CC.PixelIdxList{i};

            obj1(i,k) = mean(Ich1(ind));
            obj2(i,k) = mean(Ich2(ind));

%             obj1(i,k) = mean(Ich1(ind));
%             obj2(i,k) = mean(Ich2(ind));


        end
%     else
%         obj1(i,:) = zeros(size(channel1,3),1);
%         obj2(i,:) = zeros(size(channel1,3),1);
%    end
    xcorrobj(i,:) = xcorr(obj1(i,:),obj2(i,:),'unbiased');
%     stdobj1 = std(obj1(i,:));
%     stdobj2 = std(obj2(i,:));
    %xcorrobj(i,:) = xcorrobj(i,:)/(stdobj1*stdobj2);
end

axes(handles.axes6);
zero_lag = (size(xcorrobj,2)+1)/2
colormap(jet(256));
caxis([max(xcorrobj(:))/3 max(xcorrobj(:))])
[min(xcorrobj(:)) max(xcorrobj(:))]
lag = str2num(handles.maxlag.String);
set(handles.axes6,'xlim',[zero_lag-lag,zero_lag+lag]);
save TempoCrossCorr.mat xcorrobj  -v7.3;
Cross_p50m50 = xcorrobj(:,zero_lag-lag:zero_lag+lag);
save Cross_p50m50.mat Cross_p50m50  -v7.3;
%imagesc(Cross_p50m50);
y_cursor_cor = (1:1:CC.NumObjects);
x_cursor_cor = (size(Cross_p50m50,2)+1)/2*ones(1,CC.NumObjects);

%scatter(x_cursor_cor,y_cursor_cor,'filled','d');


imagesc(Cross_p50m50);
hold(handles.axes6,'on')
plot(x_cursor_cor,y_cursor_cor,'+g');
hold(handles.axes6,'off')




% --- Executes on button press in chan1.
function chan1_Callback(hObject, eventdata, handles)
% hObject    handle to chan1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = initialize_read(0,'Select Stacked frames of cell body');
set(hObject,'userdata',data);

h1 = findobj('Tag','slider1');
h2 = findobj('Tag','slider3');




 set(h1,'Min',1/size(data,3));
% 
 set(h1,'Max',1);
set(h1,'value',1/size(data,3));
set(h1,'sliderstep',[1/size(data,3),3/size(data,3)]);

set(h2,'Min',1/size(data,3));
% 
 set(h2,'Max',1);
set(h2,'value',1/size(data,3));
set(h2,'sliderstep',[1/size(data,3),3/size(data,3)]);

 



% --- Executes on button press in chan2.
function chan2_Callback(hObject, eventdata, handles)
% hObject    handle to chan2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = initialize_read(0,'Select Stacked frames of cell body');
set(hObject,'userdata',data);


% --- Executes during object creation, after setting all properties.
function ROI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function tl_Callback(hObject, eventdata, handles)
% hObject    handle to tl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tl as text
%        str2double(get(hObject,'String')) returns contents of tl as a double


% --- Executes during object creation, after setting all properties.
function tl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function th_Callback(hObject, eventdata, handles)
% hObject    handle to th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of th as text
%        str2double(get(hObject,'String')) returns contents of th as a double


% --- Executes during object creation, after setting all properties.
function th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




h1 = findobj('Tag','chan1');
h2 = findobj('Tag','chan2');
h3 = findobj('Tag','crop');

h4 = findobj('Tag','tl');
h5 = findobj('Tag','th');
h6 = findobj('Tag','ROI');

tmin = str2double(get(h4,'string'));
tmax = str2double(get(h5,'string'));
rect = get(h3,'userdata');
channel1 = get(h1,'UserData');
% save('ch1_cropped_tiff.mat','channel1');

max_ch1 = max(channel1(:))
thr_ch1 = 0;
channel2 = get(h2,'UserData');
% save('ch2_cropped_tiff.mat','channel1');
max_ch2 = max(channel2(:));
thr_ch2 = 0;

% Initialization
filename = '1_result_testdata_tiff.xlsx';
ch1_binary = '2_result_ch1_binary_tiff.tif';
ch2_binary = '2_result_ch2_binary_tiff.tif';
ch1_or_ch2file = '2_result_ch1orch2_tiff.tif';
ch1_and_ch2file = '2_result_ch1andch2_tiff.tif';
ch1_correl_ch2file = '2_result_correlation2D_tiff.tif';
statistics = '1_result_statisticalAnalysis_tiff.xlsx';
data = zeros(size(channel1,3),3);
% Imch1 = zeros(size(channel1(:,:,1)));
% Imch2 = zeros(size(channel1(:,:,1)));
num_channel = zeros(size(channel1,3)+1,5);

num_channel(1,1:11) = zeros(1,11);
Ic1 = zeros(size(channel1(:,:,1)));
Ic2 = zeros(size(channel2(:,:,1)));



IcT1 = imcrop(Ic1,rect);
IcT2 = imcrop(Ic2,rect);
chan1 = zeros(size(IcT1,1),size(IcT1,2),size(channel1,3));
chan2 = zeros(size(IcT2,1),size(IcT2,2),size(channel2,3));

size(channel1,3)
xc_size = size(xcorr2(IcT1,IcT2));
crossCorr = zeros(xc_size(1),xc_size(2),size(channel1,3));
save 1_result_crosscorr_tiff.mat crossCorr -v7.3;
save 1_result_chan1_cropped_tiff.mat chan1 -v7.3;
save 1_result_chan2_cropped_tiff.mat chan2 -v7.3;
% k = round(get(hObject,'value')*size(channel1,3))
hold(handles.axes6,'off');
for k = 1:size(channel1,3)
    
%     Imch1 = zeros(size(channel1(:,:,1)));
%     Imch2 = zeros(size(channel1(:,:,1)));

    axes(handles.axes1);
    imagesc(channel1(:,:,k));

   
    
    
    hold(handles.axes1,'off')
    axes(handles.axes2);
    imagesc(channel2(:,:,k));
    hold(handles.axes2,'off')
    axes(handles.axes1);
    % reg_edge = {[rect(2),rect(2)+rect(4)],[rect(1),rect(1)+rect(3)]};
    I1 = getimage(handles.axes1);
    Ic1 = imcrop(I1,rect);
    axes(handles.axes3);
    imagesc(Ic1)
    I2 = getimage(handles.axes2);
    Ic2 = imcrop(I2,rect);
    
    Imch1 = zeros(size(Ic1(:,:)));
    Imch1(Ic1(:,:)>thr_ch1) = 1;  %%binary image threshold channel 1
    
    Imch2 = zeros(size(Ic2(:,:)));
    Imch2(Ic2(:,:)>thr_ch2) = 1; 
    
    axes(handles.axes4);
    imagesc(Ic2)
    axes(handles.axes5);
    Icorr= xcorr2((Ic1),(Ic2));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    normCorr = normxcorr2((Ic1),(Ic2));
    imagesc(normCorr);
    example = matfile('1_result_crosscorr_tiff.mat','Writable',true);
    example.crossCorr(:,:,k) = normCorr;
    
    example_ch1 = matfile('1_result_chan1_cropped_tiff.mat','Writable',true);
    %example_ch1.chan1(:,:,k) = Ic1;
    
    example_ch2 = matfile('1_result_chan2_cropped_tiff.mat','Writable',true);
    % example_ch2.chan2(:,:,k) = Ic2;
    imwrite(Icorr, ch1_correl_ch2file, 'writemode', 'append');
    hold on
    plot(size(Icorr,2)/2,size(Icorr,1)/2,'+r')
    size(Icorr)
    plot(handles.axes6,k,normCorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2)),'*r');
    hold(handles.axes6,'on')
    round(size(Icorr,1)/2)
    round(size(Icorr,2)/2)
    Icorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2));
    data(k,1)=k;
    data(k,2)= Icorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2));
    data(k,3)= normCorr(round(size(Icorr,1)/2),round(size(Icorr,2)/2));
    pause(0.001)
    hold off
    
    
    IcT1 = IcT1+Imch1;
    imwrite(Imch1, ch1_binary, 'writemode', 'append' );
     %%binary image threshold channel 2
    IcT2 = IcT2+Imch2;
    
    imwrite(Imch2, ch2_binary, 'writemode', 'append' );
    
    ch1_or_ch2 = (Imch1+Imch2);
    ch1_or_ch2(ch1_or_ch2>1)=1;
    imwrite(ch1_or_ch2, ch1_or_ch2file, 'writemode', 'append');
    ch1_and_ch2 = Imch1.*Imch2;
    imwrite(ch1_and_ch2, ch1_and_ch2file, 'writemode', 'append');
    
    %Statistical analysis
    CC_ch1 = bwconncomp(Imch1);
    CC_ch2 = bwconncomp(Imch2);
    CC_ch1_or_ch2 = bwconncomp(ch1_or_ch2);
    CC_ch1_and_ch2 = bwconncomp(ch1_and_ch2);
    
   
    num_channel(k+1,1) = CC_ch1.NumObjects;
    num_channel(k+1,2) = CC_ch2.NumObjects;
    num_channel(k+1,3) = CC_ch1_or_ch2.NumObjects;
    num_channel(k+1,4) = CC_ch1_and_ch2.NumObjects;
    num_channel(k+1,5) = num_channel(k+1,4)/num_channel(k+1,3)*100;
    num_channel(k+1,7) = length(find(Imch1==1));
    num_channel(k+1,8) = length(find(Imch2==1));
    num_channel(k+1,9) = length(find(ch1_or_ch2==1));
    num_channel(k+1,10) = length(find(ch1_and_ch2==1));
    num_channel(k+1,11) = num_channel(k+1,10)/num_channel(k+1,9);
     
end
% axes(handles.axes1);
% imagesc(Imch1);
% hold(handles.axes1,'off')
% axes(handles.axes2);
% imagesc(Imch2);
% hold(handles.axes2,'off')
xlswrite(filename,data)
xlswrite(statistics,num_channel);

IcTMax = max(IcT1,IcT2);
maxI = max(IcTMax(:));

axes(handles.axes7);
mesh(IcTMax,'FaceAlpha',0.1)

IcTFilter = IcTMax;
IcTFilter(IcTFilter<tmin) = 0;
IcTFilter(IcTFilter>tmax) = 0;
IROI = maxI*imregionalmax(IcTFilter,26);
% IcTFilter = imregionalmax(IcTFilter);
hold(handles.axes7,'on')
imagesc(IROI);
IcMxlevel = tmax*ones(size(IcTFilter));
surf(IcMxlevel,'FaceAlpha',0.1)
hold(handles.axes7,'on')
IcMnlevel = tmin*ones(size(IcTFilter));
surf(IcMnlevel,'FaceAlpha',0.1)
hold(handles.axes7,'off')
set(h6,'UserData',IROI);


col_header={'cha1','cha2','OR','AND','Percentage','','Areacha1','Areacha2','AreachaORch2','Areach1ANDch2','PercentageArea'};
xlswrite(statistics,col_header,'Sheet1','A1');     %Write column header
% xlswrite(statistics,cha2,'Sheet1','A1');     %Write column header
% xlswrite(statistics,OR,'Sheet1','A1');     %Write column header
% xlswrite(statistics,AND,'Sheet1','A1');     %Write column header
% xlswrite(statistics,Percentage,'Sheet1','A1');     %Write column header


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
h1 = findobj('Tag','chan1');
h2 = findobj('Tag','chan2');
channel1 = get(h1,'UserData');
channel2 = get(h2,'UserData');
k = round(get(hObject,'value')*size(channel1,3))
axes(handles.axes1);
imagesc(channel1(:,:,k));
hold(handles.axes1,'off')
axes(handles.axes2);
imagesc(channel2(:,:,k));
hold(handles.axes2,'off')

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function maxlag_Callback(hObject, eventdata, handles)
% hObject    handle to maxlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxlag as text
%        str2double(get(hObject,'String')) returns contents of maxlag as a double


% --- Executes during object creation, after setting all properties.
function maxlag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
