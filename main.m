%Truss Simulator: Alpha 1.0
%William Breckwoldt
function main
%% Reading the Materials file
fID = fopen('materials.txt');
reading = 1;
materialID = 0;
material = struct;
MEMBCOLOR = [];
while reading
    newLine = fgetl(fID);
    if ~isempty(strfind(newLine,'[END]'))
        reading = 0;
    elseif ~isempty(strfind(newLine,'[MATERIAL]'))
        materialID = materialID + 1;
        textData = strsplit(newLine,']');
        material(materialID).name = textData{2};
    elseif ~isempty(strfind(newLine,'[YOUNGS MODULUS]'))
        textData = strsplit(newLine,']');
        material(materialID).youngsModulus = str2num(textData{2});
    elseif ~isempty(strfind(newLine,'[YIELD STRENGTH]'))
        textData = strsplit(newLine,']');
        material(materialID).yieldStrength = str2num(textData{2});
    elseif ~isempty(strfind(newLine,'[TENSILE STRENGTH]'))
        textData = strsplit(newLine,']');
        material(materialID).tensileStrength = str2num(textData{2});
    elseif ~isempty(strfind(newLine,'[REDUCTION]'))
        textData = strsplit(newLine,']');
        material(materialID).reduction = str2num(textData{2});
    elseif ~isempty(strfind(newLine,'[ELONGATION]'))
        textData = strsplit(newLine,']');
        material(materialID).elongation = str2num(textData{2});
    elseif ~isempty(strfind(newLine,'[DENSITY]'))
        textData = strsplit(newLine,']');
        material(materialID).density = str2num(textData{2});
    elseif ~isempty(strfind(newLine,'[R]'))
        textData = strsplit(newLine,']');
        MEMBCOLOR(materialID,1) = str2num(textData{2})/255;
    elseif ~isempty(strfind(newLine,'[G]'))
        textData = strsplit(newLine,']');
        MEMBCOLOR(materialID,2) = str2num(textData{2})/255;
    elseif ~isempty(strfind(newLine,'[B]'))
        textData = strsplit(newLine,']');
        MEMBCOLOR(materialID,3) = str2num(textData{2})/255;
    end
end
fclose(fID);

%% Initializing the rest of the stuff
NAME = 'Truss Simulator Alpha';
BKRDCOLOR = 0.1 * [1 1 1];
TEXTCOLOR = 1 * [1 1 1];
GRIDCOLOR = 0.7 * [1 1 1];
NODECOLOR = 0.4 * [1 2 1];
RNODECOLOR = NODECOLOR;
% MEMBCOLOR = [0,0.4470,0.7410];
MEMBWIDTH = 2;
GHSTCOLOR = 0.5 * [1 1 1];
n1 = 9;n2 = 10;T = 11;
SIZE = [16, 16]*3;%Size of the map in small squares

leftBound = 1;%Minimum X value where a node can be placed, always an integer, always at least one
rightBound = SIZE(1) - 1;%Maximum X value where a node can be placed, always an integer, always greater than leftBound
bottomBound = 1;%Minimum Y value where a node can be placed, always an integer, always at least one
topBound = SIZE(2) - 1;%Maximum Y value where a node can be placed, always an integer, always greater than bottomBound

R = [10,15,20;...
    10,10,10];
for x = leftBound:rightBound%Over the width
    for y = bottomBound:topBound%Over the height
        n(x,y).test = 0;%The test number telling how many members attach to this point
        n(x,y).connection = [];%The (X,Y) coordinates of the nodes this node attaches to
        n(x,y).primary = [];%Boolean indicating if this node is the primary node of it's attached members
        n(x,y).shape = [];%Material of the connected member
        n(x,y).material = [];%Shape of the connected member
        n(x,y).member = [];%ID number of the member that forms the connection
    end
end

mTest = [];%whatever
m = struct;%Somehow this works for initializing m

aspectRatio = 0;
minZoom = 0;
maxZoom = 0;
up = 0;
right = 0;
type = 1;%Type of axial member
deleteKey = 'shift';
modKey = 'space';
tempNode = 0;
tempNodePos = [0,0];
ROOT = groot;
state = 'initialize';%State of the program
global running
global simming
simming = 0;
mainFigure = figure(...
    'Visible','off',...
    'Name', NAME,...
    'Menubar', 'none',...
    'NumberTitle', 'off',...
    'Units', 'normalized',...
    'Color', BKRDCOLOR,...
    'WindowScrollWheelFcn', @Zoom,...
    'WindowButtonMotionFcn', @mouseMove,...
    'ResizeFcn', @Resize,...
    'KeyPressFcn', @KeyPress,...
    'KeyReleaseFcn', @KeyRelease,...
    'CloseRequestFcn', @Close,...
    'OuterPosition', [.25 .25 .75 .75]);
set(mainFigure, 'Units', 'pixels');

mainAxes = axes(...
    'color', BKRDCOLOR,...
    'DataAspectRatioMode', 'Manual',...
    'PlotBoxAspectRatioMode', 'Manual',...
    'CameraViewAngleMode', 'Manual',...
    'CameraPositionMode', 'Manual',...
    'CameraUpVectorMode', 'Manual',...
    'DataAspectRatio', [1 1 1],...
    'PlotBoxAspectRatio', [1 1 1],...
    'Position', [0 0 1 1],...
    'ButtonDownFcn', @(~,~)MouseClick('A', 0),...
    'GridLineStyle', ':',...
    'XGrid', 'on',...
    'YGrid', 'on',...
    'XColor', GRIDCOLOR,...
    'YColor', GRIDCOLOR,...
    'TickDir', 'out');
set(mainAxes, 'Units', 'pixels');
set(mainAxes,...
    'PlotBoxAspectRatio', [SIZE,1],...
    'CameraViewAngle', [90],...
    'CameraPosition', [0 0 20],...
    'CameraTarget', [0 0 0],...
    'CameraUpVector', [0 1 0],...
    'XTick', [0:1:SIZE(1)],...
    'YTick', [0:1:SIZE(2)],...
    'XLim', [0,SIZE(1)],...
    'YLim', [0,SIZE(1)]);
Resize

set(mainFigure, 'Visible', 'on')
hold on
RScatter = scatter(mainAxes, R(1,:), R(2,:),...
    'filled',...
    'MarkerEdgeColor', RNODECOLOR*.5,...
    'MarkerFaceColor', RNODECOLOR,...
    'HitTest', 'off',...
    'LineWidth',5);
NScatter = scatter(mainAxes, -1, -1,...
    'filled',...
    'MarkerEdgeColor', NODECOLOR*.5,...
    'MarkerFaceColor', NODECOLOR,...
    'ButtonDownFcn', @(~,~)MouseClick('N', 0),...
    'LineWidth',1);
tempMember = line(-1, -1,...
    'Color', GHSTCOLOR,...
    'HitTest', 'off',...
    'LineWidth', MEMBWIDTH);
tempNode = scatter(mainAxes, 0, 0,...
    'filled',...
    'MarkerEdgeColor', NODECOLOR*.5,...
    'MarkerFaceColor', NODECOLOR,...
    'Visible', 'off',...
    'LineWidth',1);

button.simulate = uicontrol('Style','pushbutton',...
    'KeyPressFcn', @KeyPress,...
    'Position',[0,0,100,25],...
    'Callback', @(~,~)simulateAttempt,...
    'String','Simulate');
dropdown.materials = uicontrol('Style','popup',...
    'Position',[100,0,100,25],...
    'String',{});
materialString = {};
for k = 1:materialID
    materialString = {materialString{:}, material(k).name};
end
set(dropdown.materials, 'String', materialString)

dropdown.shape = uicontrol('Style','popup',...
    'Position',[200,0,100,25],...
    'String','Shape 1');

get(tempNode,'ZData')
movegui(mainFigure, 'center')
set(mainFigure, 'Visible', 'on')
running = 1;
POS = get(mainAxes, 'CameraPosition');
state = 'idle';
while running
    drawnow
    if right || up
        POS = get(mainAxes, 'CameraPosition');
        POS(1) = POS(1) + right*POS(3)*.025;
        POS(2) = POS(2) + up*POS(3)*.025;
        set(mainAxes, 'CameraPosition', POS);
        POS(3) = 0;
        set(mainAxes, 'CameraTarget', POS);
        Resize
    end
end
delete(mainFigure)
%% SUBFUNCTIONS
    function Zoom(~,evt)
        if ~strcmp(state, 'simming')
            dir = evt.VerticalScrollCount;
            pos = get(mainAxes, 'CameraPosition');%z is the distance away from the axis
            set(mainAxes, 'CameraPosition', pos + [0, 0, pos(3)*0.1]*dir);
            Resize;
        end
    end

    function Resize(~,~)
        set(mainAxes, 'Units', 'normalized');
        set(mainAxes, 'Position', [0,0,1,1]);
        set(mainAxes, 'Units', 'pixels');
        dims = get(mainAxes, 'Position');
        aspectRatio = dims(3)/dims(4);%Width/Height
        maxZoom = min((SIZE(1)/2)*aspectRatio, (SIZE(2)/2)/aspectRatio);
        minZoom = max(2*aspectRatio, 2/aspectRatio);
        pos = get(mainAxes, 'CameraPosition');
        if pos(3) > maxZoom
            pos(3) = maxZoom;
        elseif pos(3) < minZoom
            pos(3) = minZoom;
        end
        if aspectRatio >= 1
            if pos(1) - pos(3)*aspectRatio <= 0
                pos(1) = pos(3)*aspectRatio;
            end
            if pos(1) + pos(3)*aspectRatio >= SIZE(1)
                pos(1) = SIZE(1) - pos(3)*aspectRatio;
            end
            if pos(2) - pos(3) <= 0
                pos(2) = pos(3);
            end
            if pos(2) + pos(3) >= SIZE(2)
                pos(2) = SIZE(2) - pos(3);
            end
        elseif aspectRatio <= 1
            if pos(1) - pos(3) <= 0
                pos(1) = pos(3);
            end
            if pos(1) + pos(3) >= SIZE(1)
                pos(1) = SIZE(1) - pos(3);
            end
            if pos(2) - pos(3)/aspectRatio <= 0
                pos(2) = pos(3)/aspectRatio;
            end
            if pos(2) + pos(3)/aspectRatio >= SIZE(2)
                pos(2) = SIZE(2) - pos(3)/aspectRatio;
            end
        end
        set(mainAxes, 'CameraPosition', pos);
        pos(3) = 0;
        set(mainAxes, 'CameraTarget', pos);
    end

    function KeyPress(~,evt)
        if strcmp(evt.Key,'w')
            up = 1;
        elseif strcmp(evt.Key,'a')
            right = -1;
        elseif strcmp(evt.Key,'s')
            up = -1;
        elseif strcmp(evt.Key,'d')
            right = 1;
        elseif strcmp(evt.Key,deleteKey) && strcmp(state,'idle')
            disp('DELETE')
            setState('delete')
        elseif strcmp(evt.Key,modKey) && strcmp(state,'idle')
            state = 'modify';
%         elseif strcmp(evt.Key,'return') && strcmp(state,'simulate')
%             simming = 0;
%         elseif strcmp(evt.Key,'return') && strcmp(state,'idle')
%             state = 'simulate';
%             simming = 1;
%             simulate(N,M,nCheck,mCheck, MLine, NScatter, R);
%             simming = 0;
%             mm = find(mCheck ~= 0);
%             for i = 1:length(mm)
%                 set(MLine(mm(i)), 'XData',[N(M(mm(i),n1),1,1),N(M(mm(i),n2),1,1)],'YData',[N(M(mm(i),n1),1,2),N(M(mm(i),n2),1,2)], 'Visible', 'on', 'Color', MEMBCOLOR, 'LineWidth', MEMBWIDTH)
%             end
%             set(NScatter, 'XData', N(nCheck ~= 0,1,1));
%             set(NScatter, 'YData', N(nCheck ~= 0,1,2));
%             state = 'idle';
        end
    end

    function KeyRelease(~,evt)
        if strcmp(evt.Key,'w')
            up = 0;
        elseif strcmp(evt.Key,'a')
            right = 0;
        elseif strcmp(evt.Key,'s')
            up = 0;
        elseif strcmp(evt.Key,'d')
            right = 0;
        elseif strcmp(evt.Key,modKey) && strcmp(state,'modify')
            state = 'idle';
        elseif strcmp(evt.Key,deleteKey) && strcmp(state,'delete')
            state = 'idle';
        end
    end

    function Close(~, ~)
        simming = 0;
        running = 0;
    end

    function MouseClick(object, num)
        click = round(get(mainAxes, 'CurrentPoint'));%Gets click data and rounds it
        if strcmp(state, 'idle')%Click will start making a member & nodes
            if (sum(n(click(1,1),click(1,2)).test) ~= 0) || sum(click(1,1) == R(1,:) & click(1,2) == R(2,:)) %Only allows you to build on existing nodes
                set(tempNode, 'XData',click(1,1), 'YData', click(1,2), 'Visible', 'on')%Moves temp node to click position and makes temp node visible
                set(tempMember, 'XData',[click(1,1), click(2,1)],'YData',[click(1,2), click(2,2)], 'Visible', 'on')%Sets the X1, Y1, X2 and X2 of the temp member
                setState('creatingMember')%Identify waiting for second click for mouseMove (callback fcn)
            end
        elseif strcmp(state, 'creatingMember')%Click will finish creating a member & nodes
            createMember
        elseif strcmp(state, 'modify')%Click will modify a node or member
            %Changes member's material or shape
        elseif strcmp(state, 'delete')%Click will delete a node or member
            if strcmp(object, 'M')%Delete a member
                deleteMember(num)
            elseif strcmp(object, 'N')%Delete a node
                deleteNode(click(1,1),click(1,2))
            end
        end
    end

    function mouseMove(~,~)
        if strcmp(state, 'creatingMember')
            figPosRaw = get(mainFigure, 'Position');
            figPos = figPosRaw(1:2);
            POS = get(mainAxes, 'CameraPosition');
            scale = 2*POS(3)/min(figPosRaw(3:4));
            mousePos = (ROOT.PointerLocation - figPos);
            if aspectRatio >= 1
                pointerPos = mousePos*scale + POS(1:2) - [POS(3)*aspectRatio,POS(3)];
            else
                pointerPos = mousePos*scale + POS(1:2) - [POS(3),POS(3)/aspectRatio];
            end
            pointerPos = round(pointerPos);
            set(tempMember, 'XData',[get(tempNode, 'XData'), pointerPos(1)],'YData',[get(tempNode, 'YData'), pointerPos(2)])%Sets the X1, Y1, X2 and X2 of the temp member
        end
    end

    function createMember
        tempMemberXData = get(tempMember, 'XData');
        tempMemberYData = get(tempMember, 'YData');
        primaryNode = [tempMemberXData(1),tempMemberYData(1)];%Reassigning tempMember's X and Y Data into the more useable primaryNode
        secondaryNode = [tempMemberXData(2),tempMemberYData(2)];%Reassigning tempMember's X and Y Data into the more useable secondaryNode
        if ~(primaryNode(1) == secondaryNode(1) && primaryNode(2) == secondaryNode(2))%Checking if primaryNode is the same point as secondaryNode
            if (secondaryNode(1) >= leftBound && secondaryNode(1) <= rightBound) && (secondaryNode(2) >= bottomBound && secondaryNode(2) <= topBound)%Checking that secondaryNode is within the set bounds
                
                mID = find(mTest == 0,1);%Finds an avalible place to index the new member
                if isempty(mID)%if there isn't an avalible place,
                    mID = length(mTest) + 1;%adds another place to the end
                end
                
                nID = find(n(primaryNode(1),primaryNode(2)).test == 0,1);
                if isempty(nID)
                    nID = length(n(primaryNode(1),primaryNode(2)).test) + 1;
                end
                m(mID).primaryTest = nID;%sets primary node relation
                n(primaryNode(1),primaryNode(2)).test(nID) = 1;%reflects the addition of one connection to this node
                n(primaryNode(1),primaryNode(2)).connection(nID,[1,2]) = [secondaryNode(1), secondaryNode(2)];%Designates other node of connection
                n(primaryNode(1),primaryNode(2)).primary(nID) = 1;%Designates this as the primary node of this connection
                n(primaryNode(1),primaryNode(2)).member(nID) = mID;%Designates this as the primary node of this connection

                nID = find(n(secondaryNode(1),secondaryNode(2)).test == 0,1);
                if isempty(nID)
                    nID = length(n(secondaryNode(1),secondaryNode(2)).test) + 1;
                end
                m(mID).secondaryTest = nID;%sets secondary node relation
                n(secondaryNode(1),secondaryNode(2)).test(nID) = 1;%reflects the addition of one connection to this node
                n(secondaryNode(1),secondaryNode(2)).connection(nID,[1,2]) = [primaryNode(1), primaryNode(2)];%Designates other node of connection
                n(secondaryNode(1),secondaryNode(2)).primary(nID) = 0;%Designates this as the primary node of this connection
                n(secondaryNode(1),secondaryNode(2)).member(nID) = mID;%Designates this as the primary node of this connection
              
                addNode(primaryNode(1),primaryNode(2))
                addNode(secondaryNode(1),secondaryNode(2))
                
                mTest(mID) = 1;%Indicates there is a member associated with mID
                m(mID).primary = primaryNode;%sets primary node
                m(mID).secondary = secondaryNode;%sets secondary node
                m(mID).material = get(dropdown.materials, 'Value');%sets material of member
                m(mID).shape = 1;%sets shape of member
                m(mID).line = line([primaryNode(1),secondaryNode(1)],[primaryNode(2),secondaryNode(2)],...
                    'Color', MEMBCOLOR(m(mID).material,:),...
                    'Visible', 'on',...
                    'ButtonDownFcn', @(~,~)MouseClick('M', mID),...
                    'LineWidth', MEMBWIDTH(m(mID).shape));%Draws the member
                uistack(m(mID).line, 'bottom')%Puts the member behind other objects
                
                set(tempNode, 'Visible', 'off')%makes the temp node vanish
                set(tempMember, 'Visible', 'off')%Sets temp member as invisible
                setState('idle')
            else
                disp('BOUNDRY REJECTION')%node placement was out of bounds
            end
        else
            disp('MEMBER LENGTH REJECTION')%primary and second nodes were the same, so no member was created.
        end
    end

    function deleteNode(X, Y)
        for j = find(n(X,Y).test == 1)%For each active connection to the node,
            deleteMember(n(X,Y).member(j))%delete the member associated with that connection
        end
    end

    function deleteMember(mID)
        set(m(mID).line, 'Visible', 'off')%Removes member from sight
        mTest(mID) = 0;%Indicates the member is not being used
        %For the member's primary node:
        n(m(mID).primary(1),m(mID).primary(2)).test(m(mID).primaryTest) = 0;%remove the indication that this member is attached to this node
        if sum(n(m(mID).primary(1),m(mID).primary(2)).test) == 0%If there are no more indicated members attached to the node,
            removeNode(m(mID).primary(1),m(mID).primary(2))%remove the node from view
        end
        %For the member's secondary node:
        n(m(mID).secondary(1),m(mID).secondary(2)).test(m(mID).secondaryTest) = 0;%remove the indication that this member is attached to this node
        if sum(n(m(mID).secondary(1),m(mID).secondary(2)).test) == 0%If there are no more indicated members attached to the node,
            removeNode(m(mID).secondary(1),m(mID).secondary(2))%remove the node from view
        end
    end

    function addNode(X, Y)
        nodeXData = get(NScatter, 'XData');
        nodeYData = get(NScatter, 'YData');
        find(X == nodeXData & Y == nodeYData, 1)
        if isempty(find(X == nodeXData & Y == nodeYData, 1))%Only adds a node if no other nodes exist at that point
            set(NScatter, 'XData', [nodeXData, X]);
            set(NScatter, 'YData', [nodeYData, Y]);
        end
    end

    function removeNode(X, Y)
        nodeXData = get(NScatter, 'XData');
        nodeYData = get(NScatter, 'YData');
        set(NScatter, 'XData', nodeXData(~(X == nodeXData & Y == nodeYData)));
        set(NScatter, 'YData', nodeYData(~(X == nodeXData & Y == nodeYData)));
    end

    function simulateAttempt
        if simming == 0
            setState('simumating')
            set(button.simulate, 'String', 'Stop Simulation')
            nID = 0;
            nIDMatrix = [];
            for xx = leftBound:rightBound%Over the width
                for yy = bottomBound:topBound%Over the height
                    if sum(n(xx,yy).test) ~= 0
                        nID = nID + 1;
                        N(nID,1,[1,2]) = [xx,yy];
                        nIDMatrix(xx,yy) = nID;
                    end
                end
            end
            for j = find(mTest == 1)%For all members currently being used
                M(j, n1) = nIDMatrix(m(j).primary(1),m(j).primary(2));
                M(j, n2) = nIDMatrix(m(j).secondary(1),m(j).secondary(2));
                M(j, T) = m(j).material;
                mLine(j) = m(j).line;
            end
            RR = [];
            for i = 1:length(R(1,:))
                if sum(n(R(1,i),R(2,i)).test)
                    RR = [RR,nIDMatrix(R(1,i),R(2,i))];
                end
            end
            simming = 1;
            simulate(N, M, find(mTest == 1), mLine, NScatter, RR, material)
            simming =  0;
        else
            simming = 0;
            for j = find(mTest == 1)%For all members currently being used
                set(m(j).line,...
                    'Visible', 'on',...
                    'Color', MEMBCOLOR(m(j).material,:),...
                    'LineWidth', MEMBWIDTH(m(j).shape),...
                    'XData',[m(j).primary(1),m(j).secondary(1)],...
                    'YData',[m(j).primary(2),m(j).secondary(2)]);
            end
            nXData = [];
            nYData = [];
            for xx = leftBound:rightBound%Over the width
                for yy = bottomBound:topBound%Over the height
                    if sum(n(xx,yy).test) ~= 0
                        nXData = [nXData, xx];
                        nYData = [nYData, yy];
                    end
                end
            end
            set(NScatter, 'XData', nXData, 'YData', nYData)
            setState('idle')
            set(button.simulate, 'String', 'Simulate')
            axes(mainAxes)
        end
    end

    function setState(stateIn)
        state = stateIn;
    end
end

%% SIMULATION CODE

function simulate(N, M, m, mLine, nScatter, R, material)
disp(material)
disp('SIMULATION')
disp(R)
%% Beam Properties
% Area = 12700/(1000^2);
Area = 1730/(1000^2);
% I = 36.6/(1000^4)*10^6;
I = 0.912/(1000^4)*10^6;
%% Initialization
% set(gcf, 'CloseRequestFcn', @CloseSim)
global simming
n = 1:length(N(:,1,1));
% m = 1:length(M(:,1));
X = 1;Y = 2;
n1 = 9; n2 = 10;
tmax = 10;
dt = 0.0001;
fps = 20;
imaging = 1;
smax = ceil(tmax/dt);
s = 0;
c = 1;
P = 0*10*10^3;
Pax = Y;%(Y) or (X)
Pdir = -1;%(-1) or (1)
Pnode = 10;
FR = 1/fps;%Frame Rate
for j = m
    N1(j) = M(j,n1);
    N2(j) = M(j,n2);
end
S = 1;V = 2;A = 3;
N(:,[V,A],[X,Y]) = 0;
N2M = cell([max(n),1]);
N2M1 = cell([max(n),1]);
N2M2 = cell([max(n),1]);
for j = m
    N2M{N1(j)} = [N2M{N1(j)}, j];
    N2M{N2(j)} = [N2M{N2(j)}, j];
    N2M1{N1(j)} = [N2M1{N1(j)}, j];
    N2M2{N2(j)} = [N2M2{N2(j)}, j];
end
L = 1;Th = 2;DX = 3;DY = 4;F = 5;Fmax = 6;Fmin = 7;W = 8;T = 11;K = 12;Fx = 13; Fy = 14;
for j = m
    M(j,DX) = N(M(j,n2),S,X) - N(M(j,n1),S,X);
    M(j,DY) = N(M(j,n2),S,Y) - N(M(j,n1),S,Y);
    M(j,L) = sqrt(M(j,DX).^2 + M(j,DY).^2);
    % M(m,W) = M(m,L).*WpL(M(m,T))/2;
    M(j,W) = M(j,L).*material(M(j,T)).density*Area/2;
    M(j,K) = material(M(j,T)).youngsModulus.*Area;
    M(j,Fmax) = material(M(j,T)).yieldStrength*Area;
    disp(M(j,Fmax))
    M(j,Fmin) = -(pi^2).*I.*material(M(j,T)).youngsModulus./(M(j,L).^2);
    NW = zeros(length(n),1);
end
for i = n
    NW(i) = sum(M(N2M{i},W));
end
% cla
% axes('XLim', [0,10], 'YLim', [0,4], 'Position', [0,0,1,1])
tF = 0;%Time since last frame
f = 0;
%% Iterations
tic;
while simming
    s = s + 1;
    M(m,DX) = N(M(m,n2),S,X) - N(M(m,n1),S,X);
    M(m,DY) = N(M(m,n2),S,Y) - N(M(m,n1),S,Y);
    M(m,Th) = atan(M(m,DY)./M(m,DX));
    mTh = m(M(m,DX) < 0);
    M(mTh,Th) = M(mTh,Th) + pi;
    M(m,F) = M(m,K).*(sqrt(M(m,DX).^2 + M(m,DY).^2) - M(m,L));
    M(m,Fx) = M(m,F).*cos(M(m,Th));
    M(m,Fy) = M(m,F).*sin(M(m,Th));
    mFail = m(M(m,F) > M(m,Fmax) | M(m,F) < M(m,Fmin));
    if ~isempty(mFail)
        for j = mFail
            disp('FAILURE')
            for i = M(j,[n1,n2])
                N2M{i} = N2M{i}(N2M{i} ~= j);
                N2M1{i} = N2M1{i}(N2M1{i} ~= j);
                N2M2{i} = N2M2{i}(N2M2{i} ~= j);
                if isempty(N2M{i})
                    n = n(n ~= i);
                end
            end
            m = m(m ~= j);
            set(mLine(j), 'Visible', 'off')
%             set(mLine(j), 'Color', [0,0,0], 'LineWidth', 1)
        end
    end
    N(:,A,[X,Y]) = 0;
    for j = m
        N(M(j,n1),A,X) = N(M(j,n1),A,X) + M(j,Fx);
        N(M(j,n2),A,X) = N(M(j,n2),A,X) - M(j,Fx);
        N(M(j,n1),A,Y) = N(M(j,n1),A,Y) + M(j,Fy);
        N(M(j,n2),A,Y) = N(M(j,n2),A,Y) - M(j,Fy);
    end
    N(n,A,X) = N(n,A,X)/NW(i);
    N(n,A,Y) = N(n,A,Y)/NW(i) - 9.8;
%     if sum(n == Pnode)
%         N(Pnode,A,Pax) = N(Pnode,A,Pax) + Pdir*P/NW(Pnode);
%     end
    N(R,A,[X,Y]) = 0;
    N(n,V,[X,Y]) = (1-dt*c)*N(n,V,[X,Y]) + N(n,A,[X,Y])*dt;% - dt*c*N(n,V,[X,Y]).^3;
    N(n,S,[X,Y]) = N(n,S,[X,Y]) + N(n,V,[X,Y])*dt;
    if toc >= tF
        tF = toc + FR;
        f = f + 1;
        if imaging
            set(nScatter, 'XData', N(n,1,1));
            set(nScatter, 'YData', N(n,1,2));
            for j = m
                if M(j,F) > 0
                    C = abs(M(j,F)/M(j,Fmax));
                    if C <= 1 && C >= 0
                        set(mLine(j),'XData',[N(M(j,n1),S,X), N(M(j,n2),S,X)], 'YData' ,[N(M(j,n1),S,Y), N(M(j,n2),S,Y)], 'Color', [1,1-C,1-C], 'LineWidth', 4 - 3*C)
                    end
                else
                    C = abs(M(j,F)/M(j,Fmin));
                    if C <= 1 && C >= 0
                        set(mLine(j),'XData',[N(M(j,n1),S,X), N(M(j,n2),S,X)], 'YData' ,[N(M(j,n1),S,Y), N(M(j,n2),S,Y)], 'Color', [1-C,1-C,1], 'LineWidth', 4 + 3*C)
                    end
                end
            end
            disp(dt*s/toc)
            pause(0.0001)
        end
    end
end
tFinal = toc;
% if ~imaging
%     cla
%     for j = m
%         if M(j,F) > 0
%             line([N(M(j,n1),S,X), N(M(j,n2),S,X)],[N(M(j,n1),S,Y), N(M(j,n2),S,Y)], 'Color', [mod(abs(M(j,F)/M(j,Fmax)),1),0,0], 'LineWidth', 5)
%         else
%             line([N(M(j,n1),S,X), N(M(j,n2),S,X)],[N(M(j,n1),S,Y), N(M(j,n2),S,Y)], 'Color', [0,0,mod(abs(M(j,F)/M(j,Fmin)),1)], 'LineWidth', 5)
%         end
%     end
% end
fprintf('Time Elapsed: %.3fs\nTime Simulated: %.3fs\nRelative Speed: %.2f%%\nTime Step: %fs\n%.2f Average FPS\n',tFinal, tmax, 100*tmax/tFinal, dt, f/toc)

%% Subfunctions

    function CloseSim(~, ~)
        simming = 0;
        running = 0;
    end

end