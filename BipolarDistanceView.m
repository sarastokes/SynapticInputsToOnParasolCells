classdef BipolarDistanceView < handle

    properties (SetAccess = private)
        neuron
        analysis
        amacrineSynapseID
        amacrineXYZ
    end

    properties (Hidden, Access = private)
        figureHandle
        axesHandle
        
        currentSynapse
    end

    methods 
        function obj = BipolarDistanceView(neuron, amacrineSynapseID)
            obj.neuron = neuron;
            obj.amacrineSynapseID = amacrineSynapseID;

            obj.currentSynapse = 1;

            obj.doAnalysis();
            obj.createUi();
        end
    end

    methods (Access = private)
        function onKeyPress(obj, ~, evt)
            switch evt.Key 
                case 'rightarrow'
                    obj.nextPair();
                case 'leftarrow'
                    obj.lastPair();
            end
        end

        function onEditAmacrineID(obj, src, ~)
            try
                newID = str2double(src.String);
            catch
                src.String = 'error';
                return;
            end
            obj.amacrineSynapseID = newID;
            obj.doAnalysis();
            set(findobj(obj.axesHandle, 'Tag', 'Amacrine'),...
                'XData', obj.amacrineXYZ(1), 'YData', obj.amacrineXYZ(2),...
                'ZData', obj.amacrineXYZ(3));
            set(findobj(obj.figureHandle, 'Tag', 'AmacrineSynapseID'),...
                'String', num2str(obj.amacrineSynapseID));
            
            obj.currentSynapse = 1;
            obj.show();
        end

        function nextPair(obj, varargin)
            if obj.currentSynapse == height(obj.analysis)
                return;
            end
            obj.currentSynapse = obj.currentSynapse + 1;
            obj.show();
        end

        function lastPair(obj, varargin)
            if obj.currentSynapse == 1
                return;
            end
            obj.currentSynapse = obj.currentSynapse - 1;

            obj.show();
        end
        
        function onSelectedRotate(obj, src, ~)
            if strcmp(src.String, 'Rotate on')
                set(src, 'String', 'Rotate off');
                rotate3d(obj.axesHandle, 'on');
            else
                set(src, 'String', 'Rotate on');
                rotate3d(obj.axesHandle, 'off');
            end
        end
    end

    methods (Access = private)

        function doAnalysis(obj)
            volumeScale = getODataScale(obj.neuron.source)/1000;
            bipolars = obj.neuron.getSynapseNodesByType('RibbonPost', false);

            endpoint = getODataURL(obj.amacrineSynapseID, 't', 'location');
            data = readOData(endpoint);
            data = cat(1, data.value{:});

            X = volumeScale(1) * vertcat(data.VolumeX);
            Y = volumeScale(2) * vertcat(data.VolumeY);
            Z = volumeScale(3) * vertcat(data.Z);

            if numel(X) > 1
                warning('%u has %u locations!', obj.amacrineSynapseID, numel(X));
                X = mean(X); Y = mean(Y); Z = mean(Z);
            end
            obj.amacrineXYZ = [X, Y, Z];

            allDistances = fastEuclid3d([X, Y, Z], bipolars.XYZum);
            [distances, ind] = mink(allDistances, 5);
            obj.analysis = table(bipolars{ind, 'ID'}, bipolars{ind, 'ParentID'},...
                bipolars{ind, 'XYZum'}, distances,...
                'VariableNames', {'LocationID', 'BipolarID', 'BipolarXYZ', 'Distance'});
        end
        
        function show(obj)
            % Update text displays
            set(findobj(obj.figureHandle, 'Tag', 'SynapseNumber'),...
                'String', num2str(obj.currentSynapse));
            set(findobj(obj.figureHandle, 'Tag', 'DistanceText'),...
                'String', num2str(obj.analysis{obj.currentSynapse, 'Distance'}));
            set(findobj(obj.figureHandle, 'Tag', 'BipolarID'),...
                'String', num2str(obj.analysis{obj.currentSynapse, 'BipolarID'}));
            set(findobj(obj.figureHandle, 'Tag', 'LocationID'),...
                'String', num2str(obj.analysis{obj.currentSynapse, 'LocationID'}));

            
            bipolarXYZ = obj.analysis{obj.currentSynapse, 'BipolarXYZ'};
            set(findobj(obj.axesHandle, 'Tag', 'Bipolar'),...
                'XData', bipolarXYZ(1), 'YData', bipolarXYZ(2), 'ZData', bipolarXYZ(3));
            set(findobj(obj.axesHandle, 'Tag', 'Link'),...
                'XData', [obj.amacrineXYZ(1), bipolarXYZ(1)],... 
                'YData', [obj.amacrineXYZ(2), bipolarXYZ(2)],... 
                'ZData', [obj.amacrineXYZ(3), bipolarXYZ(3)]);
            axis(obj.axesHandle, 'tight');
        end


        function createUi(obj)
            LayoutManager = sbfsem.ui.LayoutManager;
            obj.figureHandle = figure(...
                'Name', 'Bipolar Distance View',...
                'Menubar', 'none',...
                'Toolbar', 'none',...
                'NumberTitle', 'off',...
                'DefaultUicontrolBackgroundColor', 'w',...
                'DefaultUicontrolFontName', 'Segoe Ui',...
                'DefaultUicontrolFontSize', 10,...
                'KeyPressFcn', @obj.onKeyPress);

            mainLayout = uix.HBox('Parent', obj.figureHandle,...
                'BackgroundColor', 'w');
            uiLayout = uix.VBox('Parent', mainLayout,...
                'BackgroundColor', 'w');
            obj.axesHandle = axes('Parent', mainLayout);
            hold(obj.axesHandle, 'on');
            
            set(mainLayout, 'Widths', [-1, -5]);
            
            uicontrol(uiLayout,...
                'Style', 'text',...
                'String', num2str(obj.amacrineSynapseID),...
                'FontWeight', 'bold',...
                'Tag', 'AmacrineSynapseID');
            uicontrol(uiLayout,...
                'Style', 'text',...
                'Tag', 'SynapseNumber',...
                'FontWeight', 'bold',...
                'String', num2str(obj.currentSynapse));
            
            navLayout = uix.HBox('Parent', uiLayout,...
                'BackgroundColor', 'w');
            uicontrol(navLayout,...
                'Style', 'push',...
                'String', '<--',...
                'Callback', @obj.lastPair);
            
            uicontrol(navLayout,...
                'Style', 'push',...
                'String', '-->',...
                'Callback', @obj.nextPair);
            uicontrol(uiLayout,...
                'Style', 'push',...
                'String', 'Rotate on',...
                'Callback', @obj.onSelectedRotate);
            uix.Empty('Parent', uiLayout,...
                'BackgroundColor', 'w');
            
            LayoutManager.verticalBoxWithLabel(...
                uiLayout, 'Amacrine ID:',...
                'Style', 'edit', 'String', '',...
                'Callback', @obj.onEditAmacrineID);
            LayoutManager.verticalBoxWithLabel(...
                uiLayout, 'Distance:',...
                'Style', 'text', 'String', '',...
                'Tag', 'DistanceText');
            LayoutManager.verticalBoxWithLabel(...
                uiLayout, 'Bipolar ID:',...
                'Style', 'text', 'String', '',...
                'Tag', 'BipolarID');
            LayoutManager.verticalBoxWithLabel(...
                uiLayout, 'Location ID:',...
                'Style', 'text', 'String', '',...
                'Tag', 'LocationID');
            
            set(uiLayout, 'Heights', [-1, -1, -1, -1, -1, -1, -1, -1, -1]);
            
            % Set up the render display
            obj.neuron.render('ax', obj.axesHandle, ...
                'FaceColor', [0.5, 0.5, 0.5]);
            axis(obj.axesHandle, 'equal');
            lightangle(obj.axesHandle, 45, 30); 
            lightangle(obj.axesHandle, 225, 30); 
            grid(obj.axesHandle, 'on');
            axis(obj.axesHandle, 'tight'); 
            drawnow;
            plot3(obj.amacrineXYZ(1), obj.amacrineXYZ(2), obj.amacrineXYZ(3), ...
                'Parent', obj.axesHandle, ...
                'Marker', 'o', ...
                'MarkerFaceColor', 'r', ...
                'MarkerEdgeColor', 'r', ...
                'Tag', 'Amacrine');
            plot3(0, 0, 0, ...
                'Parent', obj.axesHandle, ...
                'Marker', 'o', ...
                'MarkerFaceColor', 'none', ...
                'LineWidth', 1, ...
                'MarkerEdgeColor', hex2rgb('00cc4d'), ...
                'Tag', 'Bipolar');
            plot3([obj.amacrineXYZ(1), 0], [obj.amacrineXYZ(2), 0], [obj.amacrineXYZ(3), 0], ...
                'Parent', obj.axesHandle, ...
                'Color', 'b', 'LineWidth', 1.5, ...
                'Tag', 'Link');
            obj.show();

        end
    end


end