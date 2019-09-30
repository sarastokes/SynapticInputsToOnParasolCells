classdef InhibitionDistanceView < handle

    properties (SetAccess = private)
        neuron
        % AC synapses from wide-field amacrines
        wfAmacrines
        % AC synapses from other amacrines
        otherAmacrines
        % Amacrine cell synapses and associated amacrine cell IDs
        linkedAmacrines

        analysis
        bipolarSynapseID
        bipolarAnnotations
        % Single location for bipolar cell synapse (for plotting, not analysis)
        bipolarXYZ
    end

    properties (Hidden, Access = private)
        figureHandle
        axesHandle
        
        currentSynapse
    end

    properties (Hidden, Constant)
        NUM_AMACRINES = 10;
    end

    methods 
        function obj = InhibitionDistanceView(neuron, bipolarSynapseID, linkedAmacrines)
            fprintf('Initializing InhibitionDistanceView\n');
            obj.neuron = neuron;
            obj.bipolarSynapseID = bipolarSynapseID;
            obj.currentSynapse = 1;

            obj.neuron.checkSynapses();
            amacrines = neuron.getSynapseNodesByType('ConvPost', false);

            if nargin < 3
                fprintf('\tQuerying linked amacrines...\n');
                obj.linkedAmacrines = getLinkedNeurons(neuron, 'ConvPost');
            else
                obj.linkedAmacrines = linkedAmacrines;
            end

            fprintf('\tSetting up amacrine cell synapse tables\n');

            wfAmacrineList = dlmread([...
                fileparts(mfilename('fullpath')), ...
                '\data\widefield_amacrine_ids.txt']);
            wfSynapses = linkedAmacrines{...
                    ismember(linkedAmacrines.NeuronID, wfAmacrineList), 'SynapseID'};
            obj.wfAmacrines = amacrines(ismember(amacrines.ParentID, wfSynapses), :);
            obj.otherAmacrines = amacrines(~ismember(amacrines.ParentID, wfAmacrineList), :);
            

            fprintf('\tRunning initial analysis\n');
            obj.doAnalysis();
            fprintf('\tCreating user interface\n');
            obj.createUi();
        end
    end

    methods (Access = private)

        function onEditExcitatoryInput(obj, src, ~)
            try
                newID = str2double(src.String);
            catch
                src.String = 'error';
                return;
            end
            obj.bipolarSynapseID = newID;
            obj.doAnalysis();
            set(findobj(obj.axesHandle, 'Tag', 'Bipolar'),...
                'XData', obj.bipolarXYZ(1), 'YData', obj.bipolarXYZ(2), ...
                'ZData', obj.bipolarXYZ(3));
            set(findobj(obj.figureHandle, 'Tag', 'BipolarSynapseID'),...
                'String', num2str(obj.bipolarSynapseID));
            set(findobj(obj.figureHandle, 'Tag', 'NumLocations'),...
                'String', [num2str(obj.bipolarAnnotations), ' locations']);
            
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
            % amacrines = obj.neuron.getSynapseNodesByType('ConvPost', false);
            amacrines = obj.wfAmacrines;

            endpoint = getODataURL(obj.bipolarSynapseID, 't', 'location');
            data = readOData(endpoint);
            data = cat(1, data.value{:});

            X = volumeScale(1) * vertcat(data.VolumeX);
            Y = volumeScale(2) * vertcat(data.VolumeY);
            Z = volumeScale(3) * vertcat(data.Z);
            obj.bipolarAnnotations = numel(X);

            if numel(X) > 1
                fprintf('\t%u has %u locations.\n', obj.bipolarSynapseID, numel(X));
                obj.bipolarXYZ = [mean(X), mean(Y), mean(Z)];
            else
                obj.bipolarXYZ = [X, Y, Z];
            end
            
            % Get the minimum distances per BC annotation
            allDistances = []; allInd = [];
            for i = 1:numel(X)
                distances = fastEuclid3d([X(i), Y(i), Z(i)], amacrines.XYZum);
                [distances, ind] = mink(distances, obj.NUM_AMACRINES);

                allDistances = cat(1, allDistances, distances);
                allInd = cat(1, allInd, ind);
            end

            % Get the minimum distances across all BC annotations
            [distances, ind] = mink(allDistances, obj.NUM_AMACRINES);
            amacrineInd = allInd(ind);

            obj.analysis = table(...
                amacrines{amacrineInd, 'ID'}, amacrines{amacrineInd, 'ParentID'}, ...
                amacrines{amacrineInd, 'XYZum'}, distances, ...
                'VariableNames', {'LocationID', 'AmacrineSynapseID', 'AmacrineXYZ', 'Distance'});
        end
        
        function show(obj)
            % Update text displays
            set(findobj(obj.figureHandle, 'Tag', 'SynapseNumber'),...
                'String', num2str(obj.currentSynapse));
            set(findobj(obj.figureHandle, 'Tag', 'DistanceText'),...
                'String', num2str(obj.analysis{obj.currentSynapse, 'Distance'}));
            set(findobj(obj.figureHandle, 'Tag', 'LocationID'),...
                'String', num2str(obj.analysis{obj.currentSynapse, 'LocationID'}));
            
            amacrineSynapseID = obj.analysis{obj.currentSynapse, 'AmacrineSynapseID'};
            amacrineID = obj.linkedAmacrines{obj.linkedAmacrines.SynapseID == amacrineSynapseID, 'NeuronID'};
            set(findobj(obj.figureHandle, 'Tag', 'AmacrineSynapseID'), ...
                'String', num2str(amacrineSynapseID));
            set(findobj(obj.figureHandle, 'Tag', 'AmacrineID'),...
                'String', num2str(amacrineID));

            % Update data point display
            amacrineXYZ = obj.analysis{obj.currentSynapse, 'AmacrineXYZ'};
            set(findobj(obj.axesHandle, 'Tag', 'Amacrine'),...
                'XData', amacrineXYZ(1), 'YData', amacrineXYZ(2), 'ZData', amacrineXYZ(3));
            set(findobj(obj.axesHandle, 'Tag', 'Link'),...
                'XData', [obj.bipolarXYZ(1), amacrineXYZ(1)], ...
                'YData', [obj.bipolarXYZ(2), amacrineXYZ(2)], ...
                'ZData', [obj.bipolarXYZ(3), amacrineXYZ(3)]);
            axis(obj.axesHandle, 'tight');
        end


        function createUi(obj)
            LayoutManager = sbfsem.ui.LayoutManager;
            obj.figureHandle = figure(...
                'Name', 'Inhibition Distance View',...
                'Menubar', 'none',...
                'Toolbar', 'none',...
                'NumberTitle', 'off',...
                'DefaultUicontrolBackgroundColor', 'w',...
                'DefaultUicontrolFontName', 'Segoe Ui',...
                'DefaultUicontrolFontSize', 10);

            mainLayout = uix.HBoxFlex('Parent', obj.figureHandle,...
                'BackgroundColor', 'w');
            uiLayout = uix.VBox('Parent', mainLayout,...
                'BackgroundColor', 'w');
            axLayout = uipanel('Parent', mainLayout,...
                'BackgroundColor', 'w');
            obj.axesHandle = axes('Parent', axLayout);
            hold(obj.axesHandle, 'on');
            
            set(mainLayout, 'Widths', [-1, -4.5]);
            
            uicontrol(uiLayout,...
                'Style', 'text',...
                'String', num2str(obj.bipolarSynapseID),...
                'FontWeight', 'bold',...
                'Tag', 'BipolarSynapseID');
            uicontrol(uiLayout,...
                'Style', 'text',...
                'String', [num2str(obj.bipolarAnnotations), ' locations'], ...
                'Tag', 'NumLocations');
            
            navLayout = uix.HBox('Parent', uiLayout,...
                'BackgroundColor', 'w');
            uicontrol(navLayout,...
                'Style', 'push',...
                'String', '<--',...
                'Callback', @obj.lastPair);
            uicontrol(navLayout, ...
                'Style', 'text', ...
                'Tag', 'SynapseNumber', ...
                'FontWeight', 'bold', ...
                'String', num2str(obj.currentSynapse));
            uicontrol(navLayout,...
                'Style', 'push',...
                'String', '-->',...
                'Callback', @obj.nextPair);
            set(navLayout, 'Widths', [-1, -0.5, -1]);
            uix.Empty('Parent', uiLayout,...
                'BackgroundColor', 'w');
            
            LayoutManager.verticalBoxWithLabel(...
                uiLayout, 'Bipolar Synapse ID:',...
                'Style', 'edit', 'String', '',...
                'Callback', @obj.onEditExcitatoryInput);
            uix.Empty('Parent', uiLayout,...
                'BackgroundColor', 'w');
            
            LayoutManager.verticalBoxWithLabel(...
                uiLayout, 'Distance:',...
                'Style', 'text', 'String', '',...
                'FontWeight', 'bold',...
                'Tag', 'DistanceText');
            LayoutManager.verticalBoxWithLabel(...
                uiLayout, 'Amacrine ID',...
                'Style', 'text', 'String', '',...
                'FontWeight', 'bold',...
                'Tag', 'AmacrineID');
            LayoutManager.verticalBoxWithLabel(...
                uiLayout, 'Amacrine Synapse ID:',...
                'Style', 'text', 'String', '',...
                'Tag', 'AmacrineSynapseID');
            LayoutManager.verticalBoxWithLabel(...
                uiLayout, 'Location ID:',...
                'Style', 'text', 'String', '',...
                'Tag', 'LocationID');
            
            uix.Empty('Parent', uiLayout,...
                'BackgroundColor', 'w');
            
            uicontrol(uiLayout, ...
                'Style', 'push', ...
                'String', 'Rotate on', ...
                'Callback', @obj.onSelectedRotate);
            drawnow;
            
            % Set up the render display
            obj.neuron.render('ax', obj.axesHandle, ...
                'FaceColor', [0.5, 0.5, 0.5]);
            axis(obj.axesHandle, 'equal');
            lightangle(obj.axesHandle, 45, 30); 
            lightangle(obj.axesHandle, 225, 30); 
            grid(obj.axesHandle, 'on');
            axis(obj.axesHandle, 'tight'); 
            drawnow;
            plot3(obj.bipolarXYZ(1), obj.bipolarXYZ(2), obj.bipolarXYZ(3), ...
                'Parent', obj.axesHandle, ...
                'Marker', 'o', ...
                'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', hex2rgb('00cc4d'), ...
                'Tag', 'Bipolar');
            plot3(0, 0, 0, ...
                'Parent', obj.axesHandle, ...
                'Marker', 'o', ...
                'MarkerFaceColor', 'none', ...
                'LineWidth', 1, ...
                'MarkerFaceColor', hex2rgb('ff4040'),...
                'MarkerEdgeColor', hex2rgb('ff4040'),...
                'Tag', 'Amacrine');
            plot3([obj.bipolarXYZ(1), 0], [obj.bipolarXYZ(2), 0], [obj.bipolarXYZ(3), 0], ...
                'Parent', obj.axesHandle, ...
                'Color', 'b', 'LineWidth', 1.5, ...
                'Tag', 'Link');
            obj.show();

        end
    end
end