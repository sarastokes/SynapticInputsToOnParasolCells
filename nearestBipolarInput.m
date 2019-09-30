function T = nearestBipolarInput(neuron, amacrineSynapseID, visualize)

    if nargin < 3
        visualize = false;
    end

    volumeScale = getODataScale(neuron.source)/1000;
    bipolars = neuron.getSynapseNodesByType('RibbonPost', false);

    endpoint = getODataURL(amacrineSynapseID, neuron.source, 'location');
    data = readOData(endpoint);

    data = cat(1, data.value{:});

    X = volumeScale(1) * vertcat(data.VolumeX);
    Y = volumeScale(2) * vertcat(data.VolumeY);
    Z = volumeScale(3) * vertcat(data.Z);

    if numel(X) > 1
        fprintf('\t%u has %u locations!\n', amacrineSynapseID, numel(X));
    end
    
    allDistances = []; allIndices = []; allLocations = [];
    for i = 1:numel(X)
        distances = fastEuclid3d([X(i), Y(i), Z(i)], bipolars.XYZum);
        [distances, ind] = mink(distances, 5);
        
        allDistances = cat(1, allDistances, distances);
        allIndices = cat(1, allIndices, ind);
        allLocations = cat(1, repmat(i, [5, 5]));
    end
    
    [distances, ind] = mink(allDistances, 5);
    bipolarIndex = allIndices(ind);
    acIndex = allLocations(ind);
    assignin('base', 'acIndex', acIndex);
    amacrineXYZ = [X(acIndex), Y(acIndex), Z(acIndex)];
    
    T = table(bipolars{bipolarIndex, 'ID'}, bipolars{bipolarIndex, 'ParentID'},...
        bipolars{bipolarIndex, 'XYZum'}, distances, amacrineXYZ,...
        'VariableNames', {'LocationID', 'BipolarID', 'BipolarXYZ', 'Distance', 'AmacrineXYZ'}); 

    
    if visualize
        neuron.render();
        bipolarXYZ = bipolars{bipolarIndex(1), 'XYZum'};

        plot3(X(acIndex(1)), Y(acIndex(1)), Z(acIndex(1)), 'Color', 'k', 'Marker', 'o', 'MarkerFaceColor', 'r');
        plot3(bipolarXYZ(1), bipolarXYZ(2), bipolarXYZ(3), 'Color', [0, 0.8, 0.3], 'Marker', 'o' , 'MarkerFaceColor', [0, 0.8, 0.3]);
        plot3([bipolarXYZ(1), X(acIndex(1))], [bipolarXYZ(2), Y(acIndex(1))], [bipolarXYZ(3), Z(acIndex(1))],... 
            'Color', 'b', 'LineWidth', 1.25);
    end