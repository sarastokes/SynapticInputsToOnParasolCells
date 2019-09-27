function T = nearestBipolarInput(neuron, amacrineSynapseID, visualize)

    if nargin < 3
        visualize = false;
    end

    volumeScale = getODataScale('t')/1000;
    bipolars = neuron.getSynapseNodesByType('RibbonPost', false);

    endpoint = getODataURL(amacrineSynapseID, 't', 'location');
    data = readOData(endpoint);

    data = cat(1, data.value{:});

    X = volumeScale(1) * vertcat(data.VolumeX);
    Y = volumeScale(2) * vertcat(data.VolumeY);
    Z = volumeScale(3) * vertcat(data.Z);

    if numel(X) > 1
        warning('%u has %u locations!', amacrineSynapseID, numel(X));
        X = mean(X); Y = mean(Y); Z = mean(Z);
    end

    allDistances = fastEuclid3d([X, Y, Z], bipolars.XYZum);
    [distances, ind] = mink(allDistances, 5);
    T = table(bipolars{ind, 'ID'}, bipolars{ind, 'ParentID'},...
        bipolars{ind, 'XYZum'}, distances,...
        'VariableNames', {'LocationID', 'BipolarID', 'BipolarXYZ', 'Distance'}); 

    
    if visualize
        neuron.render();
        bipolarXYZ = bipolars{ind(1), 'XYZum'};

        plot3(X, Y, Z, 'Color', 'k', 'Marker', 'o', 'MarkerFaceColor', 'r');
        plot3(bipolarXYZ(1), bipolarXYZ(2), bipolarXYZ(3), 'Color', [0, 0.8, 0.3], 'Marker', 'o' , 'MarkerFaceColor', [0, 0.8, 0.3]);
        plot3([bipolarXYZ(1), X], [bipolarXYZ(2), Y], [bipolarXYZ(3), Z],... 
            'Color', 'b', 'LineWidth', 1.25);
    end