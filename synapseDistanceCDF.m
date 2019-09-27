function [x, y] = synapseDistanceCDF(neuron, varargin)

    ip = inputParser();
    addParameter(ip, 'Synapse', [], @ischar);
    addParameter(ip, 'Increment', 0.1, @isnumeric);
    addParameter(ip, 'Normalize', false, @islogical);
    addParameter(ip, 'Parent', [], @ishandle);
    parse(ip, varargin{:});

    somaXYZ = neuron.getSomaXYZ();
    somaXYZ = somaXYZ(1, :);

    if isempty(ip.Results.Synapse)
        somaDistances = fastEuclid3d(somaXYZ, neuron.getCellXYZ);
        rgb = [0, 0, 0];
        synapseName = 'Cell';
    else
        synapseName = ip.Results.Synapse;
        somaDistances = fastEuclid3d(somaXYZ, neuron.getSynapseXYZ(synapseName));
        syn = sbfsem.core.StructureTypes.fromStr(synapseName);
        rgb = syn.StructureColor;
        
    end

    x = 0:ip.Results.Increment:ceil(max(somaDistances));

    y = cumsum(histcounts(somaDistances, x));
    if ip.Results.Normalize
        y = y / max(abs(y));
    end


    if isempty(ip.Results.Parent)
        ax = axes('Parent', figure()); 
    else
        ax = ip.Results.Parent;
    end

    hold(ax, 'on');
    plot(ax, x, [0, y],... 
        'LineWidth', 1.5, 'Color', rgb,...
        'DisplayName', synapseName);

    title(ax, sprintf('c%u - Cumulative Synapse Distances', neuron.ID));
    xlabel(ax, 'Distance from soma (microns)');
    ylabel(ax, 'Number of Synapses');
    legend(ax, 'Location', 'northwest', 'FontSize', 10, 'EdgeColor', 'none');