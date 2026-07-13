% Simple fallback for using GPU arrays
function x = gpuArrayTry(x)
    % Convert to gpuArrayTry if available
    try
        if gpuDeviceCount > 0
            x = gpuArray(x);
        end
    catch
        % No GPU available
    end
end
